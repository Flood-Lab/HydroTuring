/* HydroTuring driver for NOAA-OWP CFE, linked against the unmodified libcfebmi.
 *
 * ht_cfe_driver <config> <step_seconds> <forcing.bin> <output.bin>
 *
 * forcing.bin holds two doubles per row: the rain depth for the row in mm (BMI's
 * "mm h-1" is a depth per one-hour row, which is what Update divides by 1000) and
 * the potential evaporation in m/s. output.bin holds N_OUT doubles per row, in the
 * order the adapter expects, with the state before the first row in front. Upstream's
 * printf diagnostics and mass_balance_check go to stdout, which the adapter keeps.
 *
 * BMI exposes the soil and groundwater stores and the outgoing fluxes. It does not
 * expose the water held in either Nash cascade, nor the ET components and mass-balance
 * terms, so those are read from the model's own struct through upstream's header.
 * Nothing here computes hydrology.
 *
 * time_step_size is the one setting BMI has no config key for: new_bmi_cfe sets it
 * to 3600. It is set here before Initialize, as Update_until does for a fractional
 * step, because the harness's resolution probe runs every model at the next finer
 * step whatever the manifest lists.
 */
#include <stdio.h>
#include <stdlib.h>
#include "bmi.h"
#include "bmi_cfe.h"

#define N_OUT 19

static double value(Bmi *model, const char *name)
{
    double out = 0.0;
    if (model->get_value(model, name, &out) != BMI_SUCCESS) {
        fprintf(stderr, "get_value failed for %s\n", name);
        exit(3);
    }
    return out;
}

static void set(Bmi *model, const char *name, double in)
{
    if (model->set_value(model, name, &in) != BMI_SUCCESS) {
        fprintf(stderr, "set_value failed for %s\n", name);
        exit(3);
    }
}

int main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr, "usage: %s config step_seconds forcing.bin output.bin\n", argv[0]);
        return 2;
    }

    Bmi *model = (Bmi *)malloc(sizeof(Bmi));
    register_bmi_cfe(model);
    cfe_state_struct *cfe = (cfe_state_struct *)model->data;
    cfe->time_step_size = atoi(argv[2]);
    if (model->initialize(model, argv[1]) != BMI_SUCCESS) {
        fprintf(stderr, "initialize failed for %s\n", argv[1]);
        return 3;
    }

    FILE *in = fopen(argv[3], "rb");
    FILE *out = fopen(argv[4], "wb");
    if (!in || !out) {
        fprintf(stderr, "cannot open forcing or output file\n");
        return 3;
    }

    /* The state before the first row: the same columns, with no fluxes. */
    double initial[N_OUT] = {0.0};
    initial[6] = value(model, "SOIL_STORAGE");
    initial[7] = value(model, "GW_STORAGE");
    for (int k = 0; k < cfe->nash_surface_params.N_nash; k++)
        initial[8] += cfe->nash_surface_params.nash_storage[k];
    for (int k = 0; k < cfe->N_nash_subsurface; k++)
        initial[9] += cfe->nash_storage_subsurface[k];
    if (fwrite(initial, sizeof(double), N_OUT, out) != N_OUT) return 3;

    double row[2];
    while (fread(row, sizeof(double), 2, in) == 2) {
        set(model, "atmosphere_water__liquid_equivalent_precipitation_rate", row[0]);
        set(model, "water_potential_evaporation_flux", row[1]);
        model->update(model);

        double nash_surface = 0.0, nash_subsurface = 0.0;
        for (int k = 0; k < cfe->nash_surface_params.N_nash; k++)
            nash_surface += cfe->nash_surface_params.nash_storage[k];
        for (int k = 0; k < cfe->N_nash_subsurface; k++)
            nash_subsurface += cfe->nash_storage_subsurface[k];

        double record[N_OUT] = {
            value(model, "Q_OUT"),                                        /*  0 */
            value(model, "ACTUAL_ET"),                                    /*  1 */
            value(model, "DIRECT_RUNOFF"),                                /*  2 */
            value(model, "NASH_LATERAL_RUNOFF"),                          /*  3 */
            value(model, "DEEP_GW_TO_CHANNEL_FLUX"),                      /*  4 */
            value(model, "INFILTRATION_EXCESS"),                          /*  5 */
            value(model, "SOIL_STORAGE"),                                 /*  6 */
            value(model, "GW_STORAGE"),                                   /*  7 */
            nash_surface,                                                 /*  8 */
            nash_subsurface,                                              /*  9 */
            cfe->et_struct.actual_et_from_rain_m_per_timestep,            /* 10 */
            cfe->et_struct.actual_et_from_retention_depth_m_per_timestep, /* 11 */
            cfe->et_struct.actual_et_from_soil_m_per_timestep,            /* 12 */
            cfe->infiltration_depth_m,                                    /* 13 */
            cfe->nash_surface_params.runon_infiltration,                  /* 14 */
            *cfe->flux_perc_m,                                            /* 15 */
            *cfe->flux_lat_m,                                             /* 16 */
            cfe->vol_struct.volin,                                        /* 17 */
            cfe->vol_struct.volout,                                       /* 18 */
        };
        if (fwrite(record, sizeof(double), N_OUT, out) != N_OUT) return 3;
    }
    fclose(in);
    fclose(out);

    /* Upstream's own closing report, kept in run.json. */
    mass_balance_check(cfe);
    fflush(stdout);
    return 0;
}
