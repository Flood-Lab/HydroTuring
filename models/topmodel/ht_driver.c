/* A second way into TOPMODEL's own core, for HydroTuring.
 *
 * The package's c_topmodel() runs the model and hands back Q, qo, qs, S_mean, fex and Ea. It
 * keeps the root-zone deficit (Srz), the unsaturated store (Suz) and the water still in the
 * delay function to itself and frees them, and the budget cannot be closed without them. This
 * file is the same loop as c_topmodel() over the same unmodified functions (topidx_calc,
 * memory_allocation, param_init, run_topmodel, memory_free), and it copies those stores out
 * before they are freed. It computes nothing the model does not, and changes no equation.
 * check_driver.R fails the image build if it ever differs from topmodel::topmodel().
 *
 * It is compiled against the pinned package sources, with everything in metres per step.
 */
#include "topmodel.h"

void ht_topmodel(double *parameters, double *topidx, double *delay, double *rain, double *ETp,
                 int *nidxclass, int *ntimestep, int *nch,
                 double *Q, double *qo, double *qs, double *S_mean, double *fex, double *Ea,
                 double *Srz_mean, double *Suz_mean, double *ex_mean,
                 double *transit0, double *weight_sum)
{
    int i, j, n = *nidxclass, nt = *ntimestep;
    double w;

    topidx_calc(topidx, n);
    memory_allocation(*nch, nt, n);
    param_init(parameters, delay, *nch, 0, n, nt);

    /* param_init seeds Qt with the flow already on its way to the outlet when the run starts.
       That water is the delay function's initial store. */
    *transit0 = 0.0;
    for (i = 0; i < nt; i++) *transit0 += misc.Qt[i];

    for (i = 0; i < nt; i++) run_topmodel(rain, ETp, n, i, nt);

    /* The class weights are the core's own Aatb_r (core_topmodel.c), the same ones it weights
       qv and Ea with, so these means are the stores those fluxes moved. */
    *weight_sum = 0.0;
    for (j = 0; j < n; j++)
        *weight_sum += (idxstats.Aatb_r[j] + (j < n - 1 ? idxstats.Aatb_r[j + 1] : 0.0)) / 2.0;

    for (i = 0; i < nt; i++) {
        Q[i] = misc.Qt[i];
        qo[i] = misc.qo[i][n];
        qs[i] = misc.qs[i];
        S_mean[i] = misc.S_mean[i];
        fex[i] = misc.fex[i];
        Ea[i] = misc.Ea[i][n];
        Srz_mean[i] = 0.0;
        Suz_mean[i] = 0.0;
        ex_mean[i] = misc.ex[i][n];
        for (j = 0; j < n; j++) {
            w = (idxstats.Aatb_r[j] + (j < n - 1 ? idxstats.Aatb_r[j + 1] : 0.0)) / 2.0;
            Srz_mean[i] += w * misc.Srz[i][j];
            Suz_mean[i] += w * misc.Suz[i][j];
        }
    }

    memory_free(*nch, nt, n);
}
