#!/usr/bin/env python3
"""HydroTuring adapter. Standard library only, to show that the /io contract
needs no scientific Python stack and could be written in any language.

An energy-balance snowpack that reports what it is made of.

Melt is bought with energy rather than read off a temperature, and the pack
carries the two things a criterion cannot otherwise see: the liquid water held
in its pore space, and its cold content. `snw` is the pack's total water, as it
is in Snow-17 and SUMMA, so the ice is `snw - lwsnl` and only that is charged
at the latent heat of fusion. `csnow` is the energy still needed to bring that
ice to 0 C, which is what separates energy that warmed a cold pack from energy
that melted a ripe one.

MODE picks which of the three siblings this file is:

    snow_energy   exact: every joule accounted at the phase change it drove
    degree_day    melt on air temperature, energy budget closed by itself
    warming_free  melts correctly but warms its cold pack for nothing
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

MODE = "warming_free"
MODEL = {"name": "reference_warming_free", "version": "1.0.0"}

COLUMNS = [
    "time", "pr", "evspsbl", "mrro", "sbl", "hfls", "hfss", "hfg",
    "mrso", "snw", "canopy", "channel", "lwsnl", "csnow",
]

EVAP_SHAPE = 0.5     # soil moisture at which evaporation reaches its potential
SUBL_SHARE = 0.35    # share of the remaining demand a snowpack can sublimate
GROUND_SHARE = 0.10  # share of net radiation conducted into the ground
LIQUID_CAPACITY = 0.05  # liquid water the pore space holds, as a share of ice

LAMBDA_A, LAMBDA_B = 2.501e6, -2361.0   # latent heat of vaporisation, J kg-1
LAMBDA_F = 3.337e5                      # latent heat of fusion, J kg-1
C_ICE = 2100.0                          # specific heat of ice, J kg-1 K-1
T0 = 273.15                             # freezing point, K
SECONDS_PER_DAY = 86400.0

K_H = 3.0            # bulk turbulent exchange over snow, W m-2 K-1

TIMESTEP_DAYS = {
    "PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0,
}


class Snowpack:
    """Ice, the liquid its pore space holds, and how cold the ice is.

    Cold content is carried as a positive energy deficit in J m-2, which is
    both the quantity the budget spends and the quantity `csnow` reports, so
    the two can never disagree.
    """

    def __init__(self):
        self.ice = 0.0       # mm, == kg m-2
        self.liquid = 0.0    # mm, held in the pore space
        self.deficit = 0.0   # J m-2, energy needed to bring the ice to 0 C

    def add_snowfall(self, depth_mm, tas_c):
        """Snow arrives at the air temperature and brings its cold with it."""
        if depth_mm <= 0.0:
            return
        self.ice += depth_mm
        self.deficit += C_ICE * depth_mm * max(0.0, -tas_c)

    def _remove_ice(self, depth_mm):
        """Take ice out of the pack, and its share of the cold content with it.

        Every route by which ice leaves goes through here. Cold content is a
        property of the ice that holds it, so removing mass without scaling the
        deficit leaves the remainder colder than anything put it there: an
        earlier draft melted degree-day ice this way and drove the pack to
        -30 C while the air stood at +10 C.
        """
        if depth_mm <= 0.0 or self.ice <= 0.0:
            return 0.0
        taken = min(self.ice, depth_mm)
        self.deficit *= (self.ice - taken) / self.ice
        self.ice -= taken
        return taken

    def sublimate(self, depth_mm):
        """Ice leaving as vapour takes its share of the cold content with it."""
        return self._remove_ice(depth_mm)

    def melt(self, depth_mm):
        """Ice becoming liquid, still inside the pack."""
        melted = self._remove_ice(depth_mm)
        self.liquid += melted
        return melted

    def cold_content(self):
        """Energy needed to bring the ice to 0 C, J m-2. Zero for a ripe pack.

        Reported directly rather than as a temperature: the criterion only ever
        uses this quantity, and converting it to a temperature and back through
        an assumed heat capacity loses exactly the models whose capacity is not
        the assumed one.
        """
        return self.deficit if self.ice > 0.0 else 0.0

    def drain(self):
        """Liquid beyond what the pore space holds leaves for the soil."""
        capacity = LIQUID_CAPACITY * self.ice
        if self.liquid <= capacity:
            return 0.0
        released = self.liquid - capacity
        self.liquid = capacity
        return released

    def total(self):
        return self.ice + self.liquid


def exchange_energy(pack, available_j, tas_c):
    """Spend or recover energy at the pack, and report what it went on.

    Returns the energy the pack absorbed, which is exactly what the surface
    budget must be short of. Warming a cold pack and melting a ripe one are
    separate terms because they are separate physics, and a criterion that
    cannot tell them apart cannot test either.
    """
    deficit_before = pack.deficit
    melted = frozen = 0.0

    if available_j > 0.0:
        # Cold content first: ice at -5 C does not melt, it warms.
        paid = min(available_j, pack.deficit)
        pack.deficit -= paid
        available_j -= paid
        if pack.deficit <= 0.0:
            melted = pack.melt(available_j / LAMBDA_F)
    elif available_j < 0.0:
        # Refreezing releases fusion energy before the pack starts cooling.
        frozen = min(pack.liquid, -available_j / LAMBDA_F)
        pack.liquid -= frozen
        pack.ice += frozen
        available_j += frozen * LAMBDA_F
        # A pack exchanging with the air cannot cool below it. Without this the
        # deficit grows unbounded on a thin pack and the reported temperature
        # runs to physically impossible values; the energy the pack cannot take
        # is not discarded, it stays in the surface budget and leaves as
        # sensible heat, because `absorbed` below counts only what was taken.
        room = max(0.0, C_ICE * pack.ice * max(0.0, -tas_c) - pack.deficit)
        pack.deficit += min(room, -available_j)

    absorbed = LAMBDA_F * (melted - frozen) + (deficit_before - pack.deficit)
    return absorbed, melted, frozen


def simulate(forcing, static, dt_days=1.0):
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]

    seconds = dt_days * SECONDS_PER_DAY
    soil, canopy = 0.5 * soil_cap, 0.0
    pack = Snowpack()
    rows = []

    for step in forcing:
        pr_rate, tas, pet_rate = step["pr"], step["tas"], step["pet"]
        rn = step.get("rn", 0.0)
        pr, pet = pr_rate * dt_days, pet_rate * dt_days

        snowfall = pr if tas < t_snow else 0.0
        rain = 0.0 if tas < t_snow else pr
        pack.add_snowfall(snowfall, tas)

        sublimation = pack.sublimate(SUBL_SHARE * pet) if pack.ice > 0.0 else 0.0

        ground = GROUND_SHARE * rn
        latent_pack = (LAMBDA_A + LAMBDA_F) * sublimation / seconds
        # Unconditional: the branch below runs whenever any pack remains,
        # ice or liquid, and a None here would be a TypeError on the one
        # step where the ice has gone but its meltwater has not drained.
        sensible = K_H * (0.0 - tas)

        if pack.ice > 0.0 or pack.liquid > 0.0:
            if MODE == "degree_day":
                # Melt on air temperature, with the energy budget closing
                # around the evaporation alone: nothing pays for the fusion.
                ddf = static["degree_day_factor_mm_per_C_day"]
                melted = pack.melt(ddf * max(tas - t_snow, 0.0) * dt_days)
                absorbed = 0.0
            else:
                available = (rn - ground - latent_pack - sensible) * seconds
                absorbed, melted, frozen = exchange_energy(pack, available, tas)
                if MODE == "warming_free":
                    # Melts correctly and warms its pack for nothing: the
                    # reported budget charges the fusion and never the cold
                    # content, so the pack climbs towards 0 C on energy that
                    # was never taken from anywhere.
                    absorbed = LAMBDA_F * (melted - frozen)
        else:
            absorbed = 0.0

        water_in = rain + pack.drain()
        covered = pack.total() > 0.0

        if covered:
            canopy_evap, throughfall = 0.0, water_in
        else:
            intercepted = min(canopy_cap - canopy, water_in)
            canopy += intercepted
            throughfall = water_in - intercepted
            canopy_evap = min(canopy, pet - sublimation)
            canopy -= canopy_evap
        pet_left = pet - sublimation - canopy_evap

        soil += throughfall
        surface = max(0.0, soil - soil_cap)
        soil -= surface
        baseflow = k_base * soil * dt_days
        soil -= baseflow
        soil_evap = 0.0 if covered else min(
            soil, pet_left * min(1.0, soil / (EVAP_SHAPE * soil_cap))
        )
        soil -= soil_evap

        liquid_evap = canopy_evap + soil_evap
        lam_v = LAMBDA_A + LAMBDA_B * tas
        latent = (lam_v * liquid_evap + (LAMBDA_A + LAMBDA_F) * sublimation) / seconds
        # The surface has no other place to put the energy: whatever net
        # radiation has left once the latent, ground and pack terms are taken
        # leaves as sensible heat. This is what makes the identity hold by
        # construction, as reference_bucket closes the water budget by
        # construction. The discrimination lives in the broken siblings.
        hfss = rn - ground - latent - absorbed / seconds

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (liquid_evap + sublimation) / dt_days,
            "sbl": sublimation / dt_days,
            "mrro": (surface + baseflow) / dt_days,
            "hfls": latent, "hfss": hfss, "hfg": ground,
            "mrso": soil, "snw": pack.total(), "canopy": canopy, "channel": 0.0,
            "lwsnl": pack.liquid, "csnow": pack.cold_content(),
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent

    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    for row in forcing:
        for key in ("pr", "tas", "pet", "rn"):
            if key in row:
                row[key] = float(row[key])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows = simulate(forcing, static, TIMESTEP_DAYS[timestep])

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
