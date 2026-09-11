#!/usr/bin/env python3
"""Keep the conductive reference fluxes but freeze its reported temperature."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reference_soil_heat.ht_adapter import main, simulate as reference_simulate

MODEL = {"name": "reference_frozen_soil", "version": "1.0.0"}


def simulate(forcing, static, dt_days=1.0 / 24.0):
    return reference_simulate(forcing, static, dt_days, temperature_scale=0.0)


if __name__ == "__main__":
    sys.exit(main(model=MODEL, simulate_model=simulate))
