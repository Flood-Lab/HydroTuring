#!/usr/bin/env python3
"""Backward-Euler confined aquifer with the same physical storage law."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reference_confined_aquifer.ht_adapter import main

if __name__ == "__main__":
    main(method="backward_euler")
