#!/usr/bin/env python3
"""Conservative negative control: amplifies the departure of river stage from initial head by 40%."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reference_confined_aquifer.ht_adapter import main

if __name__ == "__main__":
    main(head_gain=1.4)
