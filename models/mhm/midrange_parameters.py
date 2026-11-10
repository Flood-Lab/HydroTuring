#!/usr/bin/env python3
"""An alternative parameter set for README.md's sensitivity table, written at build time.

Every parameter upstream marks as free (flag 1) in test_domain/mhm_parameter.nml
moved to the middle of upstream's own stated range. It is not a calibration and not
upstream's; it is a second documented point in parameter space, far from the first.
Usage: midrange_parameters.py input.nml output.nml
"""

import re
import sys

LINE = re.compile(r"^(\s*)(\w+)(\s*=\s*)([-+.\deE]+)(\s*,\s*)([-+.\deE]+)(\s*,\s*)([-+.\deE]+)(\s*,\s*)(\d)(\s*,\s*\d\s*)$")

out = []
for line in open(sys.argv[1]).read().splitlines():
    m = LINE.match(line)
    if m and m.group(10) == "1":
        low, high = float(m.group(4)), float(m.group(6))
        middle = 0.5 * (low + high)
        line = "".join(m.group(1, 2, 3, 4, 5, 6, 7)) + f"{middle:.6f}" + "".join(m.group(9, 10, 11))
    out.append(line)
open(sys.argv[2], "w").write("\n".join(out) + "\n")
