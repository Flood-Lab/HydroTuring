#!/usr/bin/env python3
"""Export a neuralhydrology training scaler from its pickle to JSON, at build time.

The pickle holds xarray Datasets from 2020 that no longer open on a current stack,
as the Zenodo record's own 2026 note says of its pickles. Run under the versions the
paper's environment.yml names, this writes the same numbers in a form that does not
care. Run only inside the image build, on the archive pinned in the Dockerfile.

    export_scaler.py train_data_scaler.p scaler.json
"""

import json
import pickle
import sys


def main(source, target):
    with open(source, "rb") as fh:
        scaler = pickle.load(fh)
    # This revision's name for the attribute statistics, and the older one.
    means = scaler["attribute_means"] if "attribute_means" in scaler else scaler["camels_attr_means"]
    stds = scaler["attribute_stds"] if "attribute_stds" in scaler else scaler["camels_attr_stds"]
    out = {
        "xarray_means": {str(k): float(v) for k, v in scaler["xarray_means"].data_vars.items()},
        "xarray_stds": {str(k): float(v) for k, v in scaler["xarray_stds"].data_vars.items()},
        # A list of pairs would lose nothing either, but JSON objects keep their order
        # here and that order is the one the network was trained on.
        "attribute_means": {str(k): float(v) for k, v in means.items()},
        "attribute_stds": {str(k): float(v) for k, v in stds.items()},
    }
    with open(target, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
