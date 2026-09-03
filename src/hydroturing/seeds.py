"""Seed policy.

Forcing is generated fresh at run time so that a model cannot memorise a case.
That only works if two further things are true:

1.  The seed is recorded, so any result is exactly reproducible. A benchmark
    whose failures cannot be reproduced is not usable in a dispute.
2.  A probe runs several seeds and all of them must pass, so a model cannot
    get through on a lucky draw.

The CI acceptance gate uses fixed seeds instead, derived from the probe id.
The gate must be deterministic, otherwise a flaky gate blocks a good PR.
"""

from __future__ import annotations

import hashlib
import secrets

MAX_SEED = 2**31 - 1


def gate_seeds(probe_id: str, n: int) -> list[int]:
    """Deterministic seeds for the probe acceptance gate."""
    digest = hashlib.sha256(probe_id.encode()).digest()
    base = int.from_bytes(digest[:8], "big")
    return [(base + i * 2654435761) % MAX_SEED for i in range(n)]


def eval_seeds(n: int) -> list[int]:
    """Fresh seeds for a model evaluation run. Recorded in the report."""
    return [secrets.randbelow(MAX_SEED) for _ in range(n)]
