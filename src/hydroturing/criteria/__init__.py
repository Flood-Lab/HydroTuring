"""Criteria used by probes. Importing this module registers all of them."""

from hydroturing.criteria.base import (  # noqa: F401
    CRITERIA,
    FAIL,
    PAIRED,
    PASS,
    CriterionResult,
    criterion,
    get,
    is_paired,
    make_window,
    segments,
)
from hydroturing.criteria import (  # noqa: F401,E402
    closure,
    bounds,
    degeneracy,
    regime,
    response,
    stress,
    symmetry,
)

__all__ = [
    "CRITERIA",
    "PAIRED",
    "CriterionResult",
    "PASS",
    "FAIL",
    "criterion",
    "get",
    "is_paired",
    "make_window",
    "segments",
]
