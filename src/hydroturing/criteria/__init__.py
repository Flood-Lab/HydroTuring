"""Criteria used by probes. Importing this module registers all of them."""

from hydroturing.criteria.base import (  # noqa: F401
    CRITERIA,
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    get,
    make_window,
)
from hydroturing.criteria import closure, bounds, degeneracy  # noqa: F401,E402

__all__ = ["CRITERIA", "CriterionResult", "PASS", "FAIL", "criterion", "get", "make_window"]
