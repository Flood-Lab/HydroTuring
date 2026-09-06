# mass/phase-counterfactual

The snow-season sibling of `mass/warming-response`. There, demand changed
and the sign of the response was the question. Here demand does not
change: the only difference between the variants is whether the winter's
water falls as snow or as rain. So the integrated volumes must agree. The
snowpack shifts runoff from winter to spring and lets a little more or less
evaporate on the way, and the tolerance, five percent of the rain, is set
where the exact bucket's own timing effect lands; a model without a snow
module is unaffected and passes trivially.

What it catches is water that leaves a snowpack unreported: a sublimation
term, a learned sink on cold days, a phase split that does not conserve.
Those differences exist only when it snows, which is exactly the
comparison being made.

| Criterion | Asserts |
| --- | --- |
| `phase_invariance` | integrated runoff, and evaporation where reported, differ between the snowy and the rainy record by at most 0.05 of the rain |
| `non_degenerate` | runoff varies with the weather |

Must-fail: `reference_sublimating`, the bucket losing forty percent of every
snowfall.
