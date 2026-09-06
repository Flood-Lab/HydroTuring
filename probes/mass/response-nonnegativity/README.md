# mass/response-nonnegativity

Add rain on one day and compare the two hydrographs from that day on. The
perturbed one may never be lower. Stores drain faster when fuller, so more
water in can only mean more water out at every moment; there is no
mechanism in a catchment by which a storm lowers the flow a week later.

`mass/extreme-rain` asserts the integral of this; a model whose impulse
response has a negative lobe, an overshoot followed by a dip, satisfies the
integral and fails here. Learned hydrographs do that when they have fitted
sharp recessions.

| Criterion | Asserts |
| --- | --- |
| `response_nonnegativity` | after the added storm, runoff in the perturbed run is never below the control's by more than 0.001 of the added rain per day |
| `non_degenerate` | runoff varies with the weather |

Must-fail: `reference_overshooting`, the bucket with its runoff sharpened by
a derivative term, whose recession after an added storm undershoots the
control's.
