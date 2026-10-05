# mass/groundwater-datum-invariance

A common translation of the numerical vertical datum must not alter the
river–aquifer exchange trajectory. `control`, `datum_up` and `datum_down`
share the same stage differences, storage properties and conductance; only the
absolute stage, initial head and elevations are shifted.

The paired `datum_flux_invariance` criterion compares the non-cancelling,
absolute exchange difference against the control gross exchange, with a small
absolute floor. It also requires at least 0.1 mm of activity. Groundwater
balance and signed component checks run independently for every variant.
The criterion compares the exchange only: each variant's `gw` is held to its
own balance and bounds, not to the control's.

`reference_datum_exact` is the only `must_pass` baseline, under the exception
in `docs/writing-a-probe.md` for a process none of the four physical models
has: none of them exchanges water between an aquifer and a river, so none
has a head or an elevation to shift. `models/modflow6`, a real MODFLOW 6.7.0
binary run through Docker, is the independent check. It is archived as a
submitted model rather than gated on, and passes with a worst residual of
about 2e-12 mm. Re-run it whenever the criterion changes:
`ht run --model modflow6 --probe mass/groundwater-datum-invariance --gate-seeds --csv models/result.csv`.

The allowance, 1e-6 of the control's gross exchange or 1e-5 mm, is set for
double precision. The same exact physics in single precision fails it,
because rounding at a 350 m datum is coarser than at 100 m: by 1.7e-4 to
6.1e-4 mm when the exchange is computed as conductance times the head
difference, and by 1.2e-2 to 2.7e-2 mm when it is taken from the change in
stored head. Whether the allowance should absorb single-precision rounding
is still open, so read a failure from a float32 model with that in mind.
