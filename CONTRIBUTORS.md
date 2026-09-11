# Contributors

Everyone who has shaped HydroTuring. Probe authors are also recorded in each
probe's `probe.yaml` and appear in every report that runs their probe.

See [CONTRIBUTING.md](CONTRIBUTING.md#credit) for how credit and authorship
work here.

The counts on the repository page (probes merged, models evaluated) are
generated from `probes/` and `models/` by `scripts/badges.py` every time
the site deploys, so this file and those badges cannot drift apart for
long.

## Probes

| Probe | Authors |
| --- | --- |
| `mass/catchment-closure` | Zhi Li (CU Boulder) |
| `mass/resolution-invariance` | Zhi Li (CU Boulder) |
| `mass/warming-response` | Zhi Li (CU Boulder) |
| `mass/causality` | Zhi Li (CU Boulder) |
| `mass/dry-down` | Zhi Li (CU Boulder) |
| `mass/steady-state` | Zhi Li (CU Boulder) |
| `mass/extreme-rain` | Zhi Li (CU Boulder) |
| `mass/runoff-bounds` | Zhi Li (CU Boulder) |
| `mass/area-invariance` | Zhi Li (CU Boulder) |
| `mass/response-nonnegativity` | Zhi Li (CU Boulder) |
| `mass/antecedent-monotonicity` | Zhi Li (CU Boulder) |
| `mass/phase-counterfactual` | Zhi Li (CU Boulder) |
| `mass/time-origin-invariance` | Siavash Shams (Columbia University) |
| `mass/precipitation-counterfactual` | Qingyi Yang (Politecnico di Milano) |
| `energy/pet-consistency` | Zhi Li (CU Boulder) |
| `energy/latent-heat-et-consistency` | Changming Li (SCUT) |
| `energy/evaporative-partition` | Changming Li (SCUT) |
| `energy/surface-energy-closure` | Han Wang (The Hong Kong University of Science and Technology) |
| `momentum/routing-conservation` | Zhi Li (CU Boulder) |

## Models

Proposed models, with who proposed them and who packaged them. These are
separate contributions and are recorded separately: deciding what belongs in
the benchmark is a different judgement from wrapping it in a container, and
five accepted proposals count as one probe towards authorship. Every
evaluation is archived in [models/result.csv](models/result.csv).

| Model | Proposed by | Packaged by |
| --- | --- | --- |
| `google_flood_forecast` | Zhi Li (CU Boulder), [#1](../../issues/1) | HydroTuring maintainers |
| `dhbv2` | Zhi Li (CU Boulder), [#2](../../issues/2) | HydroTuring maintainers |
| `wflow_sbm` | [@kawh1111](https://github.com/kawh1111), [#29](../../issues/29) | HydroTuring maintainers |
| `cwatm` | [@kawh1111](https://github.com/kawh1111), [#28](../../issues/28) | HydroTuring maintainers |

Physical reference models, which every compatible probe must pass when their
outputs support its criteria:

| Model | Author | Source |
| --- | --- | --- |
| `flex_lumped` | Zhi Li (CU Boulder) | [chrimerss/HydrologicModels](https://github.com/chrimerss/HydrologicModels), `lumped_model/` |
| `flex_topo` | Zhi Li (CU Boulder) | [chrimerss/HydrologicModels](https://github.com/chrimerss/HydrologicModels), `semi-distributed_model/` |
| `sacsma_snow17` | E. Anderson and NWS/HRL (model); Upstream Tech (packaging); HydroTuring maintainers (port) | [Upstream-Tech/SACSMA-SNOW17](https://github.com/Upstream-Tech/SACSMA-SNOW17) |

## Harness and infrastructure

- Zhi Li (CU Boulder), maintainer

## Scientific steering committee

Being formed. See [GOVERNANCE.md](GOVERNANCE.md). If you would be willing to
serve, or want to suggest someone, open an issue.

## Reviewers

Probe review is real work and is recorded here. Reviewers of merged probes
will be listed as they accumulate.
