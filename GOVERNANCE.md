# Governance

## Why this document exists

Outsiders write the physics tests here, and reasonable hydrologists disagree
about tolerances, denominators, and what counts as a violation. That
disagreement is the scientific content of the project, not a nuisance, but it
still needs a way to be settled.

## Roles

**Maintainer.** Zhi Li (University of Colorado Boulder). Merges pull
requests, cuts suite releases, and is responsible for the repository.

**Scientific steering committee.** Four to six people spanning hydrologic
modelling, machine learning for Earth science, and operational hydrology.
They settle disputes about probe design, approve changes to a merged probe's
tolerance, and decide when a suite version is cut.

> The committee is being formed. If you would be willing to serve, or want to
> suggest someone, open an issue or email the maintainer.

**Probe authors.** Own the physics of the probes they contribute and are
consulted before anyone changes their tolerance.

## Decisions

**Merging a probe.** Two reviews, one on the physics and one on the
implementation, plus a green acceptance gate. The maintainer merges.

**Changing a tolerance on a merged probe.** Requires the probe's authors and
the steering committee, because it silently changes every recorded result.
Any such change bumps the suite minor version.

**Cutting a suite version.** The steering committee. Scores are only
comparable within a suite version, so this is a real decision rather than a
release chore.

**Removing a probe.** Only if it is shown not to discriminate, which normally
means the acceptance gate stops separating the reference models. Recorded
results referencing it stay in the history.

## Conflicts of interest

A person may not be the sole reviewer of a probe that their own model would
be evaluated against, nor of a model they contributed. Say so in the pull
request and a second reviewer will be assigned. This will be common and is
not a problem as long as it is declared.

## Changing this document

Open a pull request. Substantive changes need steering committee agreement.
