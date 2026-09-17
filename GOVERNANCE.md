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

**Reviewers.** Contributing authors who review probes for one or more
conservation laws. They check a probe's physics and its implementation, and
approve it or request changes on the pull request. The current pool is listed
under [Reviewer pool](#reviewer-pool).

**Probe authors.** Own the physics of the probes they contribute and are
consulted before anyone changes their tolerance.

**Model proposers.** Credited for the proposal whether or not they package
the model, because deciding what belongs in the benchmark is a separate
judgement from wrapping it in a container.

## Decisions

**Accepting a proposal.** The maintainer labels a probe or model proposal
`accepted` and assigns it. For a probe the question is whether it
discriminates; for a model, whether the benchmark would learn something from
its verdict and whether it can be containerised. Acceptance is deliberately
cheap and early, because the point of proposing first is to fail fast on
paper rather than slowly in code. The issue is assigned to whoever will do the
work, which for a model is often the maintainer rather than the proposer.

**Merging a probe.** Every probe pull request is assigned two reviewers from
the pool for its conservation law. Between them, the two reviews cover the
physics and the implementation. The probe is merged only after both reviewers
approve it and the acceptance gate is green. The maintainer merges.

This rule applies to every probe merged from 16 September 2026. Probes merged
before then were accepted by the maintainer under the earlier process.

**Merging a model.** One review, on the contract rather than the physics: the
adapter honours `/io`, `emits` is honest, and the image is reproducible. The
verdict itself is never grounds for rejection. A model that fails is a result,
not a defect.

**Changing a tolerance on a merged probe.** Requires the probe's authors and
the steering committee, because it silently changes every recorded result.
Any such change bumps the suite minor version.

**Cutting a suite version.** The steering committee. Scores are only
comparable within a suite version, so this is a real decision rather than a
release chore.

**Removing a probe.** Only if it is shown not to discriminate, which normally
means the acceptance gate stops separating the reference models. Recorded
results referencing it stay in the history.

## Reviewer pool

Reviewers are drawn from the project's contributing authors. The maintainer
invites them and keeps this table current. If you have a merged probe and
would like to review, say which conservation law in an issue.

| Conservation law | Reviewers |
| --- | --- |
| Energy | Changming Li (SCUT, @licm13), Han Wang (The Hong Kong University of Science and Technology, @cehw), Xin Lan (Michigan State University, @xinlan-technology) |
| Mass | Zhi Li (CU Boulder, @chrimerss), Siddik Barbhuiya (IIT Mandi, @Barbhuiya12) |
| Momentum | Zhi Li (CU Boulder, @chrimerss), Siddik Barbhuiya (IIT Mandi, @Barbhuiya12) |

## Conflicts of interest

Nobody reviews their own probe. If a law's pool has fewer than two reviewers
free of a conflict, the maintainer assigns reviewers from another pool.

A reviewer whose own model would be evaluated by a probe declares it in the
pull request before reviewing. The maintainer then decides whether they stay
on as one of the two reviewers or are replaced under the rule above. Nobody
may be the sole reviewer of a model they contributed. Say so in the pull
request and another reviewer will be assigned. Conflicts like these will be
common and are not a problem as long as they are declared.

## Changing this document

Open a pull request. Substantive changes need steering committee agreement.
