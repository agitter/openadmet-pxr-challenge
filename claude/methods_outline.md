# Methods Outline

A step-by-step skeleton of the pipeline. Each subsection has a short
description of its purpose and what was done; fill in prose, parameters,
and figures as needed.

---

## 1. PXR structure discovery and ensemble selection

Surveyed available PXR crystal structures (including a re-refined structure
set) and analyzed the ensemble to characterize the flexibility of the
ligand-binding pocket. Selected a set of receptor conformations to use as
docking and RBFE templates, motivated by the pocket's known plasticity.

## 2. Receptor preparation

Prepared each selected receptor structure for simulation (protonation,
missing-atom and missing-residue repair, and cleanup of crystallographic
artifacts) so that all receptors were consistent and simulation-ready.

## 3. Test-set clustering

Computed molecular fingerprints for the test compounds and a pairwise
similarity matrix, then clustered the test set into groups of structurally
similar compounds. Clustering defines the units within which relative free
energy perturbations are feasible.

## 4. Training-anchor assignment

For each test cluster, identified the most structurally similar training
(and, in principle, reference) compounds with known potency to serve as
experimental anchors, so that relative free energies can be converted to
absolute predicted potencies.

## 5. Cluster-to-template mapping

Mapped each compound cluster to the crystal-structure template best suited
to it, so that docking and RBFE for a given cluster use an appropriate
receptor conformation.

## 6. GNINA cross-docking

Docked each test and anchor compound into the selected receptor
conformations with GNINA, retaining per-pose scores (including the CNN
affinity and CNN score). Docking was run in rounds, with re-docking of
compounds that failed on the first pass.

## 7. RBFE pose and template selection

Used the docking results to select starting poses and the receptor template
for each compound, providing the bound-state geometries that seed the RBFE
transformations.

## 8. RBFE network planning

Constructed a relative free energy perturbation network within each cluster
(a minimum-spanning-tree topology), generating the set of ligand-to-ligand
transformation "edges" and their atom mappings that define the simulations
to run.

## 9. Timing pilot

Ran a small pilot batch of RBFE legs to measure per-leg wall-clock time and
resource needs, which informed the scale, batching, and scheduling of the
full production campaign.

## 10. RBFE production campaign

Executed all transformation legs (complex and solvent) as GPU jobs across a
distributed high-throughput computing pool. Each job included a fast-fail
check to shed unusable slots. Note the execution regime changed during the
campaign (checkpoint/resume versus fresh-start), which is relevant to
downstream compute accounting.

## 11. Salvage campaign for failed edges

Re-attempted the transformations that failed in the production campaign
using an alternative atom-mapping strategy, to test whether the failures
were caused by mapping choices or by the underlying starting geometries.

## 12. Results gathering and connectivity analysis

Collected the free energy estimates from all completed legs, combined the
primary and salvage results with provenance tracking, and analyzed how many
test compounds remained connected to an anchor through the surviving network
(and how many required a docking-only fallback).

## 13. Failure-mode analysis

Categorized the failed legs by error type and by leg type (complex versus
solvent) to identify the dominant causes of failure and where in the
pipeline they originated.

## 14. Convergence analysis

Examined per-leg convergence diagnostics (free energy uncertainty estimates
and adjacent-window phase-space overlap) to assess the statistical quality
of the free energy estimates and how it varied with perturbation size.

## 15. pEC50 propagation

Propagated predicted potencies to each connected test compound by summing
relative free energies along the network path from its anchor, and recorded
per-compound path features (path length, accumulated error, and worst-case
overlap along the path).

## 16. Phase 1 anchor-value analysis

After Phase 1 potencies were released, assessed whether adding Phase 1
compounds as additional anchors would reconnect or improve predictions for
Phase 2 compounds, and decided how to use Phase 1 data (calibration versus
additional anchors) accordingly.

## 17. Docking-versus-RBFE method comparison

Compared docking scores and RBFE-propagated predictions against the known
Phase 1 potencies, evaluating each signal's correlation with experiment and
whether the two methods provided independent information.

## 18. Edge-reliability threshold analysis

Swept a reliability threshold on the RBFE edges (using the convergence
diagnostics) and evaluated, against Phase 1, whether restricting to more
reliable edges improved agreement with experiment at the cost of coverage.

## 19. Score calibration and ensemble model sweep

Calibrated the docking and RBFE signals to the potency scale on the Phase 1
set and compared several ways of combining them (docking only, RBFE only,
fixed blend, and reliability-gated blends), using cross-validation with the
competition's bootstrapped error metric to select among models.

## 20. Blend-weight selection and final scoring

Selected the final ensemble's blend weight from the calibration analysis,
froze the complete model, applied it to all test compounds, and validated
the resulting submission against the competition's format requirements.

## 21. Data-coverage audit

Audited prediction coverage across all compounds and sets to confirm every
test compound received a prediction from at least one method, with no gaps
in the final submission.

## 22. Compute accounting

Parsed the execution logs of the entire campaign to quantify total GPU time
consumed, its distribution across hardware types and scheduling pools, and
the overhead attributable to job churn, providing a resource-cost account of
the physics-based approach.
