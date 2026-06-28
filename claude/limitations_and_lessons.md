# RBFE Pipeline — Limitations & Lessons Learned (outline)

Notes for the Phase 2 writeup. Two themes: (A) network design choices that
didn't exploit the two-phase challenge structure, and (B) protocol/input
choices that drove the ~47% complex-leg failure rate.

---

## A. Network design limitations

### A1. Anchor selection did not anticipate Phase 1 unblinding
- MSTs were built with TRAINING compounds as the only anchors, fixed before
  Phase 1 was unblinded.
- The challenge guaranteed Phase 1 EC50s would be released before the Phase 2
  deadline, so Phase 1 compounds were always going to become high-quality
  experimental anchors — but the topology was not designed to route Phase 2
  compounds toward them.
- Phase 1 and Phase 2 compounds were treated identically as "test" nodes.
- Empirical cost: adding Phase 1 compounds as anchors post hoc reconnected
  0 additional Phase 2 compounds (only 4 got shorter paths). The benefit was
  unrealizable after the fact precisely because the topology wasn't built for it.

### A2. MST topology is maximally fragile (BIGGEST issue)
- MST uses exactly N-1 edges per cluster — every edge is a bridge.
- A single failed edge disconnects every compound downstream of it.
- Combined with the ~47% complex-leg failure rate, this fragility disconnected
  43% of test compounds (221/513) from any anchor.
- Redundant topologies (cycles, LOMAP-style networks with extra edges, or
  explicit k-edge-connectivity) would have provided alternate paths around
  failed edges.
- This compounds with A1: redundant edges to Phase 1 neighbors would have
  created exactly the alternate anchor paths the post-hoc analysis found missing.

### A3. Single anchor per cluster
- Each cluster propagated from essentially one training anchor.
- Multiple anchors (training + eventually-unblinded Phase 1) would allow:
  - Error averaging across independent paths
  - Reduced dependence on any single anchor's experimental value
  - Reduced dependence on any single accumulated-ddG path
- No cross-checking of propagated values was possible with one anchor.

### A4. Cluster boundaries fixed before network construction
- Tanimoto 0.5 clustering set membership before edges were planned.
- Structurally close Phase 1 / Phase 2 compounds in different clusters could
  never be connected by an edge.
- A clustering (or post-clustering merge) ensuring each Phase 2 compound shares
  a cluster with >=1 (eventually-known) Phase 1 compound would have guaranteed
  a high-quality local anchor.

### A5. Propagation error grows with path length
- Absolute pEC50 = anchor pEC50 + sum of ddG along path.
- Longer paths accumulate more per-edge error (and uncertainty was 0.0 from
  n_protocol_repeats=1, so we can't even quantify it well — see B4).
- MST tends to produce longer paths than a redundant network would for the
  same node set. No averaging to damp the accumulation.

---

## B. Protocol / input lessons from failed runs

### B1. Atom-mapping strictness drove the two dominant failure modes
- 54 IndexError = empty Kartograf mapping (mapped_old_atom_indices[0] on an
  empty list at HybridTopologyFactory). Kartograf's atom_max_distance=0.95 +
  map_exact_ring_matches_only=true rejected ALL atom correspondences for hard
  pairs.
- 189 SimulationNaNError = strained/element-changing hybrid topology that built
  but exploded during equilibration.
- Lesson: a post-planning validation step should have checked every
  transformation JSON for a non-empty, reasonable-coverage mapping BEFORE
  submitting 1066 GPU jobs. Empty/low-coverage mappings are detectable in
  seconds and never worth a GPU-hour.

### B2. Complex legs fail; solvent legs almost never do
- Final: solvent 503/533 done, complex 315/533. Every failed edge is missing
  its complex leg.
- NaN originates in the protein-bound environment (pocket clashes in the docked
  pose, strained hybrid atoms near the binding site), not in the mapping per se.
- LOMAP salvage with element_change=False confirmed this: 90/94 salvage complex
  legs STILL NaN'd. The mapper was not the root cause — the starting geometry was.

### B3. Docked poses were never force-field-minimized before RBFE
- Poses came straight from GNINA docking. PDBFixer prepared the protein, but
  individual ligand poses were not relaxed in the actual MD force field.
- Subtle clashes that GNINA tolerates produce NaN at lambda intermediates.
- Lesson: a short per-pose minimization (or restrained equilibration ramp) in
  the production force field, before alchemical setup, would likely have
  rescued a large fraction of the 189 NaN edges. This is the highest-value
  protocol change for a re-run.

### B4. Single protocol repeat => zero uncertainty
- n_protocol_repeats=1 gave uncertainty.magnitude = 0.0 on every leg.
- We cannot quantify per-edge statistical error, can't weight edges in
  propagation, and can't report calibrated confidence intervals.
- openfe best practice is >=3 repeats. Trade-off was throughput vs. error bars;
  with the compute we ended up having, 2-3 repeats were affordable.

### B5. Equilibration protocol may be too aggressive for docked poses
- NaN occurred during equilibrate(), before production sampling.
- Gentler options not used: more minimization steps, slower heating schedule,
  smaller initial timestep ramped up, hydrogen mass repartitioning tuning.

### B6. No element-change guard in the primary run
- Kartograf default permitted element changes (N->H, O->H seen in warnings),
  which create unstable hybrid atoms.
- The salvage used element_change=False; should have been the default for the
  primary campaign.

### B7. Sampling length not tuned to difficulty
- Fixed 2000 iterations / ~5 ns production for every edge regardless of
  perturbation size.
- early_termination_target_error (converge-and-stop) was available and unused;
  would have freed compute from easy edges to spend on repeats/harder edges.

---

## C. What actually worked (keep for next time)
- Multi-pool HTCondor (GPU Lab + backfill + campus + OSPool) gave 100+
  concurrent GPUs; full 1066-leg campaign finished in <24h wall time.
- CUDA fast-fail check correctly shed broken server slots for retry.
- Dropping checkpointing once we knew legs were short (mean 1.6h, max <12h)
  simplified the pipeline and removed the -o file-conflict bug class.
- Provenance tiering (Kartograf tier-1 vs LOMAP-salvage tier-2) keeps the
  mixed-mapper question auditable for calibration.
- Preserving all failed-edge results enables total-compute accounting and
  reproducibility.

---

## D. Net impact statement (draft for writeup)
The campaign was designed as a conventional RBFE study treating all 513 test
compounds uniformly. It was methodologically clean but did not exploit the
two-phase challenge structure (no Phase-1-aware anchors, single anchors,
fragile MST topology). Combined with a ~47% complex-leg NaN rate driven by
un-minimized docked poses and permissive atom mapping, this left 43% of
compounds without an RBFE path, forcing docking fallback. The connectivity
analysis (0 Phase-2 compounds recoverable by adding Phase 1 anchors post hoc)
shows the loss was baked into the topology, not fixable after the runs.
Highest-value changes for a future run, in order: (1) redundant network
topology, (2) per-pose force-field minimization before alchemy,
(3) >=3 protocol repeats for real uncertainty, (4) post-planning mapping
validation gate, (5) Phase-1-aware anchor/cluster design.
