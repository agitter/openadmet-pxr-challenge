## Repository Organization

A map of the repository. Paths point to where the important scripts and
results reside.

### Top level

- `README.md` — project overview and this repository map.
- `methods_outline.md` — step-by-step outline of the pipeline.
- `commands.md` — reference commands used during the project.
- `LICENSE`, `.gitignore`, `.gitmodules` — standard project files.

### `/data`

Challenge data files from
[Hugging Face](https://huggingface.co/datasets/openadmet/pxr-challenge-train-test):
training set, blinded test set, and the Phase 1 unblinded test set used for
calibration. Downloaded via `download-hf-data.py`.

### `/claude`

Work produced in an assistant-run session (structure discovery, test-set
clustering, and anchor assignment), kept separate for provenance.

- `claude/structure_discovery/` — scripts to find, analyze, and select PXR
  crystal structures and characterize the binding-pocket ensemble.
- `claude/analysis/` — test-set clustering, training-anchor identification,
  and cluster-to-crystal-template mapping.
- `claude/outputs/` — resulting artifacts: structure inventory, fingerprints
  and similarity matrix, clustered test set, per-cluster representatives and
  summaries, training-anchor assignments, and cluster-template mapping.
- `claude/limitations_and_lessons.md` — limitations and lessons-learned
  notes for the writeup.
- `claude/README.md` — notes on the assistant-run session.

### `/docking`

GNINA docking pipeline and results.

- `docking/scripts/` — receptor preparation, ligand embedding, work-unit
  construction, the docking driver, result aggregation, and RBFE
  template/pose selection.
- `docking/receptors/`, `docking/ligands/` — prepared receptor structures
  and embedded ligands.
- `docking/work_units/`, `docking/results/`, `docking/logs/` — HTCondor
  work units, raw docking outputs, and job logs.
- `docking/rbfe_inputs/` — per-cluster inputs extracted for RBFE.
- `docking/docking_analysis/`, `docking/docking_analysis_extended/` —
  aggregated docking scores; the extended analysis holds the per-compound
  best-pose scores, the per-cluster summary (including cross-receptor score
  statistics), and the Phase 1 compounds joined to their docking scores.
- `docking/submit_docking.sub` — docking job submission file.

### `/external`

Third-party inputs.

- Re-refined PXR structures as a
  [submodule](https://github.com/OpenADMET/pxr_xtal_re-refinement).
- Organizer-provided evaluation and validation code from the
  [PXR Challenge Tutorial](https://github.com/OpenADMET/PXR-Challenge-Tutorial):
  the scoring metric and evaluation script, bootstrap utilities, and the
  submission-format validator.

### `/openfe`

[OpenFE](https://github.com/OpenFreeEnergy/openfe) RBFE pipeline: planning,
production execution, salvage, all downstream analysis, and the final
submission.

- `openfe/scripts/` — the numbered pipeline and analysis scripts, in
  approximate order of execution:
  - Network inputs, receptor identification, receptor preparation, RBFE
    input extraction, and planning-input preparation (early `00`–`04`).
  - Production setup, monitoring, result consolidation, and reporting
    (`05`–`08`).
  - Connectivity analysis (`09`), salvage planning (`10`), and combined
    results gathering with provenance (`11`).
  - pEC50 propagation from anchors (`12`) and the Phase 1 anchor-value
    analysis (`13`).
  - Convergence analysis (`15`), edge-reliability threshold sweep (`16`),
    docking-versus-RBFE method comparison (`17`), and data-coverage audit
    (`18`).
  - Combined-model tuning and cross-validation (`19`), blend-weight
    selection (`20`), and final model application to produce the submission
    (`21`).
  - Compute accounting (`22`).
  - Supporting utilities: transformation inventory, result rescue, structure
    probing, and quickrun timing analysis.
- `openfe/plan_settings.yaml` — network-planning settings.
- `openfe/run_quickrun.sh`, `openfe/run_plan_network.sh`,
  `openfe/run_quickrun_timing.sh` — execution wrappers for production legs,
  network planning, and the timing pilot.
- `openfe/submit_*.sub`, `openfe/openmm_test.*` — HTCondor submission files
  for planning, production, salvage, the timing pilot, and an OpenMM
  environment test.
- `openfe/production/` — per-transform working directories with free energy
  results, salvage results, and job logs (large; retained for auditing).
- Key result files (top level of `openfe/`): the combined per-edge results
  with provenance, connectivity report, RBFE-propagated predictions with
  path features, per-leg and per-edge convergence tables, method-comparison
  metrics, the cross-validation model comparison, the frozen chosen model,
  and the final submission and its supporting detail and visualizations.
- `openfe/PRODUCTION_OPS_GUIDE.md` — operational notes for the production
  campaign.

### `/analysis`

Standalone early analysis (docking-score versus potency correlation).

### `/submissions`

Participant method-summary links and tooling to retrieve them, scraped from
the challenge's
[Hugging Face app](https://openadmet-pxr-challenge.hf.space/config):
a parser for the submission table and a browser-based downloader for the
linked writeups.
