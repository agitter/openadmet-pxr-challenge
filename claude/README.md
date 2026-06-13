# PXR Challenge — Analysis & Pipeline Bundle

Snapshot of everything produced so far for the activity-track RBFE
strategy (ensemble docking -> anchor-network RBFE -> isotonic
calibration vs Phase 1).

## analysis/ — standalone, re-runnable scripts

- **01_cluster_test_set.py** — ECFP4/Tanimoto single-linkage clustering
  of the 513-compound blinded test set. Produces `test_with_clusters.csv`,
  `cluster_summary.csv`, `test_tanimoto_sim.npy`, `test_fps.pkl`.
  Default threshold 0.5 -> 124 clusters, 54 singletons, 489/513 compounds
  in clusters of size >= 2.

- **02_find_training_anchors.py** — for each test compound, finds its
  nearest-neighbor training compound (by ECFP4 Tanimoto) and that
  compound's known pEC50. Produces `test_with_train_anchors.csv`.
  Result: 504/513 test compounds have a training analog at >=0.4
  similarity; 355/513 at >=0.5 (those cluster tightly around pEC50
  ~5.9-6.0, i.e. the original 63 hits / their close analogs).

- **03_cluster_to_crystal_template.py** — for each test cluster, finds
  the most similar bound ligand among the 64 re-refined PXR crystal
  structures. Produces `cluster_template_mapping.csv`.
  Result: best similarity across ALL clusters was ~0.12-0.21 — i.e.
  NO crystal ligand is chemically similar to the test-set chemotype.
  This ruled out rigid-template pose transfer and motivated the
  ensemble-docking approach.

## outputs/ — data files produced by the analysis scripts (already run)

- `test_with_clusters.csv` — 513 test compounds + `cluster_id`
  (threshold 0.5, 124 clusters)
- `cluster_representatives.csv` / `cluster_template_mapping.csv` — one
  row per cluster: size, representative compound, best crystal-template
  match (see note above — none are good matches)
- `test_with_train_anchors.csv` — 513 test compounds + nearest training
  analog + that analog's pEC50 (calibration anchor candidates)
- `test_tanimoto_sim.npy` — full 513x513 ECFP4 Tanimoto similarity matrix
  (re-load with `np.load(...)`)
- `test_fps.pkl` — pickled fingerprints + cluster assignments
- `test_valid.csv` — test compounds with successfully-parsed SMILES
- `pxr_structure_inventory.csv` / `pxr_pocket_residues.json` — from the
  64 re-refined PXR crystal structures (ligand identity, pocket
  residues, resolution etc.) — output of
  `structure_discovery/analyze_pxr_ensemble.py` run against the
  `pxr_xtal_re-refinement` submodule.

## structure_discovery/ — scripts used to locate/characterize PXR structures

- `analyze_pxr_structures.py`, `find_pxr_structures.py`,
  `fetch_pxr_rerefinement.py` — exploratory scripts used while locating
  the 64-structure re-refined PXR ensemble (superseded once the
  `pxr_xtal_re-refinement` GitHub submodule was identified — kept for
  reference / in case the 184-structure structure-track release needs
  similar triage).
- `analyze_pxr_ensemble.py` — the script actually used to produce
  `pxr_structure_inventory.csv` and `pxr_pocket_residues.json` from the
  submodule's 64 `.pdb`/`.cif` files.

## docking/ — GNINA ensemble docking pipeline (CHTC/HTCondor)

(Removed)

## Pipeline order (to reproduce from raw HF data)

1. `analysis/01_cluster_test_set.py` (needs `pxr-challenge_TEST_BLINDED.csv`)
2. `analysis/02_find_training_anchors.py` (needs TEST + `pxr-challenge_TRAIN.csv`)
3. `git submodule add https://github.com/OpenADMET/pxr_xtal_re-refinement.git external/pxr_xtal_re-refinement`
4. `structure_discovery/analyze_pxr_ensemble.py` (run from repo root, needs the submodule)
5. `analysis/03_cluster_to_crystal_template.py` (needs outputs of 1 and 4)
6. `docking/scripts/prep_receptors.py` then `docking/scripts/make_work_units.py`
7. `condor_submit docking/submit_docking.sub`
