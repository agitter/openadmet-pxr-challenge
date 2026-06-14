#!/usr/bin/env python3
"""
docking/scripts/select_rbfe_templates.py

For each of the 124 cluster-representative ligands, select a single
(receptor, pose) to use as the RBFE starting structure, from the
ensemble-docking results aggregated by aggregate_docking_results.py.

Selection rule:
  Rank candidates (cluster, receptor, pose) by CNNscore (descending).
  Take the highest-CNNscore candidate UNLESS its minimizedAffinity is
  above --max-affinity (default 0.0 kcal/mol, i.e. energetically
  unfavorable by Vina's empirical scoring - a sign of a clashy/strained
  pose that the CNN scorer can sometimes still rate confidently). In
  that case, fall back to the next-highest-CNNscore candidate that
  passes the affinity filter.

  If NO candidate for a cluster passes the affinity filter (shouldn't
  happen given the observed data, but handled defensively), fall back
  to the single best-CNNscore candidate regardless, and flag it.

Inputs:
  docking_poses_long.csv   (from aggregate_docking_results.py - needed
                             because the fallback may pick a pose other
                             than each receptor's single best-CNNscore
                             pose recorded in docking_receptor_best.csv)
  cluster_representatives.csv (cluster_id, representative_test_compound,
                             representative_smiles)

Output:
  rbfe_template_selection.csv - one row per cluster:
    cluster_id, ligand_name, smiles,
    selected_pdb_id, selected_pose_rank,
    CNNscore, CNNaffinity, minimizedAffinity,
    used_fallback (bool), fallback_reason

Usage:
    pip install pandas
    python select_rbfe_templates.py \
        --poses-long docking_analysis/docking_poses_long.csv \
        --cluster-representatives cluster_representatives.csv \
        --max-affinity 0.0 \
        --outdir docking_analysis
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poses-long", required=True,
                    help="docking_poses_long.csv from "
                         "aggregate_docking_results.py")
    ap.add_argument("--cluster-representatives", required=True,
                    help="cluster_representatives.csv "
                         "(cluster_id, representative_test_compound, "
                         "representative_smiles)")
    ap.add_argument("--max-affinity", type=float, default=0.0,
                    help="Reject candidates with minimizedAffinity above "
                         "this value (kcal/mol; default: 0.0, i.e. "
                         "exclude energetically unfavorable poses)")
    ap.add_argument("--outdir", default="docking_analysis")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    poses = pd.read_csv(args.poses_long)
    reps = pd.read_csv(args.cluster_representatives)

    print(f"Loaded {len(poses)} pose rows across "
          f"{poses['cluster_id'].nunique()} clusters")
    print(f"Loaded {len(reps)} cluster representatives")

    rows = []
    n_fallback = 0
    n_no_pass = 0

    for cluster_id, group in poses.groupby("cluster_id"):
        ranked = group.sort_values("CNNscore", ascending=False).reset_index(drop=True)

        passing = ranked[ranked["minimizedAffinity"] <= args.max_affinity]

        if len(passing) > 0:
            chosen = passing.iloc[0]
            used_fallback = chosen.name != 0  # i.e. not the global top row
            # more robust: compare to ranked.iloc[0]
            top = ranked.iloc[0]
            used_fallback = not (
                chosen["pdb_id"] == top["pdb_id"]
                and chosen["pose_rank"] == top["pose_rank"]
            )
            reason = (
                f"top CNNscore candidate ({top['pdb_id']}, pose "
                f"{top['pose_rank']}) had minimizedAffinity="
                f"{top['minimizedAffinity']:.2f} > {args.max_affinity}; "
                f"used next-best passing candidate"
                if used_fallback else ""
            )
        else:
            # nothing passes - fall back to the global best-CNNscore
            # candidate regardless, and flag it
            chosen = ranked.iloc[0]
            used_fallback = True
            reason = (
                f"NO candidate had minimizedAffinity <= {args.max_affinity}; "
                f"used global best-CNNscore candidate "
                f"(minimizedAffinity={chosen['minimizedAffinity']:.2f})"
            )
            n_no_pass += 1

        if used_fallback:
            n_fallback += 1

        rows.append({
            "cluster_id": cluster_id,
            "ligand_name": chosen["ligand_name"],
            "selected_pdb_id": chosen["pdb_id"],
            "selected_pose_rank": int(chosen["pose_rank"]),
            "CNNscore": chosen["CNNscore"],
            "CNNaffinity": chosen["CNNaffinity"],
            "minimizedAffinity": chosen["minimizedAffinity"],
            "used_fallback": used_fallback,
            "fallback_reason": reason,
        })

    selection = pd.DataFrame(rows).sort_values("cluster_id").reset_index(drop=True)

    # merge in SMILES
    selection = selection.merge(
        reps[["cluster_id", "representative_smiles"]],
        on="cluster_id", how="left"
    ).rename(columns={"representative_smiles": "smiles"})

    out_path = outdir / "rbfe_template_selection.csv"
    selection.to_csv(out_path, index=False)

    print(f"\nWrote {out_path} ({len(selection)} clusters)")
    print(f"Fallback used for {n_fallback} cluster(s)")
    if n_fallback:
        print(selection[selection["used_fallback"]][
            ["cluster_id", "ligand_name", "selected_pdb_id",
             "selected_pose_rank", "CNNscore", "minimizedAffinity",
             "fallback_reason"]
        ].to_string(index=False))
    if n_no_pass:
        print(f"\nWARNING: {n_no_pass} cluster(s) had NO candidate passing "
              f"the affinity filter at all - inspect these manually.")

    print(f"\nSelected-pose CNNscore distribution:")
    print(selection["CNNscore"].describe().to_string())
    print(f"\nSelected-pose minimizedAffinity distribution:")
    print(selection["minimizedAffinity"].describe().to_string())

    print(f"\nReceptor usage across selections:")
    print(selection["selected_pdb_id"].value_counts().head(15).to_string())


if __name__ == "__main__":
    main()
