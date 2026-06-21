# OpenFE RBFE Production Campaign — Operations Guide

## Overview

1066 RBFE transformation legs across 124 clusters, running on CHTC GPU
infrastructure (GPU Lab + backfill + campus pools + OSPool). Expected
runtime: 5-8 days at 20-30 concurrent GPUs. Each leg takes 2.5-5.5
hours depending on GPU type.

**All `condor_*` commands must run on the CHTC access point (ap2001).**
Python scripts can run on the AP or locally.

## Quick Reference

```bash
# Check progress (from repo root)
python openfe/scripts/06_monitor_production.py --production-dir openfe/production

# Check HTCondor queue (on AP)
condor_q                    # summary
condor_q -held              # held jobs
condor_q -run               # running jobs
```

---

## Phase 1: Setup and Launch

All steps from the **repo root** unless noted.

### 1.1 Create job directories

```bash
python openfe/scripts/05_setup_production.py
```

This reads `openfe/results/network_setup_*/transformations/*.json` and creates
one directory per transformation under `openfe/production/`, each containing
the transformation JSON, network_setup.json, run_quickrun.sh, an empty
quickrun_output/, and an empty result.json.

### 1.2 Verify setup

```bash
# Should show 1066
wc -l openfe/production/transform_list.txt

# Inspect a sample job directory — should contain 5 items
ls openfe/production/118/rbfe_OCNT-2310728_complex_OADMET-0006495_complex/
# Expect: run_quickrun.sh  rbfe_...json  network_setup.json  quickrun_output/  result.json
```

### 1.3 Submit (on AP)

```bash
cd openfe
mkdir -p logs
condor_submit submit_quickrun_production.sub
```

Note the HTCondor cluster ID printed (e.g., `1066 job(s) submitted to cluster 8012345`).
You will need this for `condor_rm` if you need to remove jobs later.

---

## Phase 2: Monitoring

### How often to check

| Period | Action |
|---|---|
| First 30 minutes | Check that jobs are starting: `condor_q -run` |
| First 2 hours | Check for early failures: `condor_q -held` |
| Every 4-8 hours | Run the monitor script for progress and storage |
| Once per day | Review held jobs, decide on resubmission |

The campaign runs unattended between checks. HTCondor automatically:
- Retries failed jobs (up to 10 times per submission)
- Reschedules evicted jobs (unlimited, with checkpoint transfer)
- Releases held jobs after 10 minutes (up to 20 total starts)

### 2.1 Run the monitor (from repo root)

```bash
python openfe/scripts/06_monitor_production.py --production-dir openfe/production
```

This prints:
- Overall completion count and percentage
- Complex vs solvent leg breakdown
- Storage usage on AP
- Timing statistics for completed jobs
- Remaining time estimates at different GPU counts
- Per-cluster completion status

If any jobs are incomplete, it writes `openfe/production/transform_list_retry.txt`
with resubmission instructions.

### 2.2 Interpret the status

| Status | Meaning | Action needed? |
|---|---|---|
| `COMPLETED` | result.json valid, checkpoints cleaned | No |
| `COMPLETED_NEEDS_CLEANUP` | result.json valid, large files remain | Run with `--cleanup` |
| `IN_PROGRESS_OR_FAILED` | Has checkpoint files, no valid result | See 2.3 below |
| `STARTED` | Directory exists, minimal activity | May be queued or early in run |
| `NOT_STARTED` | Empty job directory | Job hasn't run yet |

### 2.3 Distinguish running jobs from failed jobs

`IN_PROGRESS_OR_FAILED` means the job directory has checkpoint files but
no valid result. It could be actively running or permanently failed.

```bash
# On the AP: check if jobs are still in the HTCondor queue
condor_q

# If condor_q shows running/idle jobs → wait, check again in 4-8 hours
# If condor_q shows 0 jobs and monitor shows IN_PROGRESS_OR_FAILED →
#   those jobs have permanently failed. Go to Phase 3.
```

---

## Phase 3: Handling Failures

### Decision tree

```
Monitor shows incomplete jobs
│
├─ condor_q shows jobs running or idle?
│  └─ YES → Wait. Check again in 4-8 hours.
│
├─ condor_q shows held jobs?
│  │
│  ├─ A few held (< 20)?
│  │  └─ Check why: condor_q -held
│  │     ├─ Transfer failure → Will auto-release in 10 min. Wait.
│  │     ├─ Max retries exhausted → Go to 3.1 (Resubmit).
│  │     └─ Other → Check the job's condor.err file (in its initialdir).
│  │
│  └─ Many held (> 20)?
│     └─ Systemic issue. Check one job's condor.err.
│        Fix the root cause, then go to 3.1 (Resubmit all failed).
│
└─ condor_q shows 0 jobs, but monitor shows incomplete?
   └─ Jobs exhausted all retries or were removed.
      Go to 3.1 (Resubmit).
```

### 3.1 Resubmit failed jobs

The monitor automatically writes `openfe/production/transform_list_retry.txt`
listing all incomplete jobs. To resubmit:

```bash
# On AP: remove any remaining held jobs from the old submission
condor_rm <old_cluster_id>

# Resubmit just the incomplete jobs (from openfe/ directory)
cd openfe
condor_submit submit_quickrun_production.sub \
    -append "queue cluster_id,transform_name from production/transform_list_retry.txt"
```

Checkpoint files are still in each job's initialdir. The `--resume` flag
picks up from the last checkpoint. The new submission gets a fresh retry
counter (10 more attempts).

### 3.2 Resubmit with server exclusions

If specific servers consistently fail (e.g., CUDA not exposed to container):

```bash
cd openfe
condor_submit submit_quickrun_production.sub \
    -append 'requirements = (Machine != "bad-server.chtc.wisc.edu")' \
    -append "queue cluster_id,transform_name from production/transform_list_retry.txt"
```

To exclude multiple servers:

```bash
cd openfe
condor_submit submit_quickrun_production.sub \
    -append 'requirements = (Machine != "server1.chtc.wisc.edu") && (Machine != "server2.chtc.wisc.edu")' \
    -append "queue cluster_id,transform_name from production/transform_list_retry.txt"
```

### 3.3 Resubmit with modified resource requests

If jobs fail due to out of memory, insufficient disk, or timeout:

```bash
cd openfe
condor_submit submit_quickrun_production.sub \
    -append "request_memory = 32GB" \
    -append "queue cluster_id,transform_name from production/transform_list_retry.txt"
```

### 3.4 Truly un-runnable transformations

Some edges may fail permanently regardless of retries (bad atom mappings,
extreme perturbations, charge-change issues). Since the network is a
minimum spanning tree, every edge is a bridge — losing one disconnects
at least one compound from the anchor.

For the competition, compounds unreachable via RBFE fall back to
docking-based predictions (CNNaffinity from the extended docking campaign).
A partial RBFE submission with docking fallbacks is better than no
submission.

---

## Phase 4: Storage Management

### 4.1 Clean up completed jobs with leftover large files

```bash
python openfe/scripts/06_monitor_production.py \
    --production-dir openfe/production --cleanup
```

Safe to run at any time. Only deletes large checkpoint files from jobs
that already have a valid result.json.

### 4.2 Clean up permanently failed jobs

**Only run after you are sure these jobs will not be retried:**

```bash
# First confirm no matching jobs remain in the queue
condor_q

# Then clean up checkpoint files from failed jobs
python openfe/scripts/06_monitor_production.py \
    --production-dir openfe/production --cleanup-failed
```

**WARNING:** This deletes checkpoint files permanently. Failed jobs
cannot be resumed after cleanup. Only use after either:
- You have resubmitted and the retries completed successfully, OR
- You have decided to accept the loss and use docking fallbacks

### Storage budget

| Item | Size | When |
|---|---|---|
| Job directories (setup) | ~1.6 GB | After step 1.1 |
| Active checkpoints | ~4.5 GB peak | During run (30 concurrent jobs) |
| Completed results | ~3.2 GB total | Accumulates as jobs finish |
| **Total peak** | **~9.3 GB** | |
| **AP quota remaining** | **~65 GB** | Comfortable margin |

---

## Phase 5: If Running Out of Time

If after 2-3 days the monitor's time estimates show you won't finish
in time, consider these options in order of preference:

### Option A: Replan with shorter sampling (recommended)

Reduces per-leg time from ~3.5h to ~2h by halving production length.

1. Update `openfe/plan_settings.yaml` with shorter simulation settings
2. Re-run the network planning jobs on CHTC CPU (~2 hours)
3. Delete `openfe/production/` directories for incomplete jobs only:
   ```bash
   # The retry list has the incomplete job names
   # Delete their directories so setup can recreate them with new JSONs
   while IFS=, read -r cid tname; do
       rm -rf "openfe/production/$cid/$tname"
   done < openfe/production/transform_list_retry.txt
   ```
4. Re-run `python openfe/scripts/05_setup_production.py` to create
   new directories with the updated transformation JSONs
5. Resubmit incomplete jobs

Completed jobs keep their existing results. Only incomplete jobs
are affected.

### Option B: Prioritize high-value clusters

Focus GPU time on clusters with the most test compounds. Singleton
clusters (1 test compound + 1 anchor) contribute 1 prediction per
2 GPU-legs. Large clusters (10+ compounds) contribute 10+ predictions.

```bash
python openfe/scripts/06_monitor_production.py --production-dir openfe/production
# Check "Clusters fully completed" and per-cluster status
```

### Option C: Accept partial results

- Clusters with all edges completed: full RBFE predictions
- Clusters with some edges completed: partial predictions (compounds
  still connected to an anchor via completed edges)
- Clusters with no completed edges: fallback to docking predictions

---

## Phase 6: After All Jobs Complete

```bash
# Final status check (from repo root)
python openfe/scripts/06_monitor_production.py --production-dir openfe/production

# Clean up any remaining large checkpoint files
python openfe/scripts/06_monitor_production.py \
    --production-dir openfe/production --cleanup

# Verify total storage
du -sh openfe/production/
```

Then proceed to results gathering, pEC50 propagation from anchor
compounds, isotonic calibration against Phase 1 unblinded data, and
final submission preparation. That workflow will be documented
separately.
