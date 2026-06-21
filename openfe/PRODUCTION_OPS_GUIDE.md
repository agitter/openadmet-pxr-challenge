# OpenFE RBFE Production Campaign — Operations Guide

## Overview

1066 RBFE transformation legs across 124 clusters, running on CHTC GPU
infrastructure (GPU Lab + backfill + campus pools + OSPool). Expected
runtime: 5-8 days at 20-30 concurrent GPUs. Each leg takes 2.5-5.5
hours depending on GPU type.

## Quick Reference

```bash
# Check progress (run from repo root)
python openfe/scripts/06_monitor_production.py --production-dir openfe/production

# Check HTCondor queue
condor_q                    # summary
condor_q -held              # held jobs
condor_q -run               # running jobs
condor_q -idle              # idle/waiting jobs
```

## Phase 1: Setup and Launch

### 1.1 Create job directories

```bash
python openfe/scripts/05_setup_production.py \
    --network-results openfe/results \
    --script openfe/scripts/run_quickrun.sh \
    --outdir openfe/production
```

### 1.2 Verify setup

```bash
# Should show 1066
wc -l openfe/production/transform_list.txt

# Inspect one job directory
ls openfe/production/118/rbfe_OCNT-2310728_complex_OADMET-0006495_complex/
# Expect: run_quickrun.sh  rbfe_...json  network_setup.json  quickrun_output/  result.json
```

### 1.3 Submit

```bash
cd openfe
condor_submit submit_quickrun_production.sub
```

Note the HTCondor cluster ID printed (e.g., `1066 job(s) submitted to cluster 8012345`).
You'll need this for `condor_rm` if you need to remove jobs later.

## Phase 2: Monitoring

### How often to check

| Period | Action |
|---|---|
| First 30 minutes | Check that jobs are starting: `condor_q -run` |
| First 2 hours | Check first completions and any early failures: `condor_q -held` |
| Every 4-8 hours | Run the monitor script for progress and storage |
| Once per day | Review held jobs, decide on resubmission |

The campaign runs unattended between checks. HTCondor automatically:
- Retries failed jobs (up to 10 times per submission)
- Reschedules evicted jobs (unlimited, with checkpoint transfer)
- Releases held jobs after 10 minutes (up to 20 total starts)

### 2.1 Run the monitor

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
- A retry list (`production/transform_list_retry.txt`) if any jobs are incomplete

### 2.2 Interpret the status

| Status | Meaning | Action |
|---|---|---|
| `COMPLETED` | Job finished, result.json valid, checkpoints cleaned | None |
| `COMPLETED_NEEDS_CLEANUP` | Job finished but large files remain | Run with `--cleanup` |
| `IN_PROGRESS_OR_FAILED` | Has checkpoint files, no valid result | Check if still running (see below) |
| `STARTED` | Directory exists but minimal activity | May be queued or very early |
| `NOT_STARTED` | Empty job directory | Job hasn't run yet |

### 2.3 Distinguish running from failed

`IN_PROGRESS_OR_FAILED` could mean the job is actively running or has
permanently failed. To tell the difference:

```bash
# Check if the job is still in the HTCondor queue
condor_q    # shows running/idle/held counts

# If condor_q shows 0 jobs and monitor shows IN_PROGRESS_OR_FAILED,
# those jobs have permanently failed (exhausted retries, removed from queue).
```

## Phase 3: Handling Failures

### Decision tree

```
Monitor shows incomplete jobs
│
├─ condor_q shows jobs still running/idle?
│  └─ YES → Wait. Jobs are still in progress. Check again in 4-8 hours.
│
├─ condor_q shows held jobs?
│  ├─ Few held jobs (< 20)?
│  │  └─ Check hold reason: condor_q -held
│  │     ├─ Transfer failure → Will auto-release (periodic_release). Wait.
│  │     ├─ Exceeded max retries → Go to "Resubmit failed jobs" below.
│  │     └─ Other reason → Investigate condor.err in the job's initialdir.
│  │
│  └─ Many held jobs (> 20)?
│     └─ Likely a systemic issue (bad container, all servers down).
│        Check one job's condor.err, fix the issue, then resubmit all.
│
└─ condor_q shows 0 jobs but monitor shows incomplete?
   └─ Jobs were removed or exhausted all retries.
      Go to "Resubmit failed jobs" below.
```

### 3.1 Resubmit failed jobs

The monitor automatically writes `production/transform_list_retry.txt`
with all incomplete jobs. To resubmit:

```bash
# Remove any remaining held jobs from the old submission
condor_rm <old_cluster_id>

# Resubmit just the incomplete jobs
cd openfe
condor_submit submit_quickrun_production.sub \
    -append "queue cluster_id,transform_name from production/transform_list_retry.txt"
```

Checkpoint files are still in each job's initialdir, so `--resume`
picks up from the last checkpoint. The resubmission gets a fresh retry
counter (10 more attempts).

### 3.2 Resubmit with server exclusions

If specific servers consistently fail (e.g., CUDA container issues):

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

### 3.3 Resubmit with modified settings

If jobs fail due to resource limits (out of memory, disk, time):

```bash
cd openfe
condor_submit submit_quickrun_production.sub \
    -append "request_memory = 32GB" \
    -append "queue cluster_id,transform_name from production/transform_list_retry.txt"
```

Or if you decide to replan with shorter sampling (see Phase 5 below).

### 3.4 Truly un-runnable transformations

Some edges may fail permanently regardless of retries (bad atom
mappings, extreme perturbations). Since the network is a tree (MST),
losing an edge disconnects a compound from the anchor. These compounds
will fall back to docking-based predictions (CNNaffinity) in the final
analysis. Accept and move on — a partial RBFE submission is better
than no submission.

## Phase 4: Storage Management

### 4.1 Clean up completed jobs with leftover large files

```bash
python openfe/scripts/06_monitor_production.py \
    --production-dir openfe/production --cleanup
```

Safe to run at any time. Only affects completed jobs.

### 4.2 Clean up permanently failed jobs

**Only run after confirming jobs are truly done retrying:**

```bash
# First confirm no matching jobs in queue
condor_q

# Then clean up
python openfe/scripts/06_monitor_production.py \
    --production-dir openfe/production --cleanup-failed
```

**WARNING:** This deletes checkpoint files. Failed jobs cannot be
resumed after cleanup. Only use after you've given up on those jobs
or resubmitted them and they completed on the retry.

### Storage budget

| Item | Size |
|---|---|
| Job directories (setup) | ~1.6 GB |
| Active checkpoints (during run) | ~4.5 GB peak |
| Completed results (all 1066) | ~3.2 GB |
| AP quota remaining | 74 GB |
| **Comfortable margin** | **~65 GB** |

## Phase 5: If Running Out of Time

If after 2-3 days, throughput is too low to finish all 1066 jobs
in time:

### Option A: Replan with shorter sampling (recommended first)

Reduces per-leg time from ~3.5h to ~2h by cutting production length.

1. Update `plan_settings.yaml` with shorter simulation settings
2. Re-run `condor_submit submit_plan_network.sub` (2h on CHTC CPU)
3. Re-run `05_setup_production.py` (seconds)
4. Resubmit incomplete jobs with new transformation JSONs

**Important:** completed jobs keep their existing results. Only
incomplete jobs need to be rerun with the new settings.

### Option B: Prioritize high-value clusters

Focus GPU time on clusters with the most test compounds. The monitor
shows per-cluster completion. Clusters with 1-2 compounds (singletons)
contribute fewer predictions per GPU-hour than large clusters.

```bash
# See which large clusters are incomplete
python openfe/scripts/06_monitor_production.py --production-dir openfe/production
# Look at "Clusters fully completed" section
```

### Option C: Accept partial results

Any cluster with ALL edges completed contributes predictions. Clusters
with some edges missing can still contribute partial predictions
(compounds reachable from an anchor via completed edges). Completely
failed clusters fall back to docking-based predictions.

## Phase 6: After All Jobs Complete

```bash
# Final status check
python openfe/scripts/06_monitor_production.py --production-dir openfe/production

# Final cleanup
python openfe/scripts/06_monitor_production.py \
    --production-dir openfe/production --cleanup

# Verify storage
du -sh openfe/production/
```

Then proceed to results gathering, pEC50 propagation, and submission
preparation (separate workflow).
