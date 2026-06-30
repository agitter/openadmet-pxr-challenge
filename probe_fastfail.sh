#!/usr/bin/env bash
# probe_fastfail.sh
# Investigate the CUDA fast-fail markers in .out/.err files so we can
# reliably tag each execution attempt as fast-fail-rejection vs real-run.
# Run from the repo root.

echo "============================================================"
echo "1. FIND A SHORT-DURATION ATTEMPT (likely fast-fail)"
echo "============================================================"
echo "gpu4002 had ~25s mean attempt duration - inspect one of its jobs."
LOG=$(grep -rl "gpu4002.chtc.wisc.edu" openfe/production --include="*.log" 2>/dev/null | head -1)
echo "Log: $LOG"
OUT="${LOG%.log}.out"
ERR="${LOG%.log}.err"

echo ""
echo "=== .out: CUDA / platform markers ==="
grep -i "cuda\|platform\|fast.fail\|not available\|confirmed\|gpu" "$OUT" 2>/dev/null | head -20

echo ""
echo "=== .out: first 60 lines (full context) ==="
head -60 "$OUT" 2>/dev/null

echo ""
echo "=== .err: CUDA / failure markers ==="
grep -i "cuda\|platform\|not available\|nan\|error\|fail" "$ERR" 2>/dev/null | head -20

echo ""
echo "============================================================"
echo "2. FREQUENCY OF CUDA MARKERS ACROSS ALL .out FILES"
echo "============================================================"
echo -n "  'CUDA platform confirmed available' (success marker): "
grep -rl "CUDA platform confirmed available" openfe/production --include="*.out" 2>/dev/null | wc -l
echo -n "  'CUDA platform not available' (fast-fail marker):     "
grep -rl "CUDA platform not available" openfe/production --include="*.out" 2>/dev/null | wc -l

echo ""
echo "=== Any other CUDA-related phrasings present (counts) ==="
grep -rh -i "cuda platform\|cuda.*not\|no cuda\|fast.fail\|platform not\|platform confirmed" \
    openfe/production --include="*.out" 2>/dev/null | sed -E 's/^[[:space:]]+//' \
    | sort | uniq -c | sort -rn | head -15

echo ""
echo "============================================================"
echo "3. SAME CHECK IN .err FILES (markers may go to stderr)"
echo "============================================================"
echo -n "  .err with 'CUDA platform not available': "
grep -rl "CUDA platform not available" openfe/production --include="*.err" 2>/dev/null | wc -l
echo -n "  .err with 'CUDA platform confirmed available': "
grep -rl "CUDA platform confirmed available" openfe/production --include="*.err" 2>/dev/null | wc -l
echo ""
echo "=== exact fast-fail-related lines from .err (samples) ==="
grep -rh -i "cuda platform\|platform not available\|platform confirmed" \
    openfe/production --include="*.err" 2>/dev/null | sed -E 's/^[[:space:]]+//' \
    | sort | uniq -c | sort -rn | head -15

echo ""
echo "============================================================"
echo "4. HOW MANY .out FILES HAVE NEITHER MARKER?"
echo "============================================================"
TOTAL=$(find openfe/production -name "*.out" 2>/dev/null | wc -l)
WITH_CONF=$(grep -rl "CUDA platform confirmed available" openfe/production --include="*.out" 2>/dev/null | wc -l)
WITH_FAIL=$(grep -rl "CUDA platform not available" openfe/production --include="*.out" 2>/dev/null | wc -l)
echo "  total .out files:           $TOTAL"
echo "  with 'confirmed available': $WITH_CONF"
echo "  with 'not available':       $WITH_FAIL"
echo "  (a single .out may contain multiple attempts, so these can overlap)"

echo ""
echo "============================================================"
echo "5. DOES A SINGLE .out CONTAIN MULTIPLE CUDA CHECKS?"
echo "============================================================"
echo "Count occurrences of the markers within one multi-attempt .out:"
MULTI=$(grep -rl "CUDA platform not available" openfe/production --include="*.out" 2>/dev/null | head -1)
echo "Sample .out: $MULTI"
echo -n "  'confirmed available' count in it: "
grep -c "CUDA platform confirmed available" "$MULTI" 2>/dev/null
echo -n "  'not available' count in it:       "
grep -c "CUDA platform not available" "$MULTI" 2>/dev/null
echo ""
echo "(If a .out aggregates all attempts of a leg, marker COUNTS tell us"
echo " how many fast-fail rejections vs confirmed runs happened per leg.)"
