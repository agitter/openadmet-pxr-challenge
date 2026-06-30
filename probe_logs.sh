#!/bin/bash
# Probe the production log layout and slot-identifying fields for the
# GPU compute-accounting tool. Run from the repo root.

echo "============================================================"
echo "1. PRODUCTION DIRECTORY STRUCTURE"
echo "============================================================"
echo "--- top-level cluster dirs (first 5) ---"
ls openfe/production/ | head -5
echo ""
FIRST=$(ls openfe/production/ | grep -E '^[0-9]+$' | head -1)
echo "Example cluster dir: $FIRST"
echo "--- transform dirs within it (first 3) ---"
ls openfe/production/$FIRST/ | head -3
echo ""
FIRST_TF=$(ls openfe/production/$FIRST/ | head -1)
echo "Example transform dir: $FIRST_TF"
echo "--- contents ---"
ls -la "openfe/production/$FIRST/$FIRST_TF/" 2>/dev/null | head -25
echo ""
echo "--- logs/ subdir if present ---"
ls "openfe/production/$FIRST/$FIRST_TF/logs/" 2>/dev/null | head -10
echo ""
echo "--- salvage/ subdir if present ---"
ls "openfe/production/$FIRST/$FIRST_TF/salvage/" 2>/dev/null | head -10

echo ""
echo "============================================================"
echo "2. FILE COUNTS"
echo "============================================================"
echo "all .log under production:   $(find openfe/production -name '*.log' 2>/dev/null | wc -l)"
echo "all .out under production:   $(find openfe/production -name '*.out' 2>/dev/null | wc -l)"
echo "all .err under production:   $(find openfe/production -name '*.err' 2>/dev/null | wc -l)"
echo "salvage .log:                $(find openfe/production -path '*salvage*' -name '*.log' 2>/dev/null | wc -l)"
echo "non-salvage .log:            $(find openfe/production -name '*.log' -not -path '*salvage*' 2>/dev/null | wc -l)"

echo ""
echo "--- naming pattern of log files (sample 5 full paths) ---"
find openfe/production -name "*.log" 2>/dev/null | head -5

echo ""
echo "============================================================"
echo "3. SAMPLE .log CONTENT (one complete file)"
echo "============================================================"
SAMPLE=$(find openfe/production -name "*.log" -not -path '*salvage*' 2>/dev/null | head -1)
echo "Sample file: $SAMPLE"
echo "--- full content (first 120 lines) ---"
head -120 "$SAMPLE" 2>/dev/null

echo ""
echo "============================================================"
echo "4. SLOT-IDENTIFYING FIELDS ACROSS MANY LOGS"
echo "============================================================"
echo "These fields may distinguish OSPool / prioritized / shared / backfill."
echo "Showing which appear and example values."
echo ""
for field in "SlotName" "GLIDEIN_Site" "GLIDEIN_ResourceName" "IsBackfill" \
             "ChtcProjects" "PrioritizedProjects" "RemoteHost" "MachineAttrGLIDEIN_Site0" \
             "Backfill" "OSG" "ospool" "GLIDEIN" "DeviceName" "Capability"; do
    echo "--- field: $field ---"
    grep -rh "$field" openfe/production --include="*.log" 2>/dev/null \
        | grep -i "$field" | sort | uniq -c | sort -rn | head -8
    echo ""
done

echo ""
echo "============================================================"
echo "5. DISTINCT EXECUTE HOSTS (from SlotName)"
echo "============================================================"
grep -rh "SlotName:" openfe/production --include="*.log" 2>/dev/null \
    | sed -E 's/.*@([^ ]+).*/\1/' | sort | uniq -c | sort -rn | head -40

echo ""
echo "============================================================"
echo "6. EXECUTE-EVENT COUNT PER JOB (multi-attempt check)"
echo "============================================================"
echo "Count 'Job executing' events per log to see retries/multi-attempt."
echo "Distribution of execute-events-per-log:"
for f in $(find openfe/production -name "*.log" -not -path '*salvage*' 2>/dev/null | head -200); do
    grep -c "Job executing on" "$f" 2>/dev/null
done | sort | uniq -c | sort -rn
echo "(left column = how many logs; right column = execute events in that log)"

echo ""
echo "--- exact phrasing of execute/terminate/evict events in sample ---"
grep -E "executing|terminated|evicted|aborted|Job was|return value" "$SAMPLE" 2>/dev/null | head -15
