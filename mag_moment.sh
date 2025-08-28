#!/usr/bin/env bash
# mag_moment.sh
# Usage:
#   ./mag_moment.sh [OUTCAR] [start_index] [end_index]
# Defaults: OUTCAR, 7, 12

file=${1:-OUTCAR}
start=${2:-7}
end=${3:-12}

awk -v start="$start" -v end="$end" '
/^[[:space:]]*magnetization \(x\)/ {
    # Start of a new block; reset for the (potentially) final one
    inBlock=1; inTable=0; sum=0; count=0; next
}
inBlock && /^-+$/ {
    # First dashed line -> table begins; second dashed line -> table ends
    if (!inTable) { inTable=1; next }
    else { inTable=0; inBlock=0; next }
}
inTable && /^[[:space:]]*[0-9]+/ {
    # Data rows: index in $1; "tot" is last field ($NF)
    idx=$1; tot=$NF
    if (idx>=start && idx<=end) { sum+=tot; count++ }
    next
}
END {
    if (count>0)
        printf("Average tot moment for atoms %d-%d: %.3f (sum=%.3f, n=%d)\n",
               start, end, sum/count, sum, count);
    else
        print "Failed to find a magnetization table or matching indices."
}
' "$file"
