#!/bin/bash

# This script checks if the FHI-aims calculations in the directories
# are complete and outputs information regarding which ones are complete
# and which ones aren't.

# usage:
# done_counter.sh <initial POSCAR> <final POSCAR> <Total number of images>
#
# - Uthpala Herath

# Checks for aims.out in all directories here.
complete_counter=0
incomplete_counter=0
fail_counter=0
notstarted_counter=0

for i in */; do
    if [ -f "$i/aims.out" ]; then
        if grep -q "Have a nice day." "$i/aims.out" 2>/dev/null; then
            echo "${i%/} : Complete"
            complete_counter=$((complete_counter + 1))
        else
            # Define error keywords
            error_patterns="Error|CANCELLED|insufficient virtual memory|Exited"
            if grep -Eiq "$error_patterns" "$i/aims.err" 2>/dev/null; then
                echo "${i%/} : Fail"
                fail_counter=$((fail_counter + 1))
            else
                echo "${i%/} : Incomplete"
                incomplete_counter=$((incomplete_counter + 1))
            fi
        fi
    else
        echo "${i%/} : Calculation not started"
        notstarted_counter=$((notstarted_counter + 1))
    fi
done

echo "--------------------------"
echo "Complete : $complete_counter"
echo "Incomplete : $incomplete_counter"
echo "Fail : $fail_counter"
echo "Not started : $notstarted_counter"
