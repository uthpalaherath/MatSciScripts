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

for i in */
do
    if [ -f "$i/aims.out" ]
    then
        done_word=$( tail -n 2 $i/aims.out | head -n 1 |awk '{print $1}')

        if [ "$done_word" == 'Have' ]
        then
            echo ${i:0:-1} ": Complete"
            complete_counter=$(( complete_counter + 1 ))
        else
            if [ -s "$i/aims.err" ] || [ -s "$i/aims.error" ]; then
                echo ${i:0:-1} ": Fail"
                fail_counter=$(( fail_counter + 1 ))
            else
                echo ${i:0:-1} ": Incomplete"
                incomplete_counter=$(( incomplete_counter + 1 ))
            fi
        fi
    else
        echo ${i:0:-1} ": Calculation not started"
        notstarted_counter=$(( notstarted_counter + 1 ))
    fi
done

echo "--------------------------"
echo "Complete : " $complete_counter
echo "Incomplete : " $incomplete_counter
echo "Fail : " $fail_counter
echo "Not started : " $notstarted_counter
