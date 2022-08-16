#!/bin/bash

# This script checks if the FHI-aims calculations in the directories
# are complete and outputs information regarding which ones are complete
# and which ones aren't.

# usage:
# done_counter.sh <initial POSCAR> <final POSCAR> <Total number of images>
#
# - Uthpala Herath

# Checks for aims.out in all directories here.
for i in */
do
    if [ -f "$i/aims.out" ]
    then
        done_word=$( tail -n 2 $i/aims.out | head -n 1 |awk '{print $1}')

        if [ "$done_word" == 'Have' ]
        then
            echo ${i:0:-1} ": Complete"
        else
            echo ${i:0:-1} ": Incomplete"
        fi
    else
        echo ${i:0:-1} ": Calculation not started"
    fi
done

