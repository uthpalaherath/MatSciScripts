#!/usr/bin/sh

# Check if VASP relaxation is obtained for batch jobs when relaxed with
# Convergence.py and relax.dat is created.
if [[ "$*" == "" ]]; then
    arg="^[0-9]+$"
else
    arg=$1
fi

rm -f unrelaxed_list.dat
folder_list=$(ls | grep -E $arg)
for i in $folder_list;
do if [[ -f "${i}/relax.dat" ]]; then
    echo $i
else
    printf "${i}\t" >> unrelaxed_list.dat
fi
done
