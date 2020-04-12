#!/bin/bash

# This script performs a linear interpolation of an initial POSCAR
# and an final POSCAR using nebmake.pl from VTST tools.
# Then it puts all the generated POSCARs into the directory
# 'linint'. Afterwards it runs DiSPy with provided INPUT.
# The perturbed files will be copied back to the respective directories.
# The default irreducible representation is set to #1 for perturbation.

# usage:
# makeNEB.sh <initial POSCAR> <final POSCAR> <Total number of images>
#
# - Uthpala Herath


# Images in VTST = Total number of images - 2
num_image=$(($3-2))

# Create interpolated images in directories 00..Total numer of images
nebmake.pl $1 $2 $num_image

# Copy images to linint
rm -rf linint
mkdir linint

for (( i=0; i<$3; i++ ))
do
    cp $(printf %02d $i)/POSCAR linint/$i.vasp
done

echo "Interpolated images copied into linint directory. Running dispy..."

# Running dispy without perturbing
# Unperturbed input = INPUT
# IMAGES = Total number of images
sed -i 's/PERTURB=.*/PERTURB=FALSE/' INPUT
sed -i "s/IMAGES=.*/IMAGES=$3/" INPUT
dispy INPUT

# Check results/output.out for the desired irreducible representation number.
# Set IRR_NUM in INPUT accordingly and PERTURB = TRUE.
# Automate this for the first number for now.
sed -i 's/PERTURB=FALSE/PERTURB=TRUE/g' INPUT

# running dispy for the pertubation
echo "Performing perturbation..."
dispy INPUT

# copying perturbed images back to directories
echo "Copying perturbed images back to directories ..."

for (( i=0; i<$3; i++ ))
do
    cp ./results/output_structures/$i $(printf %02d $i)/POSCAR
done
echo "Done. Ready to run NEB with VASP."

# Reset PERTURB=FALSE
sed -i 's/PERTURB=TRUE/PERTURB=FALSE/g' INPUT
rm -rf results_old*
rm -rf linint




