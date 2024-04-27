#!/usr/bin/sh

# Clean VASP files in current directoy and subdirectories.
# For only current directory use cleanvasp.sh
# Add -delete flag to delete.
find . \( \
    -name "CHGCAR*" -o \
    -name "OUTCAR*" -o \
    -name "CHG" -o \
    -name "DOSCAR" -o \
    -name "EIGENVAL" -o \
    -name "ENERGY" -o \
    -name "IBZKPT" -o \
    -name "OSZICAR*" -o \
    -name "PCDAT" -o \
    -name "REPORT" -o \
    -name "TIMEINFO" -o \
    -name "WAVECAR" -o \
    -name "XDATCAR" -o \
    -name "wannier90.wout" -o \
    -name "wannier90.amn" -o \
    -name "wannier90.mmn" -o \
    -name "wannier90.eig" -o \
    -name "wannier90.chk" -o \
    -name "wannier90.node*" -o \
    -name "PROCAR" -o \
    -name "*.o[0-9]*" -o \
    -name "vasprun.xml" -o \
    -name "relax.dat" -o \
    -name "CONTCAR*" \
\) -type f $1
