#!/usr/bin/env python3
"""
VASP convergence checker.

Author: Uthpala Herath

This tiny script uses pymatgen (github.com/materialsproject/pymatgen)
to read the vasp.xml file in a directory to check if the DFT calculation
has reached electronic and/or ionic convergence.

Place the script somewhere that is accessed globally.

Usage:

$ converged.sh

"""

import warnings

import pymatgen.io.vasp as vasp

warnings.filterwarnings("ignore")
vasprun = vasp.Vasprun("vasprun.xml")

print("\nElectronic convergence:%s " % (vasprun.converged_electronic))
print("Ionic convergence:%s " % (vasprun.converged_ionic))
