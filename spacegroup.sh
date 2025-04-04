#!/usr/bin/env python
"""
Spacegroup.

Author: Uthpala Herath

This tiny script uses PyChemia (github.com/MaterialsDiscovery/PyChemia)
to read the POSCAR file in a directory and returns its spacegroup.

Place the script somewhere that is accessed globally.

Usage:

$ spacegroup.sh

"""
import pychemia
import sys

args = sys.argv[1:]

if args:
    poscar = sys.argv[1]
else:
    poscar = "POSCAR"

st = pychemia.code.vasp.read_poscar(poscar)
cs = pychemia.crystal.CrystalSymmetry(st)
print(cs.get_spacegroup())
