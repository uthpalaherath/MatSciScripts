#!/usr/bin/env python

"""
fdf to POSCAR converter.

This script creates a POSCAR from a Siesta .fdf file.

Usage: fdftoposcar.py <name>.fdf

- Uthpala Herath

"""

import re
import sys
from itertools import groupby

file = open(sys.argv[1], "r")
data = file.read()
file.close()

atoms = []
lattice_constant = float(re.findall(r"LatticeConstant\s*([\d.]*)", data)[0])
lattice_vectors = re.findall(r"LatticeVectors\s*([\d.\s]*)%endblock", data)
atomic_coordinates = re.findall(
    r"AtomicCoordinatesAndAtomicSpecies\s*([\d.\sa-zA-Z]*)%endblock", data
)
atomic_coordinates_lines = atomic_coordinates[0].split("\n")
atm_coord_len = len(atomic_coordinates_lines)
for i in range(atm_coord_len - 1):
    atoms.append(atomic_coordinates_lines[i].split()[-1])
species = [i[0] for i in groupby(atoms)]
species_count = [len(list(group)) for key, group in groupby(atoms)]

f = open("POSCAR", "w")
f.write(" ".join(str(x) for x in species))
f.write("\n%f\n" % lattice_constant)
f.write("\n".join(str(x) for x in lattice_vectors[0].split("\n")))
f.write(" ".join(str(x) for x in species))
f.write("\n")
f.write(" ".join(str(x) for x in species_count))
f.write("\nDirect\n")
for i in range(atm_coord_len - 1):
    f.write(
        " ".join(
            [
                " ".join(map(str, atomic_coordinates_lines[i].split()[0:3])),
                atomic_coordinates_lines[i].split()[-1],
                "\n",
            ]
        )
    )
f.close()
