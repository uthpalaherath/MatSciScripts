#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Structure file converter.

This coverts DFT structures and parses structural
information.

- Uthpala Herath

"""


import argparse
import os
import re
import shutil
import sys
from itertools import groupby
from shutil import copyfile

import numpy as np


class Converters:
    """This class contains methods for converting and
    parsing files.
    """

    def __init__(self):
        self.structurename = None

    def fdf_to_poscar(self):
        """
	This function converts the siesta .fdf format to POSCAR for further calculations.
	"""
        # file = pychemia.code.siesta.SiestaInput(self.structurename + ".fdf")
        # self.st = file.get_structure()
        # pychemia.code.vasp.write_poscar(
        #     self.st, filepath="POSCAR", newformat=True, direct=True, comment=None
        # )
        fname = self.structurename + ".fdf"
        file = open(fname, "r")
        data = file.read()
        file.close()

        atoms = []
        lattice_constant = float(re.findall(r"LatticeConstant\s*([\d.]*)", data)[0])
        lattice_vectors = re.findall(r"LatticeVectors\s*([\d.\s]*)%endblock", data)

        # creating a numpy array with lattice vectors
        lattice_vec = np.array(lattice_vectors[0].split(), dtype="float64")
        lattice_vec = lattice_vec.reshape(3, 3)
        self.cell = lattice_constant * lattice_vec

        atomic_coordinates = re.findall(
            r"AtomicCoordinatesAndAtomicSpecies\s*([\d.\sa-zA-Z]*)%endblock", data
        )
        atomic_coordinates_lines = atomic_coordinates[0].split("\n")

        atm_coord_len = len(atomic_coordinates_lines)
        for i in range(atm_coord_len - 1):
            atoms.append(atomic_coordinates_lines[i].split()[-1])

        self.symbols = atoms

        species = [i[0] for i in groupby(atoms)]
        species_count = [len(list(group)) for key, group in groupby(atoms)]

        self.positions = np.zeros((atm_coord_len - 1, 3), dtype="float64")
        for counter in range(atm_coord_len - 1):
            self.positions[counter, :] = atomic_coordinates_lines[counter].split()[0:3]

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

    def qe_to_poscar(self):
        """Creates a POSCAR from a Quantum Espressso scf input file.
           CELL_PARAMETERS must be defined.
        """
        fname = self.structurename + ".scf.in"
        file = open(fname, "r")
        data = file.read()
        file.close()
        lattice_vec = re.findall(r"CELL_PARAMETERS\s*[a-zA-Z]*([e\d\s.+-]*)", data)
        lattice_vec = [float(x) for x in lattice_vec[0].split()]
        self.cell = np.array((lattice_vec), dtype="float64").reshape(3, 3)
        nat = int(re.findall(r"nat\s*=\s*([\d]*)", data)[0])
        raw_positions = re.findall(
            r"ATOMIC_POSITIONS\s*crystal([e\d.+\sa-zA-Z]*)\n", data
        )
        full_structure = np.array(raw_positions[0].split()).reshape(nat, 4)
        self.positions = np.zeros((nat, 3), dtype="float64")
        self.symbols = []
        for icount, i in enumerate(full_structure):
            self.positions[icount] = i[1:]
            self.symbols.append(i[0])
        species = [i[0] for i in groupby(self.symbols)]
        species_count = [len(list(group)) for key, group in groupby(self.symbols)]

        # Writing to POSCAR
        f = open("POSCAR", "w")
        f.write(" ".join(str(x) for x in species))
        f.write("\n%f\n" % 1.0)
        for i in range(len(self.cell)):
            f.write("%f %f %f\n" % (self.cell[i, 0], self.cell[i, 1], self.cell[i, 2]))
        f.write(" ".join(str(x) for x in species))
        f.write("\n")
        f.write(" ".join(str(x) for x in species_count))
        f.write("\n")
        f.write("Direct\n")

        for i in range(len(full_structure)):
            f.write(
                "%s %s %s %s\n"
                % (
                    full_structure[i, 1],
                    full_structure[i, 2],
                    full_structure[i, 3],
                    full_structure[i, 0],
                )
            )
        f.close()

    def read_poscar(self, fname="POSCAR"):
        "Reads a POSCAR and sets self.cell, self.symbols and self.positions."

        file = open(fname, "r")
        data = file.readlines()
        file.close()

        lattice_constant = float(data[1])
        lattice_vec = np.array(
            (
                [float(x) for x in data[2].split()],
                [float(x) for x in data[3].split()],
                [float(x) for x in data[4].split()],
            )
        )
        self.cell = lattice_constant * lattice_vec
        num_atoms = np.sum([int(x) for x in data[6].split()])
        full_structure = np.array((data[8 : 8 + num_atoms]), dtype="str")
        self.positions = np.zeros((num_atoms, 3), dtype="float64")
        self.symbols = []
        for icount, i in enumerate(full_structure):
            self.positions[icount] = i.split()[0:3]
            self.symbols.append(i.split()[-1])

    def win_to_poscar(self, input="aiida.win"):
        """This function generates a POSCAR from a .win file.
        Used in aiida calculations where the .win file is provided."""

        file = open(input, "r")
        data = file.read()
        file.close()

        # Unit cell
        unit_cell_data = re.findall(
            r"begin\s*unit\s*_cell_cart\s*[a-zA-Z]*([+-e\s\d.]*)end", data
        )
        unit_cell = np.array(unit_cell_data[0].split(), dtype="float64")
        unit_cell = unit_cell.reshape(3, 3)

        # coordinate type (cartesian or direct)
        coord_type = re.findall(r"begin\s*atoms_([a-z]*)", data)[0]

        # atomic coordinates
        atom_coords = np.array(
            re.findall(r"begin\s*atoms_[\sa-z]*([\sa-zA-Z\d.]*)end", data)[0].split()
        )
        atom_coords = atom_coords.reshape(int(len(atom_coords) / 4), 4)

        atoms = atom_coords[:, 0]
        species = [i[0] for i in groupby(atoms)]
        species_count = [len(list(group)) for key, group in groupby(atoms)]

        f = open("POSCAR", "w")
        f.write(" ".join(str(x) for x in species))
        f.write("\n%f\n" % 1.0)
        for i in range(len(unit_cell)):
            f.write("%f %f %f\n" % (unit_cell[i, 0], unit_cell[i, 1], unit_cell[i, 2]))
        f.write(" ".join(str(x) for x in species))
        f.write("\n")
        f.write(" ".join(str(x) for x in species_count))
        f.write("\n")
        if coord_type == "cart":
            f.write("Cartesian\n")
        elif coord_type == "frac":
            f.write("Direct\n")

        for i in range(len(atom_coords)):
            f.write(
                "%s %s %s %s\n"
                % (
                    atom_coords[i, 1],
                    atom_coords[i, 2],
                    atom_coords[i, 3],
                    atom_coords[i, 0],
                )
            )

        f.close()
