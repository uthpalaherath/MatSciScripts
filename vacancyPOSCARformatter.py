#!/usr/bin/env python
""" Vacancy POSCAR Formatter.

This script takes a POSCAR created with sod and reformats it
to create the POSCAR for the vacancy structure

Usage:
vacancyPOSCARformatter.py <infile> <outfile>


"""

import argparse

import numpy as np


def formatter(infile="POSCAR", outfile="POSCAR_new"):
    """
    This method takes in an unformatted poscar that includes a vacancy and reformats it.
    """

    file = open(infile, "r")
    POSCAR = file.readlines()

    # scale factor
    scale_factor = float(POSCAR[1])

    # coordinate type
    co_type = POSCAR[7]

    # cell
    cell_matrix = POSCAR[2:5]
    cell = np.zeros(shape=(3, 3))

    for i in range(len(cell_matrix)):
        cell_matrix0 = np.array(cell_matrix[i].split())
        cell[i, :] = (cell_matrix0.astype(np.float)) * np.array(
            POSCAR[1].split()
        ).astype(np.float)

    # positions
    atoms = np.array(POSCAR[6].split()).astype(np.int)
    positions_matrix = POSCAR[8 : 8 + sum(atoms)]
    positions = np.zeros(shape=(np.sum(atoms), 3))

    for j in range(len(positions_matrix)):
        positions_matrix0 = np.array(positions_matrix[j].split())[0:3]
        positions[j, :] = positions_matrix0.astype(np.float)

    structure_name_list = POSCAR[5].split()
    atom_type_list = POSCAR[6].split()

    # Where is X in structure_name_list
    X_index = structure_name_list.index("X")

    # number of vacancies
    num_vac = int(atom_type_list[X_index])

    # atoms before vacancy
    num_prev_atoms = sum(np.array(atom_type_list[0:X_index], dtype=int))

    # reformatted array
    new_positions = np.concatenate(
        (positions[0:num_prev_atoms], positions[num_prev_atoms + num_vac :]), axis=0
    )

    # updated structure_name_list by deleting the 'X' in the structure_name_list
    del structure_name_list[X_index]

    # updated atom_type_list by deleting the index corresponding to 'X' in the structure_name_list
    del atom_type_list[X_index]

    print("Number of vacancies = %s" % num_vac)
    # writing new POSCAR
    p_file = open(outfile, "w+")
    p_file.write("POSCAR for %d vacancies\n" % num_vac)
    p_file.write("%f \n" % scale_factor)
    for i_cell in range(len(cell)):
        p_file.write(
            "%f %f %f \n" % (cell[i_cell, 0], cell[i_cell, 1], cell[i_cell, 2])
        )
    for i_snl in range(len(structure_name_list)):
        p_file.write("%s " % structure_name_list[i_snl])
    p_file.write("\n")
    for i_atl in range(len(atom_type_list)):
        p_file.write("%d " % int(atom_type_list[i_atl]))
    p_file.write("\n")
    p_file.write("%s" % co_type)
    for i_pos in range(len(new_positions)):
        p_file.write(
            "%f %f %f \n"
            % (
                new_positions[i_pos, 0],
                new_positions[i_pos, 1],
                new_positions[i_pos, 2],
            )
        )

    p_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vacancy POSCAR Formatter")
    parser.add_argument("infile", type=str, help="Input POSCAR name")
    parser.add_argument("outfile", type=str, help="Formatted output POSCAR name")
    args = parser.parse_args()
    formatter(args.infile, args.outfile)
