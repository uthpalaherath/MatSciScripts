#!/usr/bin/env python3

"""
POTCAR generator.

Author: Uthpala Herath

This script uses PyChemia (github.com/MaterialsDiscovery/PyChemia.git)
to generate a POTCAR file from a given POSCAR.
It is assumed that the structure file is named POSCAR (otherwise use -poscar)
and the pseudopotentials are stored in ~/.vasp/PP-VASP/{potpaw_PBE,potpaw_LDA}/.
If the repeating order of atoms in the POSCAR is to be kept use the
-heterostructure flag.

Usage:

$ makePOTCAR.py
        -pspdir {potpaw_PBE,potpaw_LDA}
        -psp_options '{"key" : "value"}'
        -poscar POSCAR_SrVO3
        -heterostructure
E.g.-

$ makePOTCAR.py -pspdir potpaw_PBE -psp_options '{"Sr":"sv"}'

"""

import argparse
import json
import sys
from argparse import RawTextHelpFormatter
from itertools import groupby

import pychemia


def load_poscar(poscarname):
    """
    Returns the structure retrieved from POSCAR.
    """
    return pychemia.code.vasp.read_poscar(poscarname)


def makePOTCAR(args):
    """
    Function to create POTCAR.
    """
    st = load_poscar(args.poscar)

    pychemia.code.vasp.poscar.write_potcar(
        structure=st,
        pspdir=args.pspdir,
        options=args.psp_options,
        heterostructure=args.heterostructure,
    )
    if args.heterostructure:
        print("POTCAR generated for: ", [i[0] for i in groupby(st.symbols)])
    else:
        print("POTCAR generated for: ", pychemia.code.vasp.poscar.get_species_list(st))


if __name__ == "__main__":

    # top level parser
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "-pspdir",
        default="potpaw_PBE",
        type=str,
        help="Pseudopotential.",
        choices=["potpaw_LDA", "potpaw_PBE"],
    )

    parser.add_argument(
        "-psp_options", help="Pseudopotential options. ", default=None, type=json.loads,
    )

    parser.add_argument(
        "-poscar", default="POSCAR", type=str, help="POSCAR file.",
    )

    parser.add_argument(
        "-heterostructure",
        action="store_true",
        help="Keep repeating order of atoms in POSCAR?",
    )

    # End of sub-parsers
    args = parser.parse_args()
    makePOTCAR(args)
