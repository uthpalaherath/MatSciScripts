#!/usr/bin/env python3

"""
VASP Convergence Test.

Author: Uthpala Herath

This script is based on Pedram Tavazohi's tutorial on PyChemia
(github.com/MaterialsDiscovery/PyChemia.git) functions to perform k-grid and
energy cut-off convergence. It is assumed that the structure file is named POSCAR
and the pseudopotentials are stored in ~/.vasp/PP-VASP/{potpaw_PBE,potpaw_LDA}/.
The VASP executable is assumed to be vasp_std.

Usage:

$ ConvergenceTest.py {kgrid,encut,complete} -np <number of processors> -extra_vars '{"key" : "value"}' -pspdir {potpaw_PBE,potpaw_LDA} -update

E.g.-

$ ConvergenceTest.py encut -np 16 -extra_vars '{"LVCADER" : ".TRUE.", "VCA" : "0.3 0.7 1.0 1.0"}' -pspdir potpaw_PBE -update

"""

import argparse
import json
import os
import re
import sys
from argparse import RawTextHelpFormatter

import pychemia


def load_poscar():
    """
    Returns the structure retrieved from POSCAR.
    Also creates backups of INCAR and KPOINTS.
    """

    if os.path.exists("INCAR"):
        os.rename("INCAR", "INCAR.bak")
    else:
        print("INCAR not found! Generating with PyChemia.")

    if os.path.exists("KPOINTS"):
        os.rename("KPOINTS", "KPOINTS.bak")
    else:
        print("KPOINTS not found! Generating with PyChemia.")

    return pychemia.code.vasp.read_poscar("POSCAR")


def kgrid(args):
    """
    Function for k-grid convergence.

    """
    print("\nRunning k-grid convergence...")
    st = load_poscar()
    kpt_conv = pychemia.code.vasp.task.ConvergenceKPointGrid(
        structure=st,
        workdir=".",
        executable="vasp_std",
        pspdir=args.pspdir,
        extra_vars=args.extra_vars,
    )
    kpt_conv.run(args.np)
    print("\nk-grid coverged: ", kpt_conv.success)

    if kpt_conv.success == True:
        print("Optimal k-grid: ", kpt_conv.best_kpoints.grid)
    else:
        print("k-grid convergence failed!")
    os.rename("INCAR.bak", "INCAR")
    os.rename("KPOINTS.bak", "KPOINTS")

    # Create new  KPOINTS file
    if args.update:
        if kpt_conv.success == True:
            f = open("KPOINTS", "w")
            f.write("Automatic mesh\n")
            f.write("0\n")
            f.write("Gamma\n")
            f.write(
                "%d %d %d\n"
                % (
                    kpt_conv.best_kpoints.grid[0],
                    kpt_conv.best_kpoints.grid[1],
                    kpt_conv.best_kpoints.grid[2],
                )
            )
            f.write("0 0 0")
            f.close()
        else:
            print("Update failed!")


def encut(args):
    """
    Function for ENCUT convergence.

    """
    print("\nRunning ENCUT convergence...")
    st = load_poscar()
    encut_conv = pychemia.code.vasp.task.ConvergenceCutOffEnergy(
        structure=st,
        workdir=".",
        executable="vasp_std",
        pspdir=args.pspdir,
        extra_vars=args.extra_vars,
        energy_tolerance=1e-3,
    )
    encut_conv.run(args.np)
    print("\nENCUT coverged: ", encut_conv.success)

    if encut_conv.success == True:
        print("Optimal ENCUT: ", encut_conv.best_encut)
    else:
        print("ENCUT convergence failed!")
    os.rename("INCAR.bak", "INCAR")
    os.rename("KPOINTS.bak", "KPOINTS")

    if args.update:
        if encut_conv.success == True:
            # Update INCAR with new ENCUT
            estr = "ENCUT = " + encut_conv.best_encut
            with open("INCAR", "r") as sources:
                lines = sources.readlines()
            with open("INCAR", "w") as sources:
                for line in lines:
                    sources.write(re.sub(r"ENCUT\s*=\s*([\d.]*)", estr, line))
        else:
            print("Update failed!")


def complete(args):
    """
    Function for both k-grid and ENCUT convergence.
    """
    kgrid(args)
    encut(args)


if __name__ == "__main__":

    args = sys.argv[1:]
    if args:

        # top level parser
        parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=RawTextHelpFormatter
        )
        subparsers = parser.add_subparsers(help="sub-command help")

        # parser for k-grid
        parser_kgrid = subparsers.add_parser("kgrid", help="k-grid convergence")
        parser_kgrid.add_argument(
            "-np", default=1, type=int, help="Number of processors",
        )
        parser_kgrid.add_argument(
            "-extra_vars", default=None, type=json.loads, help="Extra INCAR parameters"
        )
        parser_kgrid.add_argument(
            "-pspdir",
            default="potpaw_PBE",
            type=str,
            help="Pseudopotential",
            choices=["potpaw_LDA", "potpaw_PBE"],
        )
        parser_kgrid.add_argument(
            "-update",
            help="Write new KPOINTS file? (Will use a Gamma centered mesh.)",
            action="store_true",
        )

        parser_kgrid.set_defaults(func=kgrid)

        # parser for encut
        parser_encut = subparsers.add_parser("encut", help="Energy cut-off convergence")
        parser_encut.add_argument(
            "-np", default=1, type=int, help="Number of processors",
        )
        parser_encut.add_argument(
            "-extra_vars", default=None, type=json.loads, help="Extra INCAR parameters"
        )
        parser_encut.add_argument(
            "-pspdir",
            default="potpaw_PBE",
            type=str,
            help="Pseudopotential",
            choices=["potpaw_LDA", "potpaw_PBE"],
        )
        parser_encut.add_argument(
            "-update", help="Update INCAR with new ENCUT?", action="store_true",
        )
        parser_encut.set_defaults(func=encut)

        # parser for both k-grid and encut
        parser_complete = subparsers.add_parser("complete", help="Complete convergence")
        parser_complete.add_argument(
            "-np", default=1, type=int, help="Number of processors",
        )
        parser_complete.add_argument(
            "-extra_vars", default=None, type=json.loads, help="Extra INCAR parameters"
        )
        parser_complete.add_argument(
            "-pspdir",
            default="potpaw_PBE",
            type=str,
            help="Pseudopotential",
            choices=["potpaw_LDA", "potpaw_PBE"],
        )
        parser_complete.add_argument(
            "-update", help="Update INCAR and KPOINTS?", action="store_true",
        )
        parser_complete.set_defaults(func=complete)

        # End of sub-parsers
        args = parser.parse_args()
        args.func(args)
    else:
        print("Usage: convergencetest.py [-h]")
        print("{kgrid, encut, complete}")
