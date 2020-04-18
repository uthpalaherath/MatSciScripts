#!/usr/bin/env python3

"""
VASP Convergence Test.

Author: Uthpala Herath

This script is based on Pedram Tavazohi's tutorial on PyChemia
(github.com/MaterialsDiscovery/PyChemia.git) functions to perform k-grid and
energy cut-off convergence. It is assumed that the structure file is named POSCAR
and the pseudopotentials are stored in ~/.vasp/PP-VASP/{potpaw_PBE,potpaw_LDA}/.
The VASP executable is assumed to be vasp_std.
Perform the ionic relaxation only after both ENCUT and k-grid convergences are done.

Usage:

$ ConvergenceTest.py {kgrid,encut,complete,relax} -np <number of processors> -extra_vars '{"key" : "value"}' -pspdir {potpaw_PBE,potpaw_LDA} -update -psp_options '{"key" : "value"}'

E.g.-

$ ConvergenceTest.py encut -np 16 -extra_vars '{"LVCADER" : ".TRUE.", "VCA" : "0.3 0.7 1.0 1.0"}' -pspdir potpaw_PBE -update -psp_options '{"Sr":"sv"}'

$ ConvergenceTest.py relax -np 16 -max_calls 30

"""

import argparse
import json
import os
import re
import shutil
import sys
from argparse import RawTextHelpFormatter

import pychemia


def restore(type):
    """
    Restores created backups of INCAR, KPOINTS and POSCAR.
    """

    if type == "kgrid":
        if os.path.exists("KPOINTS.bak"):
            shutil.copy("KPOINTS.bak", "KPOINTS")
        else:
            print(
                "KPOINTS.bak not found! KPOINTS was probably generated with PyChemia. "
            )
        if os.path.exists("INCAR.bak"):
            shutil.copy("INCAR.bak", "INCAR")
        else:
            print("INCAR.bak not found! INCAR was probably generated with PyChemia. ")

    elif type == "encut":
        if os.path.exists("INCAR.bak"):
            shutil.copy("INCAR.bak", "INCAR")
        else:
            print("INCAR.bak not found! INCAR was probably generated with PyChemia. ")


def backup(type):
    """
    Creates backups of INCAR, KPOINTS and POSCAR.
    """
    if type == "kgrid":
        if os.path.exists("KPOINTS"):
            shutil.copy("KPOINTS", "KPOINTS.bak")
        else:
            print("KPOINTS not found! Generating with PyChemia.")
        if os.path.exists("INCAR"):
            shutil.copy("INCAR", "INCAR.bak")
        else:
            print("INCAR not found! Generating with PyChemia.")

    elif type == "encut":
        if os.path.exists("INCAR"):
            shutil.copy("INCAR", "INCAR.bak")
        else:
            print("INCAR not found! Generating with PyChemia.")

    elif type == "relax":
        if os.path.exists("POSCAR"):
            shutil.copy("POSCAR", "POSCAR.bak")
        else:
            print("POSCAR not found!")
            sys.exit()

        if os.path.exists("INCAR"):
            shutil.copy("INCAR", "INCAR.bak")
        else:
            print("INCAR not found! Generating with PyChemia.")


def load_poscar():
    """
    Returns the structure retrieved from POSCAR.
    """
    return pychemia.code.vasp.read_poscar("POSCAR")


def kgrid(args):
    """
    Function for k-grid convergence.
    """

    print("\nRunning k-grid convergence...")
    st = load_poscar()
    backup("kgrid")

    kpt_conv = pychemia.code.vasp.task.ConvergenceKPointGrid(
        structure=st,
        workdir=".",
        executable="vasp_std",
        pspdir=args.pspdir,
        extra_vars=args.extra_vars,
        psp_options=args.psp_options,
        energy_tolerance=args.energy_tolerance,
    )
    kpt_conv.run(args.np)
    print("\nk-grid coverged: ", kpt_conv.success)

    if kpt_conv.success == True:
        print("Optimal k-grid: ", kpt_conv.best_kpoints.grid)
    else:
        print("k-grid convergence failed!")

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
    else:
        # restoring files
        restore("kgrid")


def encut(args):
    """
    Function for ENCUT convergence.

    """
    print("\nRunning ENCUT convergence...")
    st = load_poscar()
    backup("encut")

    encut_conv = pychemia.code.vasp.task.ConvergenceCutOffEnergy(
        structure=st,
        workdir=".",
        executable="vasp_std",
        pspdir=args.pspdir,
        extra_vars=args.extra_vars,
        energy_tolerance=args.energy_tolerance,
        psp_options=args.psp_options,
    )

    encut_conv.run(args.np)
    print("\nENCUT coverged: ", encut_conv.success)

    if encut_conv.success == True:
        print("Optimal ENCUT: ", encut_conv.best_encut)
    else:
        print("ENCUT convergence failed!")

    if args.update:
        # restoring files
        restore("encut")
        if encut_conv.success == True:
            # Update INCAR with new ENCUT
            estr = "ENCUT = " + str(encut_conv.best_encut)
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


def relax(args):
    """
    Function to relax structure after kgrid and encut convergence is complete.
    """

    print("\nRunning ionic relaxation...")

    # retrieve ENCUT from INCAR
    fi = open("INCAR", "r")
    data = fi.read()
    fi.close()
    encut_val = float(re.findall(r"ENCUT\s*=\s*([\d.]*)", data)[0].split()[0])

    # retrieve kgrid from KPOINTS
    fi = open("KPOINTS", "r")
    for i in range(3):
        fi.readline()
    gridline = fi.readline()
    fi.close()
    gridline = [int(x) for x in gridline.split()]

    st = load_poscar()
    backup("relax")
    relax_st = pychemia.code.vasp.task.IonRelaxation(
        structure=st,
        workdir=".",
        target_forces=args.target_forces,
        executable="vasp_std",
        encut=encut_val,
        kp_grid=gridline,
        pspdir=args.pspdir,
        max_calls=args.max_calls,
        extra_vars=args.extra_vars,
        # energy_tolerance=args.energy_tolerance, #Not an input argument
    )
    relax_st.run(args.np)

    # Test if relaxed
    ediffg = abs(pychemia.code.vasp.VaspInput("INCAR").EDIFFG)
    vaspout = pychemia.code.vasp.VaspOutput("OUTCAR")
    avg_force = vaspout.relaxation_info()["avg_force"]
    print("EDIFFG = %f and Average force = %f" % (ediffg, avg_force))

    if avg_force <= ediffg:
        print("Forces are converged.")
    else:
        print("Forces are not converged.")


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
            "-extra_vars", default=None, type=json.loads, help="Extra INCAR parameters."
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
        parser_kgrid.add_argument(
            "-psp_options",
            help="Pseudopotential options. ",
            default=None,
            type=json.loads,
        )
        parser_kgrid.add_argument(
            "-energy_tolerance",
            default=1e-4,
            type=float,
            help="The energy difference required for convergence per atom",
        )

        parser_kgrid.set_defaults(func=kgrid)

        # parser for encut
        parser_encut = subparsers.add_parser("encut", help="Energy cut-off convergence")
        parser_encut.add_argument(
            "-np", default=1, type=int, help="Number of processors",
        )
        parser_encut.add_argument(
            "-extra_vars", default=None, type=json.loads, help="Extra INCAR parameters."
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
        parser_encut.add_argument(
            "-psp_options",
            help="Pseudopotential options. ",
            default=None,
            type=json.loads,
        )
        parser_encut.add_argument(
            "-energy_tolerance",
            default=1e-4,
            type=float,
            help="The energy difference required for convergence per atom",
        )
        parser_encut.set_defaults(func=encut)

        # parser for both k-grid and encut
        parser_complete = subparsers.add_parser("complete", help="Complete convergence")
        parser_complete.add_argument(
            "-np", default=1, type=int, help="Number of processors",
        )
        parser_complete.add_argument(
            "-extra_vars", default=None, type=json.loads, help="Extra INCAR parameters."
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
        parser_complete.add_argument(
            "-psp_options",
            help="Pseudopotential options. ",
            default=None,
            type=json.loads,
        )
        parser_complete.add_argument(
            "-energy_tolerance",
            default=1e-4,
            type=float,
            help="The energy difference required for convergence per atom",
        )
        parser_complete.set_defaults(func=complete)

        # parser for ionic relaxation
        parser_relax = subparsers.add_parser("relax", help="Ionic relaxation")
        parser_relax.add_argument(
            "-np", default=1, type=int, help="Number of processors",
        )
        parser_relax.add_argument(
            "-extra_vars", default=None, type=json.loads, help="Extra INCAR parameters."
        )
        parser_relax.add_argument(
            "-pspdir",
            default="potpaw_PBE",
            type=str,
            help="Pseudopotential",
            choices=["potpaw_LDA", "potpaw_PBE"],
        )
        parser_relax.add_argument(
            "-target_forces", default=1e-4, type=float, help="Target force difference",
        )
        parser_relax.add_argument(
            "-max_calls",
            default=20,
            type=int,
            help="Number of calls for copying CONTCAR to POSCAR",
        )
        parser_relax.add_argument(
            "-psp_options",
            help="Pseudopotential options. ",
            default=None,
            type=json.loads,
        )

        # parser_relax.add_argument(
        #     "-energy_tolerance",
        #     default=1e-8,
        #     type=float,
        #     help="The energy difference required for convergence per atom",
        # )
        parser_relax.set_defaults(func=relax)

        # End of sub-parsers
        args = parser.parse_args()
        args.func(args)
    else:
        print("Usage: convergencetest.py [-h]")
        print("{kgrid, encut, complete, relax}")
