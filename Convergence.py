#!/usr/bin/env python3

"""
VASP Convergence testing and structural optimization.

Author: Uthpala Herath

This script is based on Pedram Tavazohi's tutorial on PyChemia
(github.com/MaterialsDiscovery/PyChemia.git) functions to perform
k-grid and energy cut-off convergence. It is assumed that the
structure file is named POSCAR and the pseudopotentials are stored
in ~/.vasp/PP-VASP/{potpaw_PBE,potpaw_LDA}/. The VASP executable is
assumed to be vasp_std. Perform the ionic relaxation only after both
ENCUT and k-grid convergences are done.
For keeping the repeating order in POSCAR when generating POTCARs
use -heterostructure flag.
To read values from INCAR use the -incar flag followed
by INCAR.
-make_potcar will automatically generate a POTCAR from POSCAR.
-auto_ibrion turns on the adaptive IBRION. If VTST Tools is available
-fire will use the FIRE algorithm for relaxation.

WARNING: INCAR and KPOINTS will be overwritten. Remember to backup originals.

Usage:

$ Convergence.py {kgrid,encut,complete,relax}
                    -np <number of processors>
                    -extra_vars '{"key" : "value"}'
                    -pspdir {potpaw_PBE,potpaw_LDA}
                    -psp_options '{"key" : "value"}'
                    -heterostructure
                    -incar INCAR
                    -make_potcar
                    -auto_ibrion
                    -fire

E.g.-

$ Convergence.py encut
                    -np 16
                    -extra_vars '{"NCORE" : "2", "ISPIN" : "1"}'
                    -pspdir potpaw_PBE
                    -psp_options '{"Sr":"sv"}'

$ Convergence.py relax
                    -np 16
                    -relax_cell
                    -incar INCAR

"""

import argparse
import json
import os
import re
import shutil
import sys
from argparse import RawTextHelpFormatter

import pychemia
from pychemia.code.vasp import incar


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

    if args.incar is None:
        extra_vars = args.extra_vars
    else:
        in_obj = incar.read_incar(args.incar)
        extra_vars = in_obj.variables

    kpt_conv = pychemia.code.vasp.task.ConvergenceKPointGrid(
        structure=st,
        workdir=".",
        executable=args.vasp_exe,
        pspdir=args.pspdir,
        extra_vars=extra_vars,
        psp_options=args.psp_options,
        energy_tolerance=args.energy_tolerance,
        heterostructure=args.heterostructure,
        make_potcar=args.make_potcar,
    )
    kpt_conv.run(args.np)
    print("\nk-grid coverged: ", kpt_conv.success)

    if kpt_conv.success is True:
        print("Optimal k-grid: ", kpt_conv.best_kpoints.grid)
        if os.path.exists("convergence.dat"):
            fi = open("convergence.dat", "a")
        else:
            fi = open("convergence.dat", "w")
        fi.write("Optimal k-grid: %s \n" % kpt_conv.best_kpoints.grid)
        fi.close()
    else:
        print("k-grid convergence failed!")

    # Create new  KPOINTS file
    if kpt_conv.success is True:
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

    return kpt_conv


def encut(args, best_kgrid=None):
    """
    Function for ENCUT convergence.

    """
    print("\nRunning ENCUT convergence...")
    st = load_poscar()

    # Remove ENCUT from INCAR
    with open("INCAR", "r") as f:
        lines = f.readlines()
    with open("INCAR", "w") as f:
        for line in lines:
            if not re.match(r"ENCUT", line):
                f.write(line)

    if args.incar is None:
        extra_vars = args.extra_vars
    else:
        in_obj = incar.read_incar(args.incar)
        extra_vars = in_obj.variables

    encut_conv = pychemia.code.vasp.task.ConvergenceCutOffEnergy(
        structure=st,
        workdir=".",
        executable=args.vasp_exe,
        pspdir=args.pspdir,
        extra_vars=extra_vars,
        energy_tolerance=args.energy_tolerance,
        psp_options=args.psp_options,
        heterostructure=args.heterostructure,
        kpoints=best_kgrid,
        increment_factor=0.1,
        make_potcar=args.make_potcar,
    )

    encut_conv.run(args.np)
    print("\nENCUT coverged: ", encut_conv.success)

    if encut_conv.success is True:
        print("Optimal ENCUT: ", encut_conv.best_encut)
        if os.path.exists("convergence.dat"):
            fi = open("convergence.dat", "a")
        else:
            fi = open("convergence.dat", "w")
        fi.write("Optimal ENCUT: %s \n" % encut_conv.best_encut)
        fi.close()
    else:
        print("ENCUT convergence failed!")

    return encut_conv


def complete(args):
    """
    Function for both k-grid and ENCUT convergence.
    """
    kpt_conv = kgrid(args)

    # Run energy convergence with best kgrid
    encut(args, best_kgrid=kpt_conv.best_kpoints)

    # Write best KPOINTS file
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

    if args.incar is None:
        extra_vars = args.extra_vars
    else:
        in_obj = incar.read_incar(args.incar)
        extra_vars = in_obj.variables

    st = load_poscar()
    relax_st = pychemia.code.vasp.task.IonRelaxation(
        structure=st,
        workdir=".",
        target_forces=args.target_forces,
        executable=args.vasp_exe,
        encut=encut_val,
        kp_grid=gridline,
        pspdir=args.pspdir,
        max_calls=args.max_calls,
        extra_vars=extra_vars,
        heterostructure=args.heterostructure,
        relax_cell=args.relax_cell,
        psp_options=args.psp_options,
        make_potcar=args.make_potcar,
        auto_ibrion=args.auto_ibrion,
        fire=args.fire
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
            "-np",
            default=1,
            type=int,
            help="Number of processors",
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
            "-psp_options",
            help="Pseudopotential options. ",
            default=None,
            type=json.loads,
        )
        parser_kgrid.add_argument(
            "-energy_tolerance",
            default=1e-3,
            type=float,
            help="The energy difference required for convergence per atom",
        )
        parser_kgrid.add_argument(
            "-heterostructure",
            help="Keep repeating order of atoms in POSCAR for POTCAR generation?",
            action="store_true",
        )
        parser_kgrid.add_argument(
            "-incar",
            default=None,
            type=str,
            help="INCAR file name. If not provided will be generated automatically.",
        )
        parser_kgrid.add_argument(
            "-vasp_exe", default="vasp_std", type=str, help="vasp executable"
        )
        parser_kgrid.add_argument(
            "-make_potcar",
            help="Flag to automatically generate POTCAR.",
            action="store_true",
        )

        parser_kgrid.set_defaults(func=kgrid)

        # parser for encut
        parser_encut = subparsers.add_parser("encut", help="Energy cut-off convergence")
        parser_encut.add_argument(
            "-np",
            default=1,
            type=int,
            help="Number of processors",
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
            "-psp_options",
            help="Pseudopotential options. ",
            default=None,
            type=json.loads,
        )
        parser_encut.add_argument(
            "-energy_tolerance",
            default=1e-3,
            type=float,
            help="The energy difference required for convergence per atom",
        )
        parser_encut.add_argument(
            "-heterostructure",
            help="Keep repeating order of atoms in POSCAR for POTCAR generation?",
            action="store_true",
        )
        parser_encut.add_argument(
            "-incar",
            default=None,
            type=str,
            help="INCAR file name. If not provided will be generated automatically.",
        )
        parser_encut.add_argument(
            "-vasp_exe", default="vasp_std", type=str, help="vasp executable"
        )
        parser_encut.add_argument(
            "-make_potcar",
            help="Flag to automatically generate POTCAR.",
            action="store_true",
        )

        parser_encut.set_defaults(func=encut)

        # parser for both k-grid and encut
        parser_complete = subparsers.add_parser("complete", help="Complete convergence")
        parser_complete.add_argument(
            "-np",
            default=1,
            type=int,
            help="Number of processors",
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
            "-psp_options",
            help="Pseudopotential options. ",
            default=None,
            type=json.loads,
        )
        parser_complete.add_argument(
            "-energy_tolerance",
            default=1e-3,
            type=float,
            help="The energy difference required for convergence per atom",
        )
        parser_complete.add_argument(
            "-heterostructure",
            help="Keep repeating order of atoms in POSCAR for POTCAR generation?",
            action="store_true",
        )
        parser_complete.add_argument(
            "-incar",
            default=None,
            type=str,
            help="INCAR file name. If not provided will be generated automatically.",
        )
        parser_complete.add_argument(
            "-vasp_exe",
            default="vasp_std",
            type=str,
            help="vasp executable",
        )
        parser_complete.add_argument(
            "-make_potcar",
            help="Flag to automatically generate POTCAR.",
            action="store_true",
        )

        parser_complete.set_defaults(func=complete)

        # parser for ionic relaxation
        parser_relax = subparsers.add_parser("relax", help="Ionic relaxation")
        parser_relax.add_argument(
            "-np",
            default=1,
            type=int,
            help="Number of processors",
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
            "-target_forces",
            default=1e-4,
            type=float,
            help="Target force difference",
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
        parser_relax.add_argument(
            "-heterostructure",
            help="Keep repeating order of atoms in POSCAR for POTCAR generation?",
            action="store_true",
        )
        parser_relax.add_argument(
            "-relax_cell",
            help="Optimize cell parameters as well.",
            action="store_true",
        )
        parser_relax.add_argument(
            "-incar",
            default=None,
            type=str,
            help="INCAR file name. If not provided will be generated automatically.",
        )
        parser_relax.add_argument(
            "-vasp_exe", default="vasp_std", type=str, help="vasp executable"
        )
        parser_relax.add_argument(
            "-make_potcar",
            help="Flag to automatically generate POTCAR.",
            action="store_true",
        )
        parser_relax.add_argument(
            "-auto_ibrion",
            help="Flag to set adaptive IBRION update.",
            action="store_true",
        )
        parser_relax.add_argument(
            "-fire",
            help="Flag to use FIRE algorithm for relaxation. Requires -auto_ibrion. ",
            action="store_true",
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
        print("Usage: ConvergenceTest.py [-h]")
        print("{kgrid, encut, complete, relax}")
