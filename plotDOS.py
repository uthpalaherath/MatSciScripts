#!/usr/bin/env python

"""DOS Plotter

This script plots Partial DOS from VASP outputs.

Author: Uthpala Herath

"""


import argparse
import collections
import re
import sys
from argparse import RawTextHelpFormatter

from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp import Vasprun


def plot_total(args):
    """
    This method plots the total DOS.

    Parameters
    ----------

    show : boolean
        Flag to show the plot on screen

    xlim : float
        Range of x axis

    ylim : float
        Range of y axis


    """
    print("Plotting Total DOS...")
    vasprun = Vasprun("vasprun.xml", parse_dos=True)
    tdos = vasprun.tdos
    plot = DosPlotter()
    plot.add_dos("Total DOS", tdos)
    if args.show:
        plot.show(xlim=args.xlim, ylim=args.ylim)
    plot.save_plot("dos_total.png", img_format="png", xlim=args.xlim, ylim=args.ylim)


def plot_dos_atoms(args):
    """
    This method plots the partial and total DOS
    for atoms.

    Parameters
    ----------

    show : boolean
        Flag to show the plot on screen

    xlim : float
        Range of x axis

    ylim : float
        Range of y axis

    total : boolean
        Plot the total DOS


    """
    print("Plotting atom projected DOS...")
    vasprun = Vasprun("vasprun.xml", parse_dos=True)
    cdos = vasprun.complete_dos
    element_dos = cdos.get_element_dos()
    plot = DosPlotter()
    plot.add_dos_dict(element_dos)
    if args.total:
        plot.add_dos("Total DOS", cdos)
        if args.show:
            plot.show(xlim=args.xlim, ylim=args.ylim)
        plot.save_plot(
            "dos_atomic_total.png", img_format="png", xlim=args.xlim, ylim=args.ylim
        )

    else:
        if args.show:
            plot.show(xlim=args.xlim, ylim=args.ylim)
        plot.save_plot(
            "dos_atomic.png", img_format="png", xlim=args.xlim, ylim=args.ylim
        )


def plot_dos_orbitals(args):
    """
    This method plots the DOS for orbitals of
    all the atoms or a specific atom.

    Parameters
    ----------

    show : boolean
        Flag to show the plot on screen

    xlim : float
        Range of x axis

    ylim : float
        Range of y axis

    total : boolean
        Plot the orbital projected DOS from all the atoms

    atom : string
        Plot the orbital projected DOS from a specified atom

    """
    print("Plotting orbital projected DOS...")
    vasprun = Vasprun("vasprun.xml", parse_dos=True)
    cdos = vasprun.complete_dos

    if args.atom:
        spd_dos_atom = cdos.get_element_spd_dos(args.atom)
        plot = DosPlotter()
        plot.add_dos_dict(spd_dos_atom)
        save_str = "dos_orbital_" + args.atom + ".png"
        if args.show:
            plot.show(xlim=args.xlim, ylim=args.ylim)
        plot.save_plot(save_str, img_format="png", xlim=args.xlim, ylim=args.ylim)

    elif args.total:
        spd_dos = cdos.get_spd_dos()
        plot = DosPlotter()
        plot.add_dos_dict(spd_dos)
        if args.show:
            plot.show(xlim=args.xlim, ylim=args.ylim)
        plot.save_plot(
            "dos_orbital_total.png", img_format="png", xlim=args.xlim, ylim=args.ylim
        )


def plot_dos_d(args):
    """
    This method plots the DOS of d-t2g and d-eg orbitals
    of an atom.
    Remember that the indexing starts with zero.

    Parameters
    ----------

    minmax : int
        The starting and ending index of the atoms to plot d orbital
        projected DOS

    list : int
        List of atoms to plot d orbital projected DOS

    atom : str
        Name of atom to plot d orbital projected DOS

    show : boolean
        Flag to show the plot on screen

    xlim : float
        Range of x axis

    ylim : float
        Range of y axis

    """

    print("Summing:")
    vasprun = Vasprun("vasprun.xml", parse_dos=True)
    plot = DosPlotter()
    data = {}

    # creating dictionary for data values
    if args.minmax:
        for count, i in enumerate(range(args.minmax[0], args.minmax[1])):
            data[
                "data{0}".format(count)
            ] = vasprun.complete_dos.get_site_t2g_eg_resolved_dos(
                vasprun.structures[0][i]
            )
            print(vasprun.structures[0][i])
    elif args.list:
        for count, i in enumerate(args.list):
            data[
                "data{0}".format(count)
            ] = vasprun.complete_dos.get_site_t2g_eg_resolved_dos(
                vasprun.structures[0][i]
            )
            print(vasprun.structures[0][i])

    elif args.atom:
        cdos = vasprun.complete_dos
        atomlist = cdos.structure.formula.split()
        matchlist = []
        for i in atomlist:
            match = re.match(r"([a-z]+)([0-9]+)", i, re.I)
            if match:
                matchlist.append(match.groups())

        matchlistdic = collections.OrderedDict()
        for count in range(len(matchlist)):
            matchlistdic[matchlist[count][0]] = matchlist[count][1]
        keys = list(matchlistdic.keys())
        values = list(matchlistdic.values())
        values = [int(i) for i in values]

        atomindex = keys.index(args.atom)
        minindex = sum(values[:atomindex])
        maxindex = sum(values[: atomindex + 1])

        for count, i in enumerate(range(minindex, maxindex)):
            data[
                "data{0}".format(count)
            ] = vasprun.complete_dos.get_site_t2g_eg_resolved_dos(
                vasprun.structures[0][i]
            )
            print(vasprun.structures[0][i])

    # initializing sum
    data_sum_t2g = data["data0"]["t2g"]
    data_sum_eg = data["data0"]["e_g"]

    for k in data.keys():
        if k != "data0":
            data_sum_t2g = data_sum_t2g + data[k]["t2g"]
            data_sum_eg = data_sum_eg + data[k]["e_g"]

    print("\nPlotting d orbital projected DOS...")
    if args.atom:
        str_t2g = args.atom + " " + "d-t$_{2g}$"
        str_eg = args.atom + " " + "d-e$_g$"
        str_save = "dos_d_" + args.atom + ".png"
        plot.add_dos(str_t2g, data_sum_t2g)
        plot.add_dos(str_eg, data_sum_eg)
        if args.show:
            plot.show(xlim=args.xlim, ylim=args.ylim)
        plot.save_plot(str_save, img_format="png", xlim=args.xlim, ylim=args.ylim)

    else:
        plot.add_dos("d-t$_{2g}$", data_sum_t2g)
        plot.add_dos("d-e$_g$", data_sum_eg)
        if args.show:
            plot.show(xlim=args.xlim, ylim=args.ylim)
        plot.save_plot("dos_d.png", img_format="png", xlim=args.xlim, ylim=args.ylim)


if __name__ == "__main__":
    args = sys.argv[1:]
    if args:
        parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=RawTextHelpFormatter,
        )
        subparsers = parser.add_subparsers(help="sub-command help")

        # parser for total dos plotter
        parser_total = subparsers.add_parser("total", help="Total DOS  ")
        parser_total.add_argument(
            "-xlim", type=float, default=[-5, 5], nargs="+", help="Range of x-axis"
        )
        parser_total.add_argument(
            "-ylim", type=float, default=None, nargs="+", help="Range of y-axis"
        )
        parser_total.add_argument(
            "-show", help="Show plot on screen", action="store_true"
        )
        parser_total.set_defaults(func=plot_total)

        # parser for atomic dos plotter
        parser_dos_atom = subparsers.add_parser("atomic", help="Atom projected DOS")
        parser_dos_atom.add_argument(
            "-xlim", type=float, default=[-5, 5], nargs="+", help="Range of x-axis"
        )
        parser_dos_atom.add_argument(
            "-ylim", type=float, default=None, nargs="+", help="Range of y-axis"
        )
        parser_dos_atom.add_argument(
            "-show", help="Show plot on screen", action="store_true"
        )
        parser_dos_atom.add_argument(
            "-total", help="Plot the total DOS", action="store_true"
        )
        parser_dos_atom.set_defaults(func=plot_dos_atoms)

        # parser for orbital dos plotter
        parser_dos_orbital = subparsers.add_parser(
            "orbital", help="Orbital projected DOS"
        )
        parser_dos_orbital.add_argument(
            "-xlim", type=float, default=[-5, 5], nargs="+", help="Range of x-axis"
        )
        parser_dos_orbital.add_argument(
            "-ylim", type=float, default=None, nargs="+", help="Range of y-axis"
        )
        parser_dos_orbital.add_argument(
            "-show", help="Show plot on screen", action="store_true"
        )
        group1 = parser_dos_orbital.add_mutually_exclusive_group(required=True)
        group1.add_argument(
            "-atom",
            type=str,
            default=None,
            help="Select atom to plot orbital projected dos",
        )
        group1.add_argument("-total", help="Plot the total DOS", action="store_true")
        parser_dos_orbital.set_defaults(func=plot_dos_orbitals)

        # parser for d orbital dos plotter
        parser_d = subparsers.add_parser("orbital_d", help="d orbital projected DOS")
        group2 = parser_d.add_mutually_exclusive_group(required=True)
        group2.add_argument(
            "-minmax",
            type=int,
            default=None,
            nargs="+",
            help="Starting and ending indexes of atoms: [start end]",
        )
        group2.add_argument(
            "-l",
            "--list",
            type=int,
            default=None,
            nargs="+",
            help="List of indexes of atoms",
        )
        group2.add_argument(
            "-atom",
            type=str,
            default=None,
            help="Select atom to plot d-orbital projected dos",
        )
        parser_d.add_argument(
            "-xlim", type=float, default=[-5, 5], nargs="+", help="Range of x-axis"
        )
        parser_d.add_argument(
            "-ylim", type=float, default=None, nargs="+", help="Range of y-axis"
        )
        parser_d.add_argument("-show", help="Show plot on screen", action="store_true")
        parser_d.set_defaults(func=plot_dos_d)

        args = parser.parse_args()
        args.func(args)

    else:
        print("Usage: plotDOS.py [-h]")
        print("{total, atomic, orbital, orbital_d}")
