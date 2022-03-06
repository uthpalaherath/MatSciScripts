#!/usr/bin/env python
""" Multiple DOS plotter for ABO3 compounds.

This script plots multiple orbitals of B and O of
ABO3 type perovskites using PyProcar.

Authors: Uthpala Herath
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re
import argparse
from argparse import RawTextHelpFormatter
import pyprocar
import os
import sys

import warnings

warnings.filterwarnings("ignore")

# Show no PyProcar output
# f = open(os.devnull, "w")
# sys.stdout = f


def plot_dos(
    savefig="DFT_DOS.pdf",
    show=False,
    elimit=None,
    ylim=None,
    atoms=None,
    degenerate=False,
):
    """This method plots dos"""

    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams["font.family"] = "Arial"
    plt.rc("font", size=22)  # controls default text sizes
    plt.rc("axes", titlesize=22)  # fontsize of the axes title
    plt.rc("axes", labelsize=22)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=22)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=22)  # fontsize of the tick labels

    # Reading POSCAR
    file = open("POSCAR", "r")
    data = file.readlines()
    file.close()
    species = data[5].split()
    species_count = data[6].split()
    species_count = [int(i) for i in species_count]

    if atoms:
        atoms_d = [i + species_count[0] - 1 for i in atoms]
    else:
        atoms_d = [*np.arange(species_count[0], species_count[0] + species_count[1])]

    # Degenerate case
    if degenerate:
        # Here the eg and t2g orbital symmetry is conserved
        # t2g orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[4, 5, 7],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x1 = line.get_xdata()
        y1 = line.get_ydata()

        # eg orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[6, 8],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x2 = line.get_xdata()
        y2 = line.get_ydata()

        # O-p orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[1, 2, 3],
            atoms=[
                *np.arange(
                    species_count[0] + species_count[1],
                    species_count[0] + species_count[1] + species_count[2],
                )
            ],
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x3 = line.get_xdata()
        y3 = line.get_ydata()

        # Plotting
        ax.plot(x1, y1, label=species[1] + "-d$_{t2g}$")
        ax.plot(x2, y2, label=species[1] + "-d$_{eg}$")
        ax.plot(x3, y3, label="O-p")

        if ylim:
            ax.set_ylim(min(y1, y2, y3), max(y1, y2, y3))
        if elimit:
            ax.set_xlim(elimit)

        ax.set_xlabel(r"E-E$_F$ (eV)")
        ax.set_ylabel(r"DOS")
        ax.set_title(r"DFT DOS")
        ax.axvline(x=0, color="black", ls="--")
        ax.grid(color="gainsboro", ls="--", lw=0.6)
        ax.legend(loc="best")

        if show:
            plt.show()

        plt.savefig("DFT_DOS.pdf")

    # Non-degenerate case
    else:
        # Here the eg and t2g orbital symmetry is broken.
        # Each "B" orbital is plot separately

        # d-xy  orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[4],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x1 = line.get_xdata()
        y1 = line.get_ydata()

        # d-yz orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[5],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x2 = line.get_xdata()
        y2 = line.get_ydata()

        # d-xz orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[7],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x3 = line.get_xdata()
        y3 = line.get_ydata()

        # d-y2 orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[6],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x4 = line.get_xdata()
        y4 = line.get_ydata()

        # d-x2-y2 orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[8],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x5 = line.get_xdata()
        y5 = line.get_ydata()

        # O-p orbital data
        fig, ax = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[1, 2, 3],
            atoms=[
                *np.arange(
                    species_count[0] + species_count[1],
                    species_count[0] + species_count[1] + species_count[2],
                )
            ],
            plot_total=False,
            plt_show=False
            # elimit=elimit,
            # labels=["O-up", "O-down"],
        )
        line = ax.lines[0]
        x6 = line.get_xdata()
        y6 = line.get_ydata()

        # Plotting
        ax.plot(x1, y1, label=species[1] + "-d$_{xy}$")
        ax.plot(x2, y2, label=species[1] + "-d$_{yx}$")
        ax.plot(x3, y3, label=species[1] + "-d$_{xz}$")
        ax.plot(x4, y4, label=species[1] + "-d$_{z^2}$")
        ax.plot(x5, y5, label=species[1] + "-d$_{x^2-y^2}$")
        ax.plot(x6, y4, label="O-p")

        if ylim:
            ax.set_ylim(min(y1, y2, y3, y4, y5, y6), max(y1, y2, y3, y4, y5, y6))
        if elimit:
            ax.set_xlim(elimit)

        ax.set_xlabel(r"E-E$_F$ (eV)")
        ax.set_ylabel(r"DOS")
        ax.axvline(x=0, color="black", ls="--")
        ax.grid(color="gainsboro", ls="--", lw=0.6)
        ax.legend(loc="best")

        if show:
            plt.show()
        plt.savefig("DFT_DOS.pdf")

    # custom legend
    # custom_lines = [
    #     Line2D([0], [0], color="red", lw=2),
    #     Line2D([0], [0], color="blue", lw=2),
    # ]
    # plt.legend(custom_lines, ["DFT", "wannier90"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-elimit", type=float, nargs=2, help="Energy axis range", default=None
    )
    parser.add_argument(
        "-ylim", type=float, nargs=2, help="DOS axis range", default=None
    )
    parser.add_argument("-show", action="store_true", help="Flag to show plot.")
    parser.add_argument(
        "-d",
        "--degenerate",
        action="store_true",
        help="Plot degenerate eg and t2g orbitals. Otherwise, each d orbital is plot separately.",
    )
    parser.add_argument(
        "-a",
        "--atoms",
        type=int,
        default=None,
        nargs="+",
        help="List of indexes of B atoms in ABO3 perovskite. Counting starts from 1 with the B atom in POSCAR. If not provided will use all.",
    )

    args = parser.parse_args()
    plot_dos(
        elimit=args.elimit,
        ylim=args.ylim,
        atoms=args.atoms,
        degenerate=args.degenerate,
        show=args.show,
    )
