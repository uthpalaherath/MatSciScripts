#!/usr/bin/env python
""" Multiple DOS plotter for ABO3 compounds.

This script plots multiple orbitals of B and O of
ABO3 type perovskites using PyProcar.

Authors: Uthpala Herath
"""

import numpy as np
import matplotlib.pyplot as plt
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
    elimit=None,
    ylimit=None,
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

    # spin polarized
    sp = False

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

    # Create new figure environment
    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111)
    fig.tight_layout()

    # Degenerate case
    if degenerate:
        # Here the eg and t2g orbital symmetry is conserved
        # t2g orbital data
        _, ax1 = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[4, 5, 7],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False,
        )

        # check if spin polarized or not
        if len(ax1.lines[0].get_ydata()) == len(ax1.lines[1].get_ydata()):
            sp = True

        if not sp:
            line = ax1.lines[0]
            x1 = line.get_xdata()
            y1 = line.get_ydata()
        else:
            lineup = ax1.lines[0]
            linedn = ax1.lines[1]
            x1 = lineup.get_xdata()
            y1_up = lineup.get_ydata()
            y1_dn = linedn.get_ydata()

        # eg orbital data
        _, ax2 = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[6, 8],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False,
        )
        if not sp:
            line = ax2.lines[0]
            x2 = line.get_xdata()
            y2 = line.get_ydata()
        else:
            lineup = ax2.lines[0]
            linedn = ax2.lines[1]
            x2 = lineup.get_xdata()
            y2_up = lineup.get_ydata()
            y2_dn = linedn.get_ydata()

        # O-p orbital data
        _, ax3 = pyprocar.dosplot(
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
            plt_show=False,
        )
        if not sp:
            line = ax3.lines[0]
            x3 = line.get_xdata()
            y3 = line.get_ydata()
        else:
            lineup = ax3.lines[0]
            linedn = ax3.lines[1]
            x3 = lineup.get_xdata()
            y3_up = lineup.get_ydata()
            y3_dn = linedn.get_ydata()

        # Plotting
        if not sp:
            ax.plot(x1, y1, label=species[1] + "-d$_{t2g}$", color="blue")
            ax.plot(x2, y2, label=species[1] + "-d$_{eg}$", color="red")
            ax.plot(x3, y3, label="O-p", color="green")

            if ylimit:
                ax.set_ylim(ylimit)
            else:
                ax.set_ylim(
                    min(min(y1), min(y2), min(y3)), max(max(y1), max(y2), max(y3))
                )

        else:
            # t2g
            ax.plot(x1, y1_up, label=species[1] + r"-d$_{t2g} \uparrow$", color="blue")
            ax.plot(
                x1,
                y1_dn,
                label=species[1] + r"-d$_{t2g} \downarrow$",
                color="blue",
                linestyle="dotted",
            )

            # eg
            ax.plot(x2, y2_up, label=species[1] + r"-d$_{eg} \uparrow$", color="red")
            ax.plot(
                x2,
                y2_dn,
                label=species[1] + r"-d$_{eg} \downarrow$",
                color="red",
                linestyle="dotted",
            )

            # O-p
            ax.plot(x3, y3_up, label=r"O-p $\uparrow$", color="green")
            ax.plot(
                x3,
                y3_dn,
                label=r"O-p $\downarrow$",
                color="green",
                linestyle="dotted",
            )

            if ylimit:
                ax.set_ylim(ylimit)
            else:
                ax.set_ylim(
                    min(min(y1_dn), min(y2_dn), min(y3_dn)),
                    max(max(y1_up), max(y2_up), max(y3_up)),
                )

    # Non-degenerate case
    else:
        # Here the eg and t2g orbital symmetry is broken.
        # Each "B" orbital is plot separately

        # d-xy  orbital data
        _, ax1 = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[4],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False,
        )

        # check if spin polarized or not
        if len(ax1.lines[0].get_ydata()) == len(ax1.lines[1].get_ydata()):
            sp = True

        if not sp:
            line = ax1.lines[0]
            x1 = line.get_xdata()
            y1 = line.get_ydata()
        else:
            lineup = ax1.lines[0]
            linedn = ax1.lines[1]
            x1 = lineup.get_xdata()
            y1_up = lineup.get_ydata()
            y1_dn = linedn.get_ydata()

        # d-yz orbital data
        _, ax2 = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[5],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False,
        )
        if not sp:
            line = ax2.lines[0]
            x2 = line.get_xdata()
            y2 = line.get_ydata()
        else:
            lineup = ax2.lines[0]
            linedn = ax2.lines[1]
            x2 = lineup.get_xdata()
            y2_up = lineup.get_ydata()
            y2_dn = linedn.get_ydata()

        # d-xz orbital data
        _, ax3 = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[7],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False,
        )
        if not sp:
            line = ax3.lines[0]
            x3 = line.get_xdata()
            y3 = line.get_ydata()
        else:
            lineup = ax3.lines[0]
            linedn = ax3.lines[1]
            x3 = lineup.get_xdata()
            y3_up = lineup.get_ydata()
            y3_dn = linedn.get_ydata()

        # d-y2 orbital data
        _, ax4 = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[6],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False,
        )
        if not sp:
            line = ax4.lines[0]
            x4 = line.get_xdata()
            y4 = line.get_ydata()
        else:
            lineup = ax4.lines[0]
            linedn = ax4.lines[1]
            x4 = lineup.get_xdata()
            y4_up = lineup.get_ydata()
            y4_dn = linedn.get_ydata()

        # d-x2-y2 orbital data
        _, ax5 = pyprocar.dosplot(
            filename="vasprun.xml",
            mode="parametric_line",
            orbitals=[8],
            atoms=atoms_d,
            plot_total=False,
            plt_show=False,
        )
        if not sp:
            line = ax5.lines[0]
            x5 = line.get_xdata()
            y5 = line.get_ydata()
        else:
            lineup = ax5.lines[0]
            linedn = ax5.lines[1]
            x5 = lineup.get_xdata()
            y5_up = lineup.get_ydata()
            y5_dn = linedn.get_ydata()

        # O-p orbital data
        _, ax6 = pyprocar.dosplot(
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
            plt_show=False,
        )
        if not sp:
            line = ax6.lines[0]
            x6 = line.get_xdata()
            y6 = line.get_ydata()
        else:
            lineup = ax6.lines[0]
            linedn = ax6.lines[1]
            x6 = lineup.get_xdata()
            y6_up = lineup.get_ydata()
            y6_dn = linedn.get_ydata()

        # Plotting
        if not sp:
            ax.plot(x1, y1, label=species[1] + "-d$_{xy}$", color="blue")
            ax.plot(x2, y2, label=species[1] + "-d$_{yz}$", color="cyan")
            ax.plot(x3, y3, label=species[1] + "-d$_{xz}$", color="lightblue")
            ax.plot(x4, y4, label=species[1] + "-d$_{z^2}$", color="red")
            ax.plot(x5, y5, label=species[1] + "-d$_{x^2-t^2}$", color="maroon")
            ax.plot(x6, y6, label="O-p", color="green")

            if ylimit:
                ax.set_ylim(ylimit)
            else:
                ax.set_ylim(
                    min(min(y1), min(y2), min(y3), min(y4), min(y5), min(y6)),
                    max(max(y1), max(y2), max(y3), max(y4), max(y5), max(y6)),
                )

        else:
            # d-xy
            ax.plot(x1, y1_up, label=species[1] + r"-d$_{xy} \uparrow$", color="blue")
            ax.plot(
                x1,
                y1_dn,
                label=species[1] + r"-d$_{xy} \downarrow$",
                color="blue",
                linestyle="dotted",
            )

            # d-yz
            ax.plot(x2, y2_up, label=species[1] + r"-d$_{yz} \uparrow$", color="cyan")
            ax.plot(
                x2,
                y2_dn,
                label=species[1] + r"-d$_{yz} \downarrow$",
                color="cyan",
                linestyle="dotted",
            )

            # d-xz
            ax.plot(
                x3, y3_up, label=species[1] + r"-d$_{xz} \uparrow$", color="lightblue"
            )
            ax.plot(
                x3,
                y3_dn,
                label=species[1] + r"-d$_{xz} \downarrow$",
                color="lightblue",
                linestyle="dotted",
            )

            # d-z2
            ax.plot(x4, y4_up, label=species[1] + r"-d$_{z^2} \uparrow$", color="red")
            ax.plot(
                x4,
                y4_dn,
                label=species[1] + r"-d$_{z^2} \downarrow$",
                color="red",
                linestyle="dotted",
            )

            # d-x2-y2
            ax.plot(
                x5, y5_up, label=species[1] + r"-d$_{x^2-y^2} \uparrow$", color="maroon"
            )
            ax.plot(
                x5,
                y5_dn,
                label=species[1] + r"-d$_{x^2-y^2} \downarrow$",
                color="maroon",
                linestyle="dotted",
            )

            # O-p
            ax.plot(x6, y6_up, label=r"O-p $\uparrow$", color="green")
            ax.plot(
                x6,
                y6_dn,
                label=r"O-p $\downarrow$",
                color="green",
                linestyle="dotted",
            )

            if ylimit:
                ax.set_ylim(ylimit)
            else:
                ax.set_ylim(
                    min(
                        min(y1_dn),
                        min(y2_dn),
                        min(y3_dn),
                        min(y4_dn),
                        min(y5_dn),
                        min(y6_dn),
                    ),
                    max(
                        max(y1_up),
                        max(y2_up),
                        max(y3_up),
                        max(y4_up),
                        max(y5_up),
                        max(y6_up),
                    ),
                )

    # Common plot elements
    if elimit:
        ax.set_xlim(elimit)

    ax.set_xlabel(r"E-E$_F$ (eV)")
    ax.set_ylabel(r"DOS")
    ax.set_title(r"DFT DOS")
    ax.axhline(y=0, color="black", ls="-")
    ax.axvline(x=0, color="black", ls="--")
    ax.grid(color="gainsboro", ls="--", lw=0.6)
    ax.legend(loc="best")
    plt.savefig("DFT_DOS.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-elim", "--elimit", type=float, nargs=2, help="Energy axis range", default=None
    )
    parser.add_argument(
        "-ylim", "--ylimit", type=float, nargs=2, help="DOS axis range", default=None
    )
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
        ylimit=args.ylimit,
        atoms=args.atoms,
        degenerate=args.degenerate,
    )
