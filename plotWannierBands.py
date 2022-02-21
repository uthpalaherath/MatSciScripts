#!/usr/bin/env python
""" Wannier90 band structure plotter.

This script plots the wannier90 band structure and compares
it to the DFT band structure.

Authors: Uthpala Herath and Andres Tellez
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


def kpoint_conversion(value, kpoints1, kpoints2):
    """Finds transformation between DFT k-points and
    wannier90 k-points. Original implementation by Andres Tellez."""

    idx = 0
    for point in kpoints1:
        if value <= point:
            idx = kpoints1.index(point)
            break

    result = value - kpoints1[idx - 1]
    try:
        result *= (kpoints2[idx] - kpoints2[idx - 1]) / (
            kpoints1[idx] - kpoints1[idx - 1]
        )
    except ZeroDivisionError:
        result *= 1.0
    result += kpoints2[idx - 1]

    return result


def plot_bands(
    outcar="OUTCAR",
    savefig="wannier90_bands.pdf",
    show=False,
    elimit=None,
    compare=False,
):
    """This method plots wannier bands"""

    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams["font.family"] = "Arial"
    plt.rc("font", size=22)  # controls default text sizes
    plt.rc("axes", titlesize=22)  # fontsize of the axes title
    plt.rc("axes", labelsize=22)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=22)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=22)  # fontsize of the tick labels

    # Reading Fermi from OUTCAR
    fi = open(outcar, "r")
    for line in fi:
        if re.search("Fermi energy", line) or re.search("E-fermi", line):
            line_fermi = line
    val = re.search(r"(\-?\d+\.?\d*)", line_fermi)
    EFERMI = float(val.group(1))
    print("Fermi energy = {:4.4f} eV".format(EFERMI))

    # Reading labels
    fi = open("wannier90_band.gnu", "r")
    data = fi.read()
    fi.close()

    label_line = re.findall(r'xtics\s*\(([\s"A-Z0-9.,|]*)\)', data)
    label_line_split = label_line[0].split(",")

    ticks = []
    knames = []

    for i in label_line_split:
        knames.append(i.split()[0])
        ticks.append(float(i.split()[1]))

    knames = [i.strip('"') for i in knames]

    # Getting wannier band data
    x = []
    y = []
    with open("wannier90_band.dat", "r") as f:
        lines = f.readlines()
        x.append([])
        y.append([])
        i = 0
        for line in lines:
            if line != "  \n":
                x[i].append(float(line.split()[0]))
                y[i].append(float(line.split()[1]))
            else:
                x.append([])
                y.append([])
                i += 1

    # There sometimes is an extra line in wannier90_band.dat.
    # We will remove it if it is there.
    if not x[-1]:
        x = np.array(x[:-1])
        y = np.array(y[:-1])
    else:
        x = np.array(x)
        y = np.array(y)

    if not compare:
        # Plotting
        fig = plt.figure(figsize=(13, 9))
        ax = fig.add_subplot(111)
        fig.tight_layout()

        for i in range(len(x)):
            ax.plot(x[i], y[i] - EFERMI, color="blue")

        ax.set_xlim(x.min(), x.max())
        if elimit:
            ax.set_ylim(elimit)

        ax.set_xticks(ticks)
        ax.set_xticklabels(knames)
        ax.set_xlabel(r"$k$-path")
        ax.set_ylabel(r"$E-E_F$ (eV)")
        ax.axhline(y=0, color="black", ls="--")
        ax.grid()
        for xc in ticks:
            ax.axvline(x=xc, color="k")

        if show:
            plt.show()
            return None, None
        else:
            plt.savefig(savefig, bbox_inches="tight")
            return fig, ax
    else:
        # Comparison with DFT bands from PyProcar
        # Input axis from PyProcar plot to get x axis ticks

        fig, ax = pyprocar.bandsplot(
            "PROCAR",
            outcar=outcar,
            kpointsfile="KPOINTS",
            elimit=elimit,
            mode="plain",
            color="red",
            show=False,
            verbose=False,
        )

        ticks_pyprocar = [ax.lines[-i].get_xdata()[0] for i in range(2, len(ticks) + 2)]
        ticks_pyprocar.sort()

        x_new = np.zeros((x.shape), dtype="float64")

        for ix in range(len(x)):
            for iix in range(x.shape[1]):
                x_new[ix, iix] = kpoint_conversion(x[ix, iix], ticks, ticks_pyprocar)

        # Plotting
        for i in range(len(x_new)):
            ax.plot(x_new[i], y[i] - EFERMI, color="blue", linewidth=2)

        ax.set_xlim(x_new.min(), x_new.max())
        if elimit:
            ax.set_ylim(elimit)

        ax.set_xticks(ticks)
        ax.set_xticklabels(knames)
        ax.set_xlabel(r"$k$-path")
        ax.set_ylabel(r"$E-E_F$ (eV)")
        ax.axhline(y=0, color="black", ls="--")
        ax.grid()
        # for xc in ticks:
        #     ax.axvline(x=xc, color="k")

        # legend
        custom_lines = [
            Line2D([0], [0], color="red", lw=2),
            Line2D([0], [0], color="blue", lw=2),
        ]
        plt.legend(custom_lines, ["DFT", "wannier90"])

        # replot with PyProcar with new axis
        if show:
            savefig = None
        pyprocar.bandsplot(
            "PROCAR",
            outcar=outcar,
            kpointsfile="KPOINTS",
            elimit=elimit,
            mode="plain",
            color="red",
            show=show,
            ax=ax,
            savefig=savefig,
            verbose=False,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument("-outcar", type=str, help="SCF OUTCAR file.", default="OUTCAR")
    parser.add_argument(
        "-elimit", type=float, nargs=2, help="Energy axis range", default=None
    )
    parser.add_argument("-show", action="store_true", help="Flag to show plot.")
    parser.add_argument(
        "-compare",
        action="store_true",
        help="Flag to compare wannier90 bands with DFT bands (Requires PyProcar with PROCAR, KPOINTS and SCF OUTCAR).",
    )
    args = parser.parse_args()
    plot_bands(
        outcar=args.outcar, elimit=args.elimit, show=args.show, compare=args.compare
    )
