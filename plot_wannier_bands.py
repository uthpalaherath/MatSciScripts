#!/usr/bin/env python

""" Wannier90 band structure plotter.

This script plots the wannier90 band structure.

Author: Uthpala Herath
"""

import numpy as np
import matplotlib.pyplot as plt
import re
import argparse
import os
import sys
from argparse import RawTextHelpFormatter


def plot_wannier_bands(
    outcar="OUTCAR.scf", savefig="wannier90_bands.png", show=True, elimit=None
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

    # plotting wannier bands
    x = []
    y = []
    with open("wannier90_band.dat", "r") as f:
        lines = f.readlines()
        x.append([])
        y.append([])
        i = 0
        for line in lines[0:-1]:
            if line != "\n":
                x[i].append(float(line.split()[0]))
                y[i].append(float(line.split()[1]))
            else:
                x.append([])
                y.append([])
                i += 1

    x = np.array(x)
    y = np.array(y)

    # Plotting
    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111)

    for i in range(len(x)):
        ax.plot(x[i], y[i] - EFERMI, color="blue")

    ax.set_xticks(ticks)
    ax.set_xticklabels(knames)
    ax.set_xlabel(r"$k$-path")
    ax.set_ylabel(r"$E-E_F$ (eV)")
    ax.axhline(y=0, color="black", ls="--")
    for xc in ticks:
        ax.axvline(x=xc, color="k")

    ax.set_xlim(x.min(), x.max())
    if elimit:
        ax.set_ylim(elimit)

    fig.tight_layout()

    if show:
        plt.show()
        return None, None
    else:
        plt.savefig(savefig, bbox_inches="tight")
        return fig, ax


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-outcar", type=str, help="SCF OUTCAR file.", default="OUTCAR.scf"
    )
    parser.add_argument(
        "-elimit", type=float, nargs=2, help="Energy axis range", default=None
    )
    args = parser.parse_args()
    plot_wannier_bands(outcar=args.outcar, elimit=args.elimit)
