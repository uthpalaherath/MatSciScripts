#!/usr/bin/env python

"""Matplotlib plotter template.

A template for creating plots with matplotlib.

- Uthpala Herath
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import matplotlib.pyplot as plt
from cycler import cycler
from my_plot import set_size
import numpy as np
import glob
import sys

# ------------------------- INITIALIZATION -------------------------

# arguments for y-limit
if len(sys.argv) > 1:
    ymin, ymax = float(sys.argv[1]), float(sys.argv[2])
else:
    ymin = -15
    ymax = 30

# Stylesheet
plt.style.use("~/dotfiles/matplotlib/prb.mplstyle")

# Span in two columns with (2,3) grid
# width = (510 * 2 - 18) / 3
width = 768
fraction = 1

# Line and marker cyclers
NUM_COLORS = 20
cm = plt.get_cmap("tab20")

# global axis
fig, ax = plt.subplots(1, 1, figsize=set_size(width, fraction))
ax.set_prop_cycle(color=[cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])

# Line styles
linestyle = ["-", "--", "-.", ":", "-", "--", "-."]

# ------------------------- DATA + PLOTTING -------------------------

# user input
kpoints_per_segment = 21
segments = 10
nbands = 20
knames = ["$\Gamma$", "X", "W", "K", "$\Gamma$", "L", "U", "W", "L", "K", "X"]

k = np.arange(0, segments * kpoints_per_segment)
k_per_segment = np.arange(0, kpoints_per_segment)

bands = np.zeros((kpoints_per_segment, nbands), dtype="float64")
filelist = sorted(glob.glob("GW_band*"))
hsp_coords = []

for fi in filelist:
    # Text file import
    with open(fi, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        raw_line = [float(x) for x in line.split()[4:][1::2]]
        bands[i, :] = raw_line

    for i in range(len(bands.T)):
        ax.plot(
            k_per_segment, bands[:, i], marker="o", markersize=1, linestyle=linestyle[0]
        )
    k_per_segment_init = k_per_segment[0]
    hsp_coords.append(k_per_segment[0])
    k_per_segment += kpoints_per_segment - 1

hsp_coords.append(k_per_segment_init)

ax.set_title(r"GW bandstructure")
ax.set_xlabel(r"$k-path$")
ax.set_ylabel(r"$E$-$E_F$ (eV)")
ax.set_xlim([0, hsp_coords[-1]])
ax.set_ylim([ymin, ymax])
ax.axhline(y=0, color="gray", linestyle="--")
ax.set_xticks(hsp_coords)
ax.set_xticklabels(knames)
ax.grid()
fig.savefig("gw_bands.pdf")
