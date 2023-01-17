#!/usr/bin/env python

"""Matplotlib plotter template.

A template for creating plots with matplotlib.

- Uthpala Herath
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import matplotlib.pyplot as plt
from cycler import cycler
from my_plot import set_size

# import pandas as pd
import numpy as np
import sys

# ------------------------- INITIALIZATION -------------------------

# Stylesheet
plt.style.use("~/dotfiles/matplotlib/prb.mplstyle")

# Line and marker cyclers
line_cycler = cycler(
    color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]
) + cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."])

marker_cycler = (
    cycler(
        color=[
            "#E69F00",
            "#56B4E9",
            "#009E73",
            "#0072B2",
            "#D55E00",
            "#CC79A7",
            "#F0E442",
        ]
    )
    + cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"])
    + cycler(marker=["4", "2", "3", "1", "+", "x", "."])
)


# Line styles
linestyle = ["-", "--", "-.", ":", "-", "--", "-."]

# Span in two columns with (2,3) grid
# width = (510 * 2 - 18) / 3
width = 768
fraction = 1

# arguments for y-limit
if len(sys.argv) > 1:
    ymin, ymax = float(sys.argv[1]), float(sys.argv[2])
else:
    ymin = -15
    ymax = 30

    # ------------------------- DATA -------------------------

k = np.arange(0, 21)
bands = np.zeros((21, 20), dtype="float64")

# Text file import
with open("./GW_band1001.out", "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    raw_line = [float(x) for x in line.split()[4:][1::2]]
    bands[i, :] = raw_line


# ------------------------- PLOTTING -------------------------

fig, ax = plt.subplots(1, 1, figsize=set_size(width, fraction))
# ax.set_prop_cycle(line_cycler)
for i in range(len(bands.T)):
    ax.plot(k, bands[:, i], marker="o", markersize=3, linestyle=linestyle[0])

# ax.set_prop_cycle(line_cycler)
ax.set_title(r"GW bandstructure")
ax.set_xlabel(r"$k-point$")
ax.set_ylabel(r"$E$-$E_F$ (eV)")
ax.set_xlim([0, k[-1]])
ax.set_ylim([ymin, ymax])
ax.axhline(y=0, color="gray", linestyle="--")
ax.set_xticks(k)
ax.grid()
# lgd = ax.legend()

# legend sizes
# lgd.legendHandles[0]._sizes = [30]
# lgd.legendHandles[1]._sizes = [30]

fig.savefig("bands.pdf")
