#!/usr/bin/env python

""" Matplotlib plotter template.

A template for creating plots with matplotlib.

- Uthpala Herath
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import matplotlib.pyplot as plt
from cycler import cycler
from my_plot import set_size
import pandas as pd
import numpy as np

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
width = 512
fraction = 1

# ------------------------- DATA -------------------------

## Pandas import
# df = pd.read_excel("Excel.xlsx", skiprows=1, usecols="AF:AJ")
# df = df.dropna()
# y = []
# y.append(df["column"])

## Text file import
# with open("file.out", "r") as f:
#     lines = f.readlines()
#     x = [float(line.split()[0]) for line in lines]
#     y = [float(line.split()[1]) for line in lines]


# ------------------------- PLOTTING -------------------------

fig, ax = plt.subplots(1, 1, figsize=set_size(width, fraction))
# ax.set_prop_cycle(line_cycler)

ax.plot(x, y, linestyle=linestyle[0], label="$label$")

# ax.set_prop_cycle(line_cycler)
ax.set_title(r"TITLE")
ax.set_xlabel(r"XLABEL")
ax.set_ylabel(r"YLABEL")
ax.set_ylim([0, 1.5])
# ax.axvline(x=0, color="gray", linestyle="--")
ax.grid()
lgd = ax.legend()

# legend sizes
# lgd.legendHandles[0]._sizes = [30]
# lgd.legendHandles[1]._sizes = [30]

fig.savefig("plot.pdf")
