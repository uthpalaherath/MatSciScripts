#!/usr/bin/env python
#
#  Script to plot and/or compare) band structures from FHI-aims.
#  If nPlots=2, it computes RMSE and optionally plots the difference as well.
#
#  Usage:
#    compare_bands.py N_PLOTS DIRECTORY TITLE ENERGY_OFFSET ...
#                       [yMin yMax] --diffplot
#
#    * N_PLOTS:  number of band structures to plot (1 to 7).
#    * Then for each band structure, supply:
#         DIRECTORY TITLE ENERGY_OFFSET
#    * Finally, optionally specify yMin and yMax for the plot.
#
#  Example:
#    compare_bands.py 2 ./CalcA A 0.0 ./CalcB B 1.0 -8 8 --diffplot
#
# Author: Uthpala Herath
# Based on the aimsplot_compare.py script in FHIaims/utilities

import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.lines as mlines

##############################
# 1) Global Options / Parsing
##############################

# Style file
plt.style.use("~/dotfiles/matplotlib/prb.mplstyle")

should_spline = False
spline_factor = 10
output_x_axis = True
color_array = "rgbycmk"

# Detect a '--diffplot' flag
DO_PLOT_DIFF = False
if len(sys.argv) > 1 and sys.argv[-1] == "--diffplot":
    DO_PLOT_DIFF = True
    sys.argv = sys.argv[:-1]

if len(sys.argv) < 2:
    print(
        """Usage:
  aimsplot_compare.py N_PLOTS DIR TITLE OFFSET ... [yMin yMax] [--diffplot]

  e.g.:
    aimsplot_compare.py 2 CalcA A 0.0 CalcB B 1.0 -5 5
    aimsplot_compare.py 2 CalcA A 0.0 CalcB B 1.0 -5 5 --diffplot
"""
    )
    sys.exit(1)

try:
    nPlots = int(sys.argv[1])
except ValueError:
    print("Error: First arg must be integer (N_PLOTS).")
    sys.exit(1)

if not (1 <= nPlots <= 7):
    print("Error: N_PLOTS must be in [1..7].")
    sys.exit(1)

needed_args = 1 + 3 * nPlots
if len(sys.argv) < needed_args:
    print(
        f"Error: not enough arguments for nPlots={nPlots}. Need >= {needed_args} total."
    )
    sys.exit(1)

# Gather the band-structure data: (directory, title, offset)
directory = []
plotName = []
energy_offset = []

idx = 2
for i in range(nPlots):
    directory.append(sys.argv[idx])
    plotName.append(sys.argv[idx + 1])
    energy_offset.append(float(sys.argv[idx + 2]))
    idx += 3

# Possibly read yMin, yMax
CUSTOM_YLIM = False
ylim_lower, ylim_upper = -20.0, 20.0
if len(sys.argv) >= 1 + 3 * nPlots + 2:
    CUSTOM_YLIM = True
    ylim_lower = float(sys.argv[1 + 3 * nPlots + 1])
    ylim_upper = float(sys.argv[1 + 3 * nPlots + 2])

print(f"Number of band structures: {nPlots}")
for i in range(nPlots):
    print(f"  {i+1}) {directory[i]} -> '{plotName[i]}'  offset={energy_offset[i]}")
if CUSTOM_YLIM:
    print(f"Using custom y-range: [{ylim_lower}, {ylim_upper}]")
else:
    print(f"Using default y-range: [{ylim_lower}, {ylim_upper}]")

##############################
# 2) Data Structures
##############################
band_data = [dict() for _ in range(nPlots)]
max_spin_channel = [1] * nPlots
PLOT_SOC = [False] * nPlots
PLOT_GW = [False] * nPlots
latvec = [[] for _ in range(nPlots)]
rlatvec = [[] for _ in range(nPlots)]
band_segments = [[] for _ in range(nPlots)]
band_totlength = [0.0 for _ in range(nPlots)]

##############################
# 3) Build the main plot
##############################
fig_bands, ax_bands = plt.subplots(figsize=(6, 4))
if output_x_axis:
    ax_bands.axhline(0.0, color="k", linestyle=":")


def nice_label(lab):
    return r"$\Gamma$" if lab.lower() == "gamma" else lab


##############################
# 4) Read geometry.in, control.in for each
##############################
for i in range(nPlots):
    # read geometry
    geo_file = os.path.join(directory[i], "geometry.in")
    lat_temp = []
    with open(geo_file) as gf:
        for line in gf:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            if parts[0] == "lattice_vector":
                lat_temp.append(list(map(float, parts[1:4])))
    if len(lat_temp) != 3:
        raise ValueError(
            "Need exactly 3 lattice_vector lines in geometry.in for " + directory[i]
        )
    latvec[i] = np.array(lat_temp)
    vol = np.dot(latvec[i][0], np.cross(latvec[i][1], latvec[i][2]))
    r1 = 2 * math.pi * np.cross(latvec[i][1], latvec[i][2]) / vol
    r2 = 2 * math.pi * np.cross(latvec[i][2], latvec[i][0]) / vol
    r3 = 2 * math.pi * np.cross(latvec[i][0], latvec[i][1]) / vol
    rlatvec[i] = np.array([r1, r2, r3])

    # read control
    ctrl_file = os.path.join(directory[i], "control.in")
    with open(ctrl_file) as cf:
        for line in cf:
            line = line.split("#")[0].strip()
            if not line:
                continue
            if line.startswith("spin collinear"):
                max_spin_channel[i] = 2
            if any(
                line.startswith(x)
                for x in (
                    "calculate_perturbative_soc",
                    "include_spin_orbit",
                    "include_spin_orbit_sc",
                )
            ):
                PLOT_SOC[i] = True
                max_spin_channel[i] = 1
            if line.startswith("qpe_calc"):
                PLOT_GW[i] = True
            if line.startswith("output band"):
                parts = line.split()
                if len(parts) < 9:
                    raise ValueError("Bad output band line: " + line)
                start = np.array(list(map(float, parts[2:5])))
                end = np.array(list(map(float, parts[5:8])))
                npt = int(parts[8])
                sname = parts[9] if len(parts) > 9 else ""
                ename = parts[10] if len(parts) > 10 else ""
                length = np.linalg.norm(
                    np.dot(rlatvec[i], end) - np.dot(rlatvec[i], start)
                )
                band_segments[i].append((start, end, length, npt, sname, ename))
                band_totlength[i] += length

    if PLOT_SOC[i]:
        max_spin_channel[i] = 1

##############################
# 5) Plot each band structure; store data in band_data
##############################
line_handle_for_calc2 = None  # reference for the second structure's line

for i in range(nPlots):
    gap = band_totlength[i] / 30.0
    xpos = 0.0
    labels = []
    prev_end = None
    ccolor = color_array[i % len(color_array)]

    for seg_index, (start, end, length, npoint, sname, ename) in enumerate(
        band_segments[i], start=1
    ):
        if prev_end is not None and not np.allclose(start, prev_end):
            xpos += gap
        xvals = xpos + np.linspace(0, length, npoint)
        labels.append((xvals[0], sname))
        labels.append((xvals[-1], ename))
        prev_end = end
        xpos = xvals[-1]

        for spin in range(1, max_spin_channel[i] + 1):
            if PLOT_GW[i]:
                fname = os.path.join(directory[i], f"GW_band{spin}{seg_index:03}.out")
            else:
                fname = os.path.join(directory[i], f"band{spin}{seg_index:03}.out")
            # read energies
            energies_list = []
            with open(fname) as fb:
                for ln in fb:
                    w = ln.split()
                    # index,kx,ky,kz, E1_occ,E1, E2_occ,E2,...
                    e_vals = [float(xx) - energy_offset[i] for xx in w[5::2]]
                    energies_list.append(e_vals)
            energies_list = np.array(energies_list)  # shape (npt, nBands)
            xvals_plot = xvals
            # optionally spline
            if should_spline:
                from scipy.interpolate import make_interp_spline

                x_dense = np.linspace(xvals.min(), xvals.max(), npoint * spline_factor)
                nB = energies_list.shape[1]
                new_energies = []
                for b_ in range(nB):
                    spl = make_interp_spline(xvals, energies_list[:, b_], k=3)
                    new_energies.append(spl(x_dense))
                new_energies = np.array(new_energies).T
                xvals_plot = x_dense
                energies_list = new_energies

            band_data[i][(spin, seg_index)] = {
                "xvals": xvals_plot,
                "energies": energies_list,
            }

            # Plot
            n_bands_here = energies_list.shape[1]
            for bI in range(n_bands_here):
                if bI == 0 and seg_index == 1:
                    # label once
                    if max_spin_channel[i] == 2:
                        # up/dn
                        if spin == 1:
                            lbl = plotName[i] + " (up)"
                            (line_,) = ax_bands.plot(
                                xvals_plot,
                                energies_list[:, bI],
                                color=ccolor,
                                label=lbl,
                            )
                        else:
                            lbl = plotName[i] + " (dn)"
                            (line_,) = ax_bands.plot(
                                xvals_plot,
                                energies_list[:, bI],
                                color=ccolor,
                                linestyle="--",
                                label=lbl,
                            )
                    else:
                        lbl = plotName[i]
                        (line_,) = ax_bands.plot(
                            xvals_plot, energies_list[:, bI], color=ccolor, label=lbl
                        )
                        if i == 1:  # second calculation
                            line_handle_for_calc2 = line_
                else:
                    style_ = "--" if (max_spin_channel[i] == 2 and spin == 2) else "-"
                    ax_bands.plot(
                        xvals_plot, energies_list[:, bI], color=ccolor, linestyle=style_
                    )

    # ticks
    usedpos = set()
    for xx, lab in labels:
        if xx not in usedpos:
            ax_bands.axvline(xx, color="k", linestyle=":")
            usedpos.add(xx)
    if labels:
        ax_bands.set_xlim(labels[0][0], labels[-1][0])
        tickx = [l[0] for l in labels]
        tickl = [nice_label(l[1]) for l in labels]
        ax_bands.set_xticks(tickx)
        ax_bands.set_xticklabels(tickl)

if CUSTOM_YLIM:
    ax_bands.set_ylim(ylim_lower, ylim_upper)
else:
    ax_bands.set_ylim(-20, 20)

ax_bands.set_xlabel("Wave Vector")
ax_bands.set_ylabel("Energy (eV)")
legend_main = ax_bands.legend()
if legend_main:
    legend_main.get_frame().set_linewidth(1.0)

##############################
# 6) Compute RMSE if nPlots=2
##############################
def filter_bands_by_energy_range(E, y_min, y_max):
    # E shape: (nK, nB)
    keep = []
    nK, nB = E.shape
    for b in range(nB):
        arr = E[:, b]
        if arr.max() < y_min or arr.min() > y_max:
            continue
        keep.append(b)
    if len(keep) > 0:
        return E[:, keep]
    else:
        return np.zeros((nK, 0))


if nPlots == 2:
    sum_sq = 0.0
    nvals = 0

    # We'll do a naive pass over band_segments[0], spin=1 => match with #1 in second structure
    # ignoring spin=2.  Adjust if you prefer spin=2 or something else.
    for seg_index, segdata in enumerate(band_segments[0], start=1):
        (start, end, length, npoint, sname, ename) = segdata
        keyA = (1, seg_index)
        keyB = (1, seg_index)
        if keyA not in band_data[0] or keyB not in band_data[1]:
            continue
        E1 = band_data[0][keyA]["energies"]
        E2 = band_data[1][keyB]["energies"]
        if E1.shape[0] != E2.shape[0]:
            # mismatch # k-points
            continue
        F1 = filter_bands_by_energy_range(E1, ylim_lower, ylim_upper)
        F2 = filter_bands_by_energy_range(E2, ylim_lower, ylim_upper)
        if F1.size == 0 or F2.size == 0:
            continue
        nb1 = F1.shape[1]
        nb2 = F2.shape[1]
        ncomm = min(nb1, nb2)
        d = F2[:, :ncomm] - F1[:, :ncomm]
        sum_sq += np.sum(d * d)
        nvals += d.size

    if nvals > 0:
        rmse_val = math.sqrt(sum_sq / nvals)
        rmse_str = f"RMSE={rmse_val:.3f} eV"
        print(f"RMSE in range [{ylim_lower}, {ylim_upper}] eV = {rmse_val:.4f} eV.")

        # Create an invisible line2D instance:
        rmse_handle = mlines.Line2D(
            [], [], color="none", marker="", linestyle="none", label=rmse_str
        )

        # Grab the existing legend handles/labels from the main figure:
        handles, labels = ax_bands.get_legend_handles_labels()

        # Append the RMSE entry
        handles.append(rmse_handle)
        labels.append(rmse_str)

        # Re‐draw the legend with this extra line
        ax_bands.legend(handles, labels)

    else:
        print("No overlapping band data => no RMSE to compute.")

##############################
# 7) Optionally do a 2nd figure with difference
##############################
if nPlots == 2 and DO_PLOT_DIFF:
    fig_diff, ax_diff = plt.subplots(figsize=(6, 4))
    ax_diff.axhline(0, color="k", linewidth=1)
    ax_diff.set_xlabel("Wave Vector")
    ax_diff.set_ylabel("Diff (Calc2 - Calc1) [eV]")
    gap_diff = band_totlength[0] / 30.0
    xdiff = 0.0
    prev_end = None
    labels_diff = []
    for seg_index, segdata in enumerate(band_segments[0], start=1):
        (start, end, length, npoint, sname, ename) = segdata
        k1 = (1, seg_index)
        k2 = (1, seg_index)
        if k1 not in band_data[0] or k2 not in band_data[1]:
            continue
        if prev_end is not None and not np.allclose(start, prev_end):
            xdiff += gap_diff
        E1 = band_data[0][k1]["energies"]
        E2 = band_data[1][k2]["energies"]
        if E1.shape[0] != E2.shape[0]:
            continue
        f1 = filter_bands_by_energy_range(E1, ylim_lower, ylim_upper)
        f2 = filter_bands_by_energy_range(E2, ylim_lower, ylim_upper)
        if f1.size == 0 or f2.size == 0:
            continue
        ncommon = min(f1.shape[1], f2.shape[1])
        e1use = f1[:, :ncommon]
        e2use = f2[:, :ncommon]
        diffmat = e2use - e1use

        local_x = np.linspace(0, length, f1.shape[0]) + xdiff
        for b_ in range(ncommon):
            ax_diff.plot(local_x, diffmat[:, b_], color="b")

        labels_diff.append((local_x[0], sname))
        labels_diff.append((local_x[-1], ename))
        prev_end = end
        xdiff = local_x[-1]

    usedx = set()
    for xx, lab in labels_diff:
        if xx not in usedx:
            ax_diff.axvline(xx, color="k", linestyle=":")
            usedx.add(xx)
    if labels_diff:
        tx = [ld[0] for ld in labels_diff]
        tl = [nice_label(ld[1]) for ld in labels_diff]
        ax_diff.set_xticks(tx)
        ax_diff.set_xticklabels(tl)
        ax_diff.set_xlim(tx[0], tx[-1])
    ax_diff.set_title("Difference Plot (Calc2 - Calc1)")
    fig_diff.savefig("bands_diff.pdf")

# Save main figure
fig_bands.savefig("bands_compare.pdf")
