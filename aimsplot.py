#!/usr/bin/env python
#
#  Script to plot band structure and DOS calculated with FHI-aims.
#  Requires the control.in/geometry.in files as well as the output
#  of the calculation to be in the directory from which the script is called.
#
#  This plot will be created in the same directory as a file aimsplot.png if all goes well.
#
#  This version requires Python3. For Python2 (not recommended) go to the Python2 folder.
#
#  You will also need to install the corresponding matplotlib package for python.
#
#  Specifically, the script allows to create a plot of
#
#    - Energy band structures                  ("output band ...")
#
#    - Density of states                       ("output dos ...")
#
#    - Species-projected densities of states   ("output species_proj_dos ...")
#      as well as, optionally, a decomposition of the DOS into angular momentum components
#
#    - and, alternatively, the tetrahedron-integrated versions of DOS and species-projected DOS
#      (much better resolution)
#
#  This script can be called simply as "aimsplot.py", in which case a default range will be used for
#  for the energy range covered on the y axis.
#
#  There are several options that allow one to customize the energy range, type of output,
#  legend placement etc.
#
#    aimsplot.py --help
#
#  provides an overview of the available options.
#
#  For example,
#
#    aimsplot.py --Emin -20. --Emax 10. --legend_x 1. --legend_y 0.2
#
#  will customize both the y axis range shown, as well as the placement of the legend in
#  the graph.
#
#  To achieve labelling of the special points along the band structure plot,
#  add two arguments to the "output band"
#  command in the control.in, using the following syntax:
#
#    output band <start> <end> <npoints> <starting_point_name> <ending_point_name>
#
#  Example: To plot a band with 20 points from Gamma to half way along one of the
#           reciprocal lattice vectors, write (in control.in)
#
#    output band 0.0 0.0 0.0 0.5 0.0 0.0 20 Gamma <End_point_name>
#
#  It is important to note that the graph can easily be further customized by editing
#  this script, using the documentation available online for matplotlib.

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from optparse import OptionParser
from os.path import exists

###########
# OPTIONS #
###########

print_resolution = 250  # The DPI used for printing out images
default_line_width = (
    1  # Change the line width of plotted bands and k-vectors, 1 is default
)
font_size = 12  # Change the font size.  12 is the default.
should_spline = (
    False  # Turn on spline interpolation for band structures NOT VERY WELL TESTED!
)
output_x_axis = True  # Whether to output the x-axis (e.g. the e=0 line) or not
spline_factor = 10  # If spline interpolation turned on, the sampling factor (1 is the original grid)
maxdos_output = -1  # The maximum value of the DOS axis (a.k.a. x-axis) in the DOS
# For zero or negative values, the script will use its default value, the maximum
# value for the DOS in the energy window read in
########################

print("Plotting bands for FHI-aims!")
print("============================")
print()

print("Reading lattice vectors from geometry.in ...")

matplotlib.rcParams["lines.linewidth"] = default_line_width

latvec = []

CUSTOM_YLIM = False
FERMI_OFFSET = False

energy_offset = 0.0
ylim_lower = -100000000.0
ylim_upper = -100000000.0

parser = OptionParser()

parser.description = "aimsplot.py produces a graph of available band structures and densities of states in the directory from which it was called."

parser.add_option(
    "--Emin",
    dest="Emin",
    type="float",
    default=-100000000.0,
    help="Minimum energy value on the y axis of the plot(s).",
)

parser.add_option(
    "--Emax",
    dest="Emax",
    type="float",
    default=-100000000.0,
    help="Maximum energy value on the y axis of the plot(s).",
)

parser.add_option(
    "--Eoffset",
    dest="Eoffset",
    type="float",
    default=0.0,
    help="An energy offset that may be aplied to shift the values on the y axis of the plot(s).",
)

parser.add_option(
    "--legend_x",
    dest="legend_x",
    type="float",
    default=1.5,
    help="A x offset on the canvas that allows one to shift the legend horizontally.",
)

parser.add_option(
    "--legend_y",
    dest="legend_y",
    type="float",
    default=1.0,
    help="A y offset on the canvas that allows one to shift the legend vertically.",
)

parser.add_option(
    "--no_legend",
    action="store_true",
    dest="no_legend",
    help="No legend will be printed into the plot.",
)

parser.add_option(
    "--show_l_components",
    action="store_true",
    dest="show_l",
    help="L components of the species_dos will be plotted if available.",
)

(options, args) = parser.parse_args()

ylim_lower = options.Emin
ylim_upper = options.Emax
energy_offset = options.Eoffset
legend_x_offset = options.legend_x
legend_y_offset = options.legend_y

if options.no_legend:
    SHOW_LEGEND = False
else:
    SHOW_LEGEND = True

if options.show_l:
    SHOW_L = True
else:
    SHOW_L = False

# Verify implicitly which options were set

if ylim_lower == -100000000.0 and ylim_upper == -100000000.0:
    CUSTOM_YLIM = False
    ylim_lower = -20.0
    ylim_upper = 5.0
else:
    CUSTOM_YLIM = True
    if ylim_lower == -100000000.0:
        ylim_lower = -20.0
    if ylim_upper == -100000000.0:
        ylim_upper = 5.0

if energy_offset != 0.0:
    FERMI_OFFSET = True

for line in open("geometry.in"):
    line = line.split("#")[0]
    words = line.split()
    if len(words) == 0:
        continue
    if words[0] == "lattice_vector":
        if len(words) != 4:
            raise Exception("geometry.in: Syntax error in line '" + line + "'")
        latvec += [np.array(list(map(float, words[1:4])))]

if len(latvec) != 3:
    raise Exception("geometry.in: Must contain exactly 3 lattice vectors")

latvec = np.asarray(latvec)

print("Lattice vectors:")
for i in range(3):
    print(latvec[i, :])
print()

# Calculate reciprocal lattice vectors
rlatvec = []
volume = np.dot(latvec[0, :], np.cross(latvec[1, :], latvec[2, :]))
rlatvec.append(np.array(2 * np.pi * np.cross(latvec[1, :], latvec[2, :]) / volume))
rlatvec.append(np.array(2 * np.pi * np.cross(latvec[2, :], latvec[0, :]) / volume))
rlatvec.append(np.array(2 * np.pi * np.cross(latvec[0, :], latvec[1, :]) / volume))
rlatvec = np.asarray(rlatvec)

# rlatvec = inv(latvec) Old way to calculate lattice vectors
print("Reciprocal lattice vectors:")
for i in range(3):
    print(rlatvec[i, :])
print()

########################

print("Reading information from control.in ...")

PLOT_BANDS = False
PLOT_DOS = False
PLOT_DOS_TETRAHEDRON = False
PLOT_DOS_SPECIES = False
PLOT_DOS_SPECIES_TETRAHEDRON = False
PLOT_DOS_ATOM = False
PLOT_DOS_ATOM_TETRAHEDRON = False
PLOT_SOC = False  # This is needed because there will only be one "spin" channel output,
# but collinear spin may (or may not) be turned on, so the "spin
# collinear" setting needs to be overridden
PLOT_DOS_REVERSED = False

species = []

max_spin_channel = 1
band_segments = []
band_totlength = 0.0  # total length of all band segments

for line in open("control.in"):
    words = line.split("#")[0].split()
    nline = " ".join(words)

    if nline.startswith("spin collinear") and not PLOT_SOC:
        max_spin_channel = 2

    if (
        nline.startswith("calculate_perturbative_soc")
        or nline.startswith("include_spin_orbit")
        or nline.startswith("include_spin_orbit_sc")
    ):
        PLOT_SOC = True
        max_spin_channel = 1

    if nline.startswith("output band "):
        if len(words) < 9 or len(words) > 11:
            raise Exception("control.in: Syntax error in line '" + line + "'")
        PLOT_BANDS = True
        start = np.array(list(map(float, words[2:5])))
        end = np.array(list(map(float, words[5:8])))
        length = np.linalg.norm(np.dot(rlatvec, end) - np.dot(rlatvec, start))
        band_totlength += length
        npoint = int(words[8])
        startname = ""
        endname = ""
        if len(words) > 9:
            startname = words[9]
        if len(words) > 10:
            endname = words[10]
        band_segments += [(start, end, length, npoint, startname, endname)]

    if nline.startswith("output dos "):
        PLOT_DOS = True

    if nline.startswith("output dos_tetrahedron"):
        PLOT_DOS_TETRAHEDRON = True

    if nline.startswith("output species_proj_dos "):
        PLOT_DOS_SPECIES = True

    if nline.startswith("output species_proj_dos_tetrahedron"):
        PLOT_DOS_SPECIES_TETRAHEDRON = True

    if nline.startswith("output atom_proj_dos "):
        PLOT_DOS_ATOM = True

    if nline.startswith("output atom_proj_dos_tetrahedron"):
        PLOT_DOS_ATOM_TETRAHEDRON = True

    if nline.startswith("species"):
        if len(words) != 2:
            raise Exception("control.in: Syntax error in line '" + line + "'")
        species += [words[1]]

#######################

if PLOT_SOC:
    max_spin_channel = 1

if PLOT_BANDS and (
    PLOT_DOS or PLOT_DOS_TETRAHEDRON or PLOT_DOS_SPECIES or PLOT_DOS_SPECIES_TETRAHEDRON
):
    ax_bands = plt.axes([0.1, 0.1, 0.6, 0.8])
    ax_dos = plt.axes([0.72, 0.1, 0.18, 0.8], sharey=ax_bands)
    ax_dos.set_title("DOS")
    plt.setp(ax_dos.get_yticklabels(), visible=False)
    ax_bands.set_ylabel("E [eV]")
    PLOT_DOS_REVERSED = True
elif PLOT_BANDS:
    ax_bands = plt.subplot(1, 1, 1)
elif (
    PLOT_DOS or PLOT_DOS_TETRAHEDRON or PLOT_DOS_SPECIES or PLOT_DOS_SPECIES_TETRAHEDRON
):
    ax_dos = plt.subplot(1, 1, 1)
    ax_dos.set_title("DOS")
    PLOT_DOS_REVERSED = False

#######################

if PLOT_BANDS:
    print("Plotting %i band segments..." % len(band_segments))

    if output_x_axis:
        ax_bands.axhline(0, color=(1.0, 0.0, 0.0), linestyle=":")

    prev_end = band_segments[0][0]
    distance = (
        band_totlength / 30.0
    )  # distance between line segments that do not coincide

    iband = 0
    xpos = 0.0
    labels = [(0.0, band_segments[0][4])]

    for start, end, length, npoint, startname, endname in band_segments:
        iband += 1

        if any(start != prev_end):
            xpos += distance
            labels += [(xpos, startname)]

        xvals = xpos + np.linspace(0, length, npoint)
        xpos = xvals[-1]

        labels += [(xpos, endname)]

        prev_end = end
        prev_endname = endname

        for spin in range(1, max_spin_channel + 1):
            # Check if either bandxxx or GW_bandxxx files exit.
            # UKH
            fname1 = "band%i%03i.out" % (spin, iband)
            fname2 = "GW_band%i%03i.out" % (spin, iband)
            idx = []
            kvec = []
            band_energies = []
            band_occupations = []

            if exists(fname1):
                for line in open(fname1):
                    words = line.split()
                    idx += [int(words[0])]
                    kvec += [list(map(float, words[1:4]))]
                    band_occupations += [list(map(float, words[4::2]))]
                    band_energies += [list(map(float, words[5::2]))]
                    # Apply energy offset if specified to all band energies just read in
                    band_energies[-1] = [x - energy_offset for x in band_energies[-1]]

            elif exists(fname2):
                for line in open(fname2):
                    words = line.split()
                    idx += [int(words[0])]
                    kvec += [list(map(float, words[1:4]))]
                    band_occupations += [list(map(float, words[4::2]))]
                    band_energies += [list(map(float, words[5::2]))]
                    # Apply energy offset if specified to all band energies just read in
                    band_energies[-1] = [x - energy_offset for x in band_energies[-1]]

            else:
                print("Neither bandxxx or GW_bandxxx files found!")
                sys.exit(0)

            assert (npoint) == len(idx)
            band_energies = np.asarray(band_energies)
            # Now perform spline interpolation on band structure if requested
            # if should_spline == True:
            #    xvals_smooth = np.linspace(xvals.min(),xvals.max(),spline_factor*len(xvals) ) # Interpolated x axis for spline smoothing
            #    new_band_energies = []
            #    for b in range(band_energies.shape[1]): # Spline every band, one by one
            #        new_band_energies.append(spline(xvals, band_energies[:,b], xvals_smooth))
            #    band_energies = np.asarray(new_band_energies).transpose() # recombine the bands back into the original data format
            #    xvals = xvals_smooth # and use the interpolated x axis
            for b in range(band_energies.shape[1]):
                ax_bands.plot(xvals, band_energies[:, b], color=" br"[spin])

    tickx = []
    tickl = []
    for xpos, l in labels:
        ax_bands.axvline(xpos, color="k", linestyle=":")
        tickx += [xpos]
        if len(l) > 1:
            if l == "Gamma":
                l = "$\\" + l + "$"
        tickl += [l]
    for x, l in zip(tickx, tickl):
        print("| %8.3f %s" % (x, repr(l)))

    ax_bands.set_xlim(labels[0][0], labels[-1][0])
    ax_bands.set_xticks(tickx)
    ax_bands.set_xticklabels(tickl)

#######################


def smoothdos(dos):
    dos = np.asarray(dos)
    # JW: Smoothing is actually done within FHI-aims...
    return dos


if PLOT_DOS or PLOT_DOS_TETRAHEDRON or PLOT_DOS_SPECIES or PLOT_DOS_SPECIES_TETRAHEDRON:

    print("Plotting DOS")

else:
    ax_bands.set_ylim(
        ylim_lower, ylim_upper
    )  # just some random default -- definitely better than the full range including core bands


if PLOT_DOS_SPECIES or PLOT_DOS_SPECIES_TETRAHEDRON:
    spinstrs = [""]
    if max_spin_channel == 2:
        spinstrs = ["_spin_up", "_spin_dn"]

    species_energy = []
    tdos = []
    ldos = []
    maxdos = 0.0

    for s in species:
        val_s = []
        for ss in spinstrs:
            if PLOT_DOS_SPECIES:
                f = open(s + "_l_proj_dos" + ss + ".dat")
            else:
                f = open(s + "_l_proj_dos" + ss + "_tetrahedron.dat")
            f.readline()
            f.readline()
            mu = float(f.readline().split()[-2])
            f.readline()
            val_ss = []
            for line in f:
                val_ss += [list(map(float, line.split()))]
            val_s += [val_ss]
        val_s = np.asarray(val_s).transpose(1, 0, 2)
        # Here val_s is a NumPy data structures, so to apply offset
        # we don't need to use list comprehension
        species_energy += [val_s[:, :, 0] - energy_offset]
        tdos += [smoothdos(val_s[:, :, 1])]
        ldos += [smoothdos(val_s[:, :, 2:])]
        maxdos = max(maxdos, tdos[-1].max())
    for e in species_energy:
        for i in range(e.shape[1]):
            assert all(e[:, i] == species_energy[0][:, 0])
    species_energy = species_energy[0][:, 0]

if PLOT_DOS or PLOT_DOS_TETRAHEDRON:
    if PLOT_DOS:
        f = open("KS_DOS_total.dat")
    else:
        f = open("KS_DOS_total_tetrahedron.dat")
    f.readline()
    mu = float(f.readline().split()[-2])
    f.readline()
    energy = []
    dos = []
    if max_spin_channel == 1:
        for line in f:
            if not line.startswith("#"):
                e, d = line.split()
                energy += [float(e)]
                energy[-1] = energy[-1] - energy_offset
                dos += [(1 * float(d),)]
    else:
        for line in f:
            if not line.startswith("#"):
                e, d1, d2 = line.split()
                energy += [float(e)]
                # Apply energy offset if specified to all DOS energies just read in
                energy[-1] = energy[-1] - energy_offset
                dos += [(float(d1), float(d2))]
    energy = np.asarray(energy)
    dos = smoothdos(dos)
    maxdos = dos.max()

    spinsgn = [1.0]
    if max_spin_channel == 2:
        spinsgn = [1.0, -1.0]

    if PLOT_DOS_REVERSED:
        ax_dos.axhline(0, color="k", ls="--")
        ax_dos.axvline(0, color=(0.5, 0.5, 0.5))

        if PLOT_DOS_SPECIES or PLOT_DOS_SPECIES_TETRAHEDRON:
            for sp in range(len(species)):
                for ispin in range(max_spin_channel):
                    ax_dos.plot(
                        tdos[sp][:, ispin] * spinsgn[ispin],
                        species_energy,
                        linestyle="-",
                        label="%s %s" % (species[sp], ["up", "down"][ispin]),
                    )
                    if SHOW_L:
                        for l in range(ldos[sp].shape[2]):
                            ax_dos.plot(
                                ldos[sp][:, ispin, l] * spinsgn[ispin],
                                species_energy,
                                linestyle="--",
                                label="%s (l=%i) %s"
                                % (species[sp], l, ["up", "down"][ispin]),
                            )

        if PLOT_DOS or PLOT_DOS_TETRAHEDRON:
            for ispin in range(max_spin_channel):
                ax_dos.plot(dos[:, ispin] * spinsgn[ispin], energy, color="kr"[ispin])

        if maxdos_output > 0:
            # If the user has specified a maximum DOS value, use it
            ax_dos.set_xlim(
                np.array([min(spinsgn[-1], 0.0) - 0.05, 1.00]) * maxdos_output
            )
        else:
            # Otherwise use the maximum DOS value read in
            ax_dos.set_xlim(np.array([min(spinsgn[-1], 0.0) - 0.05, 1.05]) * maxdos)
        if CUSTOM_YLIM:
            ax_dos.set_ylim(ylim_lower, ylim_upper)
        else:
            if PLOT_DOS or PLOT_DOS_TETRAHEDRON:
                ax_dos.set_ylim(energy[0], energy[-1])
            else:
                if PLOT_DOS_SPECIES or PLOT_DOS_SPECIES_TETRAHEDRON:
                    ax_dos.set_ylim(species_energy[0], species_energy[-1])
    else:
        ax_dos.axvline(0, color="k", ls="--")
        ax_dos.axhline(0, color=(0.5, 0.5, 0.5))

        if PLOT_DOS_SPECIES or PLOT_DOS_SPECIES_TETRAHEDRON:
            for sp in range(len(species)):
                for ispin in range(max_spin_channel):
                    ax_dos.plot(
                        energy,
                        tdos[sp][:, ispin] * spinsgn[ispin],
                        color="br"[ispin],
                        linestyle="-",
                        label="%s %s" % (species[sp], ["up", "down"][ispin]),
                    )
                    for l in range(ldos[sp].shape[2]):
                        ax_dos.plot(
                            energy,
                            ldos[sp][:, ispin, l] * spinsgn[ispin],
                            color="br"[ispin],
                            linestyle="--",
                            label="%s (l=%i) %s"
                            % (species[sp], l, ["up", "down"][ispin]),
                        )

        if PLOT_DOS or PLOT_DOS_TETRAHEDRON:
            for ispin in range(max_spin_channel):
                ax_dos.plot(energy, dos[:, ispin] * spinsgn[ispin], color="br"[ispin])

        ax_dos.set_xlim(energy[0], energy[-1])
        if CUSTOM_YLIM:
            ax_dos.set_xlim(ylim_lower, ylim_upper)
        else:
            if maxdos_output > 0:
                # If the user has specified a maximum DOS value, use that instead
                ax_dos.set_ylim(
                    np.array([min(spinsgn[-1], 0.0) - 0.05, 1.00]) * maxdos_output
                )
            else:
                # Otherwise use the maximum DOS value read in
                ax_dos.set_ylim(np.array([min(spinsgn[-1], 0.0) - 0.05, 1.05]) * maxdos)
        ax_dos.set_xlabel(r"$\varepsilon - \mu$ (eV)")

    if PLOT_DOS_SPECIES or PLOT_DOS_SPECIES_TETRAHEDRON:
        if SHOW_LEGEND:
            ax_dos.legend(bbox_to_anchor=(legend_x_offset, legend_y_offset))

#######################

matplotlib.rcParams["savefig.dpi"] = print_resolution
matplotlib.rcParams["font.size"] = font_size

print()
print(
    "The resolution for saving figures is set to ",
    matplotlib.rcParams["savefig.dpi"],
    " dpi.",
)

if should_spline:
    print()
    print(
        "Spine interpolation has been used on the band structure, with an interpolation factor of ",
        spline_factor,
    )
    print(
        "You should check this band structure against the un-interpolated version, as splining may cause some small artifacts not present in the original band structure."
    )


def on_q_exit(event):
    if event.key == "q":
        sys.exit(0)


plt.connect("key_press_event", on_q_exit)
plt.savefig("aimsplot.png")
#plt.show()
