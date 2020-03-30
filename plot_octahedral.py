#!/usr/bin/env python3
""" This script plots octahedral tilt angle with time."""

import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# number of abinit iterations
iterations = int(sys.argv[1])

# name of system
name = sys.argv[2]

# timestep
dt = float(sys.argv[3])
time = np.arange(0, iterations * dt, dt)
######################################

# creating zero matrices
xcart = np.zeros(iterations)
r19_mag = np.zeros(iterations)
cos_alpha = np.zeros(iterations)
acell1 = np.zeros(iterations)
acell2 = np.zeros(iterations)
acell3 = np.zeros(iterations)
cos_beta = np.zeros(iterations)
cos_gamma = np.zeros(iterations)
etot = np.zeros(iterations)
volume = np.zeros(iterations)
pressure = np.zeros(iterations)

rp1_mod = np.zeros(iterations)
rp2_mod = np.zeros(iterations)
rp3_mod = np.zeros(iterations)

for i in range(iterations):
    # open .out file
    fileop = open("./" + name + str(i) + "/" + name + str(i) + ".out", "r")
    read = fileop.read()
    fileop.close()

    # open .in file to read rprim
    filein = open("./" + name + str(i) + "/" + name + str(i) + ".in", "r")
    readin = filein.read()
    filein.close()

    # constraint values
    acell_init = re.findall(r"acell([0-9E+.\s]*)Bohr", read)[0].split()
    rprim_init = re.findall(r"rprim([0-9E+-.\s]*)", readin)[0].split()
    xcart_init = re.findall(r"xcart([0-9E+-.\s]*)xred", read)[0].split()

    # position vectors of atoms
    r1 = np.array(xcart_init[0:3], float)  # Fe1
    r2 = np.array(xcart_init[3:6], float)
    r3 = np.array(xcart_init[6:9], float)
    r4 = np.array(xcart_init[9:12], float)
    r5 = np.array(xcart_init[12:15], float)
    r6 = np.array(xcart_init[15:18], float)
    r7 = np.array(xcart_init[18:21], float)
    r8 = np.array(xcart_init[21:24], float)
    r9 = np.array(xcart_init[24:27], float)  # O1
    r10 = np.array(xcart_init[27:30], float)
    r11 = np.array(xcart_init[30:33], float)  # O3
    r12 = np.array(xcart_init[33:36], float)
    r13 = np.array(xcart_init[36:39], float)
    r14 = np.array(xcart_init[39:42], float)
    r15 = np.array(xcart_init[42:45], float)
    r16 = np.array(xcart_init[45:48], float)
    r17 = np.array(xcart_init[48:51], float)  # O9

    # cell parameter constraint
    acell1[i] = acell_init[0]
    acell2[i] = acell_init[1]
    acell3[i] = acell_init[2]

    # rprim vectors
    rp1 = np.array([rprim_init[0], rprim_init[3], rprim_init[6]], float)
    rp2 = np.array([rprim_init[1], rprim_init[4], rprim_init[7]], float)
    rp3 = np.array([rprim_init[2], rprim_init[5], rprim_init[8]], float)

    # cell angle constraint
    a = float(acell1[i]) * rp1
    b = float(acell2[i]) * rp2
    c = float(acell3[i]) * rp3

    # Octahedral angle is the angle between the vector of O9 and Fe1 and the b axis
    r_O9Fe1 = r17 - r1
    cos_gamma[i] = (np.dot(r_O9Fe1, b)) / (
        (np.linalg.norm(r_O9Fe1)) * (np.linalg.norm(b))
    )  # angle between O9Fe1 and b vectors

# PLOTTING #

font = {  # 'family' : 'normal',
    # 'weight' : 'bold',
    "size": 15
}

matplotlib.rc("font", **font)

plt.figure(1)
plt.plot(time, cos_gamma)
plt.title("Octahedral tilt angle")
plt.xlabel("Time (fs)")
plt.ylabel("cosine")
# plt.xlim(0,0.05)
plt.tight_layout()
plt.savefig("octahedral_tilting.pdf")
# plt.show()
