#!/usr/bin/env python

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

infile = sys.argv[1]
column_number = int(sys.argv[2])


file = open(infile, "r")
data = file.readlines()
file.close()


ydata = []
for i in range(len(data)):
    ydata.append(float(data[i].split()[column_number]))

plt.plot(ydata)
plt.show()
