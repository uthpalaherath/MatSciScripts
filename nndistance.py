#!/usr/bin/env python

"""
Nearest Neighbour Distance

This script reads a poscar and prints a list of the minimum distance
to the nearest neibour for an atom, for all the atoms in the structure.

- Uthpala Herath

"""

import pychemia
from heapq import nsmallest

st = pychemia.code.vasp.read_poscar("POSCAR")
ca = pychemia.analysis.cluster.ClusterAnalysis(st)
da = ca.distance_matrix()


def second_smallest(numbers):
    return nsmallest(2, numbers)[-1]


for i in range(len(da)):
    print(second_smallest(da[:, i]))
