#!/usr/bin/env python
""" Plot data from csv files.

This script generates plots from data retrieved
from csv or similar text based files.

Author: Uthpala Herath

Usage:
    Enter plotdata.py -h for help.

"""
import argparse
import os
import sys
from argparse import RawTextHelpFormatter

import matplotlib.pyplot as plt
import numpy as np


def dataplotter(args):

    infile = args.infile
    xcolumn = args.xcolumn
    ycolumn = args.ycolumn
    headerlines = args.headerlines
    xlabel = args.xlabel
    ylabel = args.ylabel
    title = args.title
    savefig = args.savefig

    file = open(infile, "r")
    if headerlines > 0:
        for i in range(headerlines):  # skipping header lines
            file.readline()
    data = file.readlines()
    file.close()

    xdata = []
    ydata = []
    for i in range(len(data)):
        if xcolumn is not None:
            xdata.append(float(data[i].split()[xcolumn - 1]))
        ydata.append(float(data[i].split()[ycolumn - 1]))
    if xcolumn is not None:
        plt.plot(xdata, ydata)
    else:
        plt.plot(ydata)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(savefig)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument("-infile", type=str, help="Input data file.")
    parser.add_argument(
        "-xcolumn",
        type=int,
        help="Column number for x value. Leave empty to plot only y values.",
        default=None,
    )
    parser.add_argument(
        "-ycolumn", type=int, help="Column number for y value.", default=2
    )
    parser.add_argument(
        "-headerlines", type=int, help="Number of header lines.", default=0
    )
    parser.add_argument("-title", type=str, help="Plot title.", default="X vs Y plot")
    parser.add_argument("-xlabel", type=str, help="xlabel.", default="X")
    parser.add_argument("-ylabel", type=str, help="ylabel.", default="Y")
    parser.add_argument(
        "-savefig", type=str, help="Save figure name.", default="plot.png"
    )
    args = parser.parse_args()
    dataplotter(args)
