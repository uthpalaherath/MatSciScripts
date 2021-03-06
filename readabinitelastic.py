#!/usr/bin/env python

"""
Parses relaxed elatic tensor from abinit output file and
prints in units of GPa.
"""

import re
import sys
import numpy as np

fi = open(sys.argv[1], "r")
data = fi.read()
fi.close()

ET = re.findall(
    "Elastic\s*Tensor\s*\(relaxed\s*ion\)\s*[\sA-Za-z:\d\^()]*\n([-\s0-9.]*)", data
)
ET2 = np.array([float(x) for x in ET[0].split()])
elastic_tensor = ET2.reshape(6, 6)

# in GPa
elastic_tensor = elastic_tensor * 100

row = elastic_tensor.shape[0]
col = elastic_tensor.shape[1]
for i in range(row):
    for j in range(col):
        print("{:>10.4f} ".format(elastic_tensor[i, j]), end=" ")
        if j == (col - 1):
            print(" ")
