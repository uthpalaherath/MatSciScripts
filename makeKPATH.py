#!/usr/bin/env python

"""
KPATH generator.

This script generates a KPOINTS file for a k-path
based on a provided POSCAR file.
It is based on PyProcar.

- Uthpala Herath

"""

import pyprocar

pyprocar.kpath("POSCAR")
