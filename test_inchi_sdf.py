#!/usr/bin/env python3

# name:    test_inchi_sdf.py
# author:  nbehrnd@yahoo.com
# license: MIT, 2022
# date:    2022-07-02 (YYYY-MM-DD)
# edit:    2022-09-02 (YYYY-MM-DD)
#

"""Identify systematic errors in inchi_sdf.py in a round trip.

Script inchi_sdf.py aims to monitor the round trip of .sdf for which InChI
trust's reference executable is going to assign a primary InChI string.  Based
on this reduced representation, InChI trust's reference executable then
reconstructs a .sdf for which again a -- then secondary -- InChI string is
assigned.  Thus, except for the generation of the initial .sdf file, this is
round trip independent of interaction by OpenBabel, or RDKit.  This script
test_inchi_sdf.py shall detect problems in the implementation of the intended
round trip, well ahead of working on libraries of molecules.

Deployment with pytest for Python 3 which may -- depending on the Linux
distribution used -- may be accessible by pytest, or pytest-3

pytest -v test_inchi_sdf.py

The execution depends on the simultaneous presence of this script, InChI trust's
reference binary for Linux (version 1.06), and inchi_sdf.py."""
import os
from subprocess import getstatusoutput, getoutput

import pytest

PROGRAM = str("./inchi_sdf.py")
TFILE = str("testfile.sdf")

def write_testfile(model=""):
    """Provide a file with a single input structure."""
    with open(TFILE, mode="w") as newfile:
        newfile.write(str(model))


def make_tester(model=""):
    """Provide the frame to perform a mono-model test."""
    register = []
    check = ""
    write_testfile(model=model)

    test = getoutput(f"python3 {PROGRAM} {TFILE}") # export

    try:
        with open("success_structures.log", mode="r") as source:
            register = source.readlines()
        check = str(register[2])
    except OSError:
        print("File 'success_structures.log' inaccessible.  Exit.")

    assert str(check[0]) == str("1")


def test_dimethylether():
    """Check a structure not prone to tautomerism."""
    model=str("""
 OpenBabel020722

  9  8  0  0  0  0  0  0  0  0999 V2000
    1.0496    0.0176   -0.0708 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4706    0.0409   -0.0805 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9681    1.2713   -0.5890 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7233   -0.9468    0.3277 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6579    0.1267   -1.0868 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6579    0.8122    0.5717 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.0606    1.2342   -0.5737 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6377    2.1060    0.0369 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6377    1.4205   -1.6216 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  2  3  1  0  0  0  0
  3  7  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
M  END
$$$$
""")
    make_tester(model=model)


def test_2hydroxypyridine():
    """Check on a structure subject to tautomerism, 1/2."""
    model=str("""
 OpenBabel020922

 12 12  0  0  0  0  0  0  0  0999 V2000
    1.1010    0.0362    0.2292 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4532    0.1277    0.0818 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0800    1.3588   -0.0232 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4624    1.3765   -0.1859 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1572    0.1724   -0.2356 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4357   -1.0046   -0.1155 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0966   -1.0509    0.0435 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7330    0.9279    0.1589 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5246    2.2874    0.0199 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.9939    2.3198   -0.2738 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.2335    0.1504   -0.3648 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.9303   -1.9708   -0.1450 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  8  1  0  0  0  0
  2  3  1  0  0  0  0
  2  7  2  0  0  0  0
  3  4  2  0  0  0  0
  3  9  1  0  0  0  0
  4  5  1  0  0  0  0
  4 10  1  0  0  0  0
  5  6  2  0  0  0  0
  5 11  1  0  0  0  0
  6  7  1  0  0  0  0
  6 12  1  0  0  0  0
M  END
$$$$
""")


def test_2pyridone():
    """Check on a structure subject to tautomerism, 2/2."""
    model=str("""
 OpenBabel020922

 12 12  0  0  0  0  0  0  0  0999 V2000
    2.1606    0.0524    0.0004 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9356    0.0285    0.0010 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1384    1.2834    0.0009 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2014    1.2181   -0.0004 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8769   -0.0594   -0.0011 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1608   -1.1905    0.0010 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2101   -1.1399    0.0021 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.6873    2.2152    0.0021 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8036    2.1203   -0.0002 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9616   -0.0773   -0.0031 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6167   -2.1750    0.0015 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7533   -1.9942    0.0027 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  2  7  1  0  0  0  0
  3  4  2  0  0  0  0
  3  8  1  0  0  0  0
  4  5  1  0  0  0  0
  4  9  1  0  0  0  0
  5  6  2  0  0  0  0
  5 10  1  0  0  0  0
  6  7  1  0  0  0  0
  6 11  1  0  0  0  0
  7 12  1  0  0  0  0
M  END
$$$$
""")
