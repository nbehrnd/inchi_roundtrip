#!/bin/usr/env python3

# name:    test_smiles_obsdf.py
# author:  nbehrnd@yahoo.com
# license: GPL v3, 2022
# date:    2022-02-07 (YYYY-MM-DD)
# edit:
#

"""Provide basic testing about script smiles_obsdf.py's round trip.

For the initial run with the set 100_smiles.txt, 57/100 entries did not pass
successfully the round trip Openbabel SMILES -> OpenBabel .sdf -> InChI string
-> InChI .sdf -> OpenBabel SMILES.  To identify systematic errors in own
programming early, this script tests the processing with pytest when calling

pytest -v test_smiles_obsdf.py

Proper execution of this test script depends on the presence of smiles_obsdf.py
in the same directory as this script, test_smiles_obsdf.py.  It equally requires
the InChI trust reference executable for Linux, an working installation of the
non-standard Python libraries about OpenBabel and Pytest (for Python3).
Depending on the Linux distribution used, Pytest (for Python3) might be called
by either pytest, or explicit pytest-3."""
import os
from subprocess import getstatusoutput, getoutput

import pytest

PROGRAM = str("./smiles_obsdf.py")
TFILE = str("testfile.smi")


def write_testfile(SMILES=""):
    """Provide a file with the input structure."""
    with open(TFILE, mode="w") as newfile:
        newfile.write(str(SMILES))


def make_tester(structure=""):
    """Provide the frame to perform tests on varying SMILES strings."""
    smiles = ""
    smiles = str(structure)
    write_testfile(SMILES=smiles)
    assert os.path.isfile(TFILE)

    test = getoutput(f"python3 {PROGRAM} {TFILE}")
    assert os.path.isfile("failing_structures.log")
    assert os.path.isfile("success_structures.log")

    register = []
    success_reader = ""
    with open("success_structures.log", mode="r") as source:
        register = source.readlines()
    success_reader = str(register[1]).split()[0]
    assert success_reader == str("1")

    try:
        for file in ["testfile.smi",
                     "failing_structures.log", "sucess_structures.log"]:
            os.remove(file)
    except OSError:
        print(f"Remove of file '{file}' failed.")


def test_dimethylether():
    """Check a structure not prone to tautomerism."""
    make_tester(structure="COC")


def test_2hydroxypyridine():
    """Check a structure prone to tautomerism, 1/2."""
    make_tester(structure="Oc1ccccn1")


def test_2pyridone():
    """Check a structure prone to tautomerism, 2/2."""
    make_tester(structure="O=c1cccc[nH]1")
