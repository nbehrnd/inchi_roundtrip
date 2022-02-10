#!/bin/usr/env python3

# name:    test_smiles_rdkit.py
# author:  nbehrnd@yahoo.com
# license:
# date:    2022-02-10 (YYYY-MM-DD)
# edit:
#

"""Elementar testing about script smiles_rdkit.py."""
import os                                         
from subprocess import getstatusoutput, getoutput    #+end_src
                                                  
import pytest

PROGRAM = str("./smiles_rdkit.py")
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
        success_reader = register[1].split()[0]
    assert success_reader == str("1")

    try:
        for file in ["testfile.smi",
                     "failing_structures.log", "success_structures.log"]:
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
