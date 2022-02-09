#!/usr/bin/env python3

# name:    smiles_obsdf.py
# author:  nbehrnd@yahoo.com
# license: MIT 2022
# date:    2022-01-29 (YYYY-MM-DD)
# edit:    2022-02-09 (YYYY-MM-DD)

"""Monitor a round trip SMILES -> .sdf -> INCHI -> .sdf -> SMILES.

The aim is to monitor how reliable the reconstruction of .sdf from an InChI
string actually is.  It is assumed that a successful round trip (SMILES at start
matching SMILES at the end) requires InChI with fixed H-layer to account for
tautomerism.  However, it is not evident if this suffices for any organic
structure submitted as this; axial chirality (the motif of 1,1'-biphenyl,
TADDOL, BINAP, etc.) possibly present a difficulty here.

Anticipated input: a list of SMILES (e.g. by a DataWarrior library)
Anticipated output: a report about SMILES passing/failing this test.

This script relays some work to the nonstandard libraries of OpenBabel and
RDKit.  The assignment of InChI as well as the regeneration of .sdf requires the
reference InChI executable distributed by InChI trust (v. 1.06); here, the
version for Linux is anticipated."""

import argparse
import os
import subprocess

import openbabel
from openbabel import pybel
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        usage="""Check round-trip SMILES -> .sdf -> INCHI -> .sdf -> SMILES.

The anticipated input file is a listing of SMILES to process (the file
extension does not matter).  Keep the inchi-1 executable (v 1.06) for
Linux by InChI trust (add the executable bit) in the same folder as
this script and provide with OpenBabel's Python libraries.

If an entry's canonical SMILES prior and after the round trip match
each other, the structure enters file success_structures.log.  Else,
the SMILES prior and after the round trip are recorded in the file
failing_structures.log.  The criterion currently deployed is OpenBabel's
canonical SMILES about the intermediate .sdf written.""")


    parser.add_argument("source_file",
                        metavar="FILE",
                        help="Input file containing a list of SMILES strings.")

    return parser.parse_args()

def split(input_file=""):
    """Read the SMILES into a list"""
    input_list = []
    with open(input_file, mode="r") as newfile:
        for entry in newfile:
            input_list.append(str(entry).strip())

    return input_list

def smiles2obabel(initial_smiles=""):
    """Convert SMILES into OpenBabel's canonical SMILES."""
    mol = pybel.readstring("smi", initial_smiles)
    obabel_smiles = str(mol.write("can"))

    return obabel_smiles


def smiles2rdkit(initial_smiles=""):
    """Convert SMILES into RDKit's SMILES."""
    mol = Chem.MolFromSmiles(initial_smiles)
    rdkit_smiles = Chem.MolToSmiles(mol,isomericsmiles=False)

    return rdkit_smiles

def sdf_obabel(raw_smiles=""):
    """Generate a .sdf with OpenBabel."""
    mol = pybel.readstring("smi", raw_smiles)
    mol.make3D()
    molecule = mol.write("sdf")

    with open("test_file.sdf", mode="w") as newfile:
        newfile.write(molecule)

def sdf_rdkit(raw_smiles=""):
    """Generate a .sdf with RDKit."""
    mol = Chem.MolFromSmiles(raw_smiles)
    with_hydrogens = Chem.AddHs(mol)
    molecule = Chem.MolToMolBlock(with_hydrogens)

    with open("test_file.sdf", mode="w") as newfile:
        newfile.write(molecule)

def assign_inchi(initial_sdf=""):
    """Assign InChI on the initial .sdf.

    Input:   test_file.sdf
    Output:  inchi.txt"""
    process=subprocess.Popen(["./inchi-1",  "-fixedH",
                              "test_file.sdf", "inchi.txt"],
                              shell=False)
    process.communicate()

    for file in os.listdir("."):
        if (file.endswith(".sdf") or
            file.endswith(".log") or file.endswith(".prb")):
            os.remove(file)


def assign_inchi_auxiliary():
    """Generate an auxiliary for a structure recovery.

    Input:  inchi.txt
    Output: auxiliary.txt"""
    process=subprocess.Popen(["./inchi-1", "-InChI2Struct",
                              "inchi.txt", "auxiliary.txt"],
                              shell=False)
    process.communicate()

    for file in os.listdir("."):
        if (file.endswith(".log") or file.endswith(".prb")):
            os.remove(file)
    os.remove("inchi.txt")


def generate_inchi_sdf():
    """Let InChI generate a .sdf.

    Input:  auxiliary.txt
    Output: output.sdf"""
    process=subprocess.Popen(["./inchi-1", "-OutputSDF",
                              "auxiliary.txt", "output.sdf"],
                             shell=False)
    process.communicate()

    for file in os.listdir("."):
        if (file.endswith(".log") or file.endswith(".prb")):
            os.remove(file)
    os.remove("auxiliary.txt")

def trim_sdf_file():
    """Remove the superfluous leading lines inchi-1 wrote in the .sdf."""
    register = []

    with open("output.sdf", mode="r") as newfile:
        register = newfile.readlines()
        register = register[1:]

    with open("output.sdf", mode="w") as newfile:
        for line in register:
            newfile.write(f"{line}")


def obabel_newsmiles():
    """Assign the canonical SMILES by OpenBabel on the new structure."""
    new_smiles = ""
    for mol in pybel.readfile("sdf", "output.sdf"):
        new_smiles = mol.write("can")

    return new_smiles


def rdkit_smiles():
    """Assign the SMILES by RDKit on the new structure."""
    new_smiles = ""
    mol = Chem.MolFromMolFile("output.sdf")
    new_smiles = Chem.MolToSmiles(mol, isomericsmiles=False)

    return new_smiles

def main():
    """Join the functions."""
    args = get_args()
    input_file = args.source_file

    success = []
    failing = []

    listed = split(input_file)
    for entry in listed:
        raw_smiles = ""
        raw_smiles = str(smiles2obabel(entry))
        raw_smiles = raw_smiles.split()[0]

        sdf_obabel(raw_smiles)

        assign_inchi("test_file.sdf")
        assign_inchi_auxiliary()
        generate_inchi_sdf()

        trim_sdf_file()

        new_smiles = str(obabel_newsmiles()).strip().split()[0]

        if str(raw_smiles) == str(new_smiles).split()[0]:
            success.append(raw_smiles)
        else:
            retain = "\t".join([raw_smiles, new_smiles])
            failing.append(retain)
    os.remove("output.sdf")

    print("\n---- ----\n")
    print("Brief report:")
    print(f"success structures: {len(success)}")
    with open("success_structures.log", mode="w") as newfile:
        for entry in success:
            newfile.write(f"{entry}\n")
        newfile.write("END")

    print(f"failing structures: {len(failing)}")
    with open("failing_structures.log", mode="w") as newfile:
        newfile.write("SMILES (prior)\tSMILES (after) round trip:\n")
        for entry in failing:
            newfile.write(f"{entry}\n")
        newfile.write("\nEND")

    print("\nCheck file 'success_structures.log' and 'failing_structures.log'.")


if __name__ == "__main__":
    main()
