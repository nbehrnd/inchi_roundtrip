#!/usr/bin/env python3

# name:    smiles_rdkit.py
# author:  nbehrnd@yahoo.com
# license: MIT 2022
# date:    2022-01-09 (YYYY-MM-DD)
# edit:

"""Monitor a round trip SMILES > .sdf > InChI > .sdf > SMILES (rdkit).

Starting with a list of SMILES, RDKit will assign a unique primary SMILES
and generate a .sdf for an InChI assignment by InChI trust's reference
executable.  InChI trust's reference executable recreates a new .sdf for
RDKit's assignment of a secondary SMILES.

The round trip is successful if the two SMILES strings match each other.

For a successful execution, deposit this script with inchi-1 and the
.sdf to process in the same folder.  Provide inchi-1 the executable bit.
RDKit is not part of Python's standard library."""

import argparse
import os
import subprocess

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
this script and provide with RDKit's Python libraries.

If an entry's isomeric SMILES prior and after the round trip match
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

def smiles2rdkit(initial_smiles=""):
    """Convert SMILES into RDKit's SMILES."""
    mol = Chem.MolFromSmiles(initial_smiles)
    rdkit_smiles = Chem.MolToSmiles(mol)

    return rdkit_smiles

def sdf_rdkit(raw_smiles=""):
    """Generate a .sdf with RDKit."""
    mol = Chem.MolFromSmiles(raw_smiles)
    with_hydrogens = Chem.AddHs(mol)
    molecule = Chem.MolToMolBlock(with_hydrogens)

    with open("test_file.sdf", mode="w") as newfile:
        newfile.write(molecule)

def assign_inchi():
    """Assign InChI on the initial .sdf.

    Input:   test_file.sdf
    Output:  inchi.txt"""
    register = []
    inchi = ""
    
    process=subprocess.Popen(["./inchi-1",  "-fixedH",
                              "test_file.sdf", "inchi.txt"],
                              shell=False)
    process.communicate()

    with open("inchi.txt", mode="r") as source:
        register = source.readlines()
        inchi = str(register[0])
        
    try:
        for file in ["test_file.sdf",
                     "test_file.sdf.log", "test_file.sdf.prb"]:
            os.remove(file)
    except OSError:
        print(f"Remove of '{file}' failed.")

    return inchi


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

def rdkit_smiles():
    """Assign the SMILES by RDKit on the new structure."""
    new_smiles = ""
    mol = Chem.MolFromMolFile("output.sdf")
    new_smiles = Chem.MolToSmiles(mol)

    return new_smiles

def main():
    """Join the functions."""
    args = get_args()
    input_file = args.source_file

    success = []
    failing = []
    counter = int(1)
    
    listed = split(input_file)
    for entry in listed:
        raw_smiles = ""
        raw_smiles = str(smiles2rdkit(initial_smiles=entry))  #.strip()))
        raw_smiles = raw_smiles.split()[0]

        sdf_rdkit(raw_smiles=raw_smiles)
        primary_inchi = assign_inchi()
        assign_inchi_auxiliary()
        generate_inchi_sdf()
        trim_sdf_file()
        secundary_smiles = rdkit_smiles()

        if raw_smiles == secundary_smiles:
            success.append(f"{counter}\t{raw_smiles}\t{secundary_smiles}")
        else:
            failing.append(f"{counter}\t{raw_smiles}\t{secundary_smiles}")
        counter += int(1)
        os.remove("output.sdf")

    print("\n---- ----\n")
    print("Brief report:")
    print(f"success structures: {len(success)}")
    with open("success_structures.log", mode="w") as newfile:
        newfile.write("SMILES (prior)\tSMILES (after) round trip:\n")
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
