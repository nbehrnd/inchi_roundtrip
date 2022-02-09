#!/usr/bin/env python3

# name:    inchi_sdf.py
# author:  nbehrnd@yahoo.com
# license: MIT 2022
# date:    2022-02-01 (YYYY-MM-DD)
# edit:    2022-02-09 (YYYY-MM-DD)

"""Monitor a round trip .sdf -> InChI -> .sdf, InChI executable.

Starting from a multi-model .sdf non-zero coordinates, this script relies on
modules of the Python standard library and InChI trusts' reference executable
(version 1.06).  This script, InChI trust's executable (with the added
executable bit), and the data to process are expected to reside in the same
folder.

Anticipated input:  multi-model .sdf (e.g., by DataWarrior)
Anticipated output: a report about structures passing/failing this test.

The criterion for passing the round trip is an invariant InChI string."""

import argparse
import os
import subprocess


def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        usage="""Check round-trip .sdf -> INCHI -> .sdf with InChI v1.06.

The anticipated input file is a multi-model .sdf file to process (the file
extension does not matter).  Keep the inchi-1 executable (v 1.05) for Linux by
InChI trust (add the executable bit) in the same folder as this script and
provide with OpenBabel's Python libraries.

If an entry's non-standard InChI (fixed H-layer) prior and after the round trip
match each other, the structure enters file success_structures.log.  Else, the
entry is reported in file failing_structures.log.""")


    parser.add_argument("source_file",
                        metavar="FILE",
                        help="Input .sdf file containing a list molecules.")

    return parser.parse_args()


def model_lister(input_file=""):
    """Return the model data of a .sdf as a listing.

This counts (x + 1) entries and looses the terminal `$$$$` string of each model."""
    all_model_data = ""
    register = []

    with open(input_file, mode="r") as source:
        for line in source:
            all_model_data += "".join([str(line)])
            
  #  del register[0]
    register = all_model_data.split("$$$$\n")
    return register

def assign_primary_inchi(input_model=""):
    """Assign the primary InChI to the original datum.

    Input:  the primary .sdf string
    Output: the primary InChI string."""
    register = []
    primary_inchi = ""

    with open("testfile.sdf", mode="w") as newfile:
        newfile.write(input_model)
        newfile.write("$$$$\n")
        
    process=subprocess.Popen(["./inchi-1", "-FixedH", "-AuxNone",
                              "-NoLabels",
                              "testfile.sdf", "primary_inchi.txt"],
                             shell=False)
    process.communicate()

    try:
        with open("primary_inchi.txt", mode="r") as source:
            register = source.readlines()
            primary_inchi = str(register[0]).strip()
    except OSError:
        print("No access to 'primary_inchi.txt'.")

    for file in ["testfile.sdf",
                 "testfile.sdf.prb", "testfile.sdf.log",
                 "primary_inchi.txt"]:
        try:
            os.remove(file)
        except OSError:
            print(f"Remove of file '{file}' was unsuccessful.")

    return primary_inchi

def assign_inchi_auxiliary(inchi_string=""):
    """Generate an auxiliary file for a structure recovery.

    Input:  the primary InChI string (cf. function assign_primary_inchi)
    Output: a temporary auxiliary.txt"""
    with open("testfile.txt", mode="w") as newfile:
        newfile.write(inchi_string)

    process=subprocess.Popen(["./inchi-1", "-InChI2Struct",
                              "testfile.txt", "auxiliary.txt"],
                              shell=False)
    process.communicate()

    for file in ["testfile.txt",
                 "testfile.txt.log", "testfile.txt.prb"]:
        try:
            os.remove(file)
        except OSError:
            print(f"Remove of file '{file}' was unsuccessful.")


def generate_inchi_sdf():
    """Based on an auxiliary, let InChI generate a .sdf.

    Input:  auxiliary.txt
    Output: output.sdf"""

    process=subprocess.Popen(["./inchi-1", "-OutputSDF",
                              "auxiliary.txt", "output.sdf"],
                             shell=False)
    process.communicate()

    for file in ["auxiliary.txt",
                 "auxiliary.txt.log", "auxiliary.txt.prb"]:
        try:
            os.remove(file)
        except OSError:
            print(f"Remove of file '{file}' was unsuccessful.")


def correct_secondary_sdf():
    """As already filed, there is a superfluous heading line in InChI's .sdf"""
    register = []
    try:
        with open("output.sdf", mode="r") as source:
            for line in source:
                register.append(str(line).rstrip())
        del register[0]

        with open("output.sdf", mode="w") as newfile:
            for entry in register:
                newfile.write(f"{entry}\n")
    except OSError:
        print(f"correction of InChI's .sdf failed.")

def assign_secondary_inchi():
    """Assign InChI on the newly generated .sdf.

    Input:   output.sdf
    Output:  secondary InChI string"""
    register = []
    secondary_inchi = ""
    
    process=subprocess.Popen(["./inchi-1", "-fixedH", "-AuxNone",
                              "-NoLabels",
                              "output.sdf", "secondary_inchi.txt"],
                              shell=False)
    process.communicate()

    try:
        with open("secondary_inchi.txt", mode="r") as source:
            register = source.readlines()
            secondary_inchi = str(register[0]).strip()
    except OSError:
        print(f"Assignment secondary InChI failed.")

    try:
        for file in ["output.sdf",
                    "output.sdf.log", "output.sdf.prb",
                    "secondary_inchi.txt"]:
            os.remove(file)
    except OSError:
        print(f"Removal of file '{file}' failed.")

    return secondary_inchi

def main():
    """Join the functions.

Reading a .sdf file which may contain one datum, or multiple model data,
the round trip's InChI strings per model are compared with each other.  If
they match pairwise, the round trip was successful; else, it failed.  The
script is going to write two new .sdf files according to these categories."""
    success = []
    failing = []
    counter = int(1)

    args = get_args()
    input_file = args.source_file
    list_of_models = model_lister(input_file=input_file)

    for entry in list_of_models[:-1]:
        primary_inchi = ""
        secondary_inchi = ""

        primary_inchi = assign_primary_inchi(input_model=entry)

        assign_inchi_auxiliary(inchi_string=primary_inchi)
        generate_inchi_sdf()
        correct_secondary_sdf()

        secondary_inchi = assign_secondary_inchi()

        if primary_inchi == secondary_inchi:
            success.append(f"{counter}\t{primary_inchi}\t{secondary_inchi}")
        else:
            failing.append(f"{counter}\t{primary_inchi}\t{secondary_inchi}")
        counter += int(1)

    print("\n---- ----\n")
    print("Brief report:")
    print(f"success structures: {len(success)}")
    with open("success_structures.log", mode="w") as newfile:
        newfile.write("successful round trips\n")
        newfile.write("\t".join(['counter',
                                 'primary InChI', 'secondary InChI\n']))
        for entry in success:
            newfile.write(f"{entry}\n")
        newfile.write("END")

    print(f"failing structures: {len(failing)}")
    with open("failing_structures.log", mode="w") as newfile:
        newfile.write("failing round trips:\n")
        newfile.write("\t".join(['counter',
                                 'primary InChI', 'secondary InChI\n']))
        for entry in failing:
            newfile.write(f"\n{entry}")
        newfile.write("\nEND")

    print("\nCheck file 'success_structures.log' and 'failing_structures.log'.")


if __name__ == "__main__":
    main()
