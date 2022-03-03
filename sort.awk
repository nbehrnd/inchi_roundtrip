#!/usr/bin/awk -f

# name:    sort.awk
# author:  nbehrnd@yahoo.com
# license: GPLv2, 2022
# date:    2022-03-03 (YYYY-MM-DD)
# edit:

# Jmol's report assigns the round trip either successful (i.e., both
# SMILES and InChI string are invariant), or failing.  By calling this
# awk script by
#
# awk -f sort.awk report.log
#
# `failing.txt` and `success.txt` written provide input for a visual
# inspection by cdkdepict.[1]  The name of this script's input file
# need not be `report.log`.
#
# [1] https://www.simolecule.com/cdkdepict/depict.html

BEGIN {print "Script sorts entries in Jmol's log according to their status."};

{if($2 == "failing"){print $3 "\n" $4 "\n" > "failing.txt"}};
{if($2 == "success"){print $3 > "success.txt"}};

END {print "See CDKdepict: https://www.simolecule.com/cdkdepict/depict.html"};
