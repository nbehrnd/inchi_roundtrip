#!/usr/bin/bash

# name:   sdf2inchi2inchi.sh
# author: nbehrnd@yahoo.com
# date:   2022-02-11 (YYYY-MM-DD)
# edit:

# Compact representation of round trip checks moderated with Python 3
# cf. https://github.com/nbehrnd/inchi_roundtrip
#
# Written for Linux Debian 12/bookworm (branch testing), requires Python 3,
# bash, awk, OpenBabel (version 3.1.1) -- all as provided by the Debian's
# repositories; and InChI trust's reference exectutable (version 1.06)
# with exectutable bit added (once `chmod u+x inchi-1`).

# generation of the primary model:
# obabel -:"COC" -h --gen3d -O model_a.sdf  # dimethyl ether (InChI strings match)
# obabel -:"OC1=NC=CC=C1" -h --gen3d -O model_a.sdf  # 2-hydroxy pyridine (InChI strings match)
# obabel -:"O=C1C=CC=CN1" -h --gen3d -O model_a.sdf  # 2-pyridone (InChI strings match)

# obabel -:"C[C@@H](O)CC" -h --gen3d -O model_a.sdf  # (R)-butan-2-ol (InChI strings do not match)
obabel -:"C/C=C/C" -h --gen3d -O model_a.sdf  # (E)-but-2-ene (InChI strings do not match)

# assign a primary InChI string:
./inchi-1 model_a.sdf inchi_a.txt -AuxNone -NoLabels -fixedH -SUCF on
rm model_a.sdf.log
rm model_a.sdf.prb

# assign an InChI auxillary:
./inchi-1 inchi_a.txt auxillary.txt -InChI2Struct 
rm inchi_a.txt.log
rm inchi_a.txt.prb

# let inchi-1 write a .sdf:
./inchi-1 auxillary.txt model_b.sdf -OutputSDF
rm auxillary.txt.log
rm auxillary.txt.prb

# remove the first line in inchi-1's .sdf written:
awk 'NR>=2 {print}' model_b.sdf > short.sdf
mv short.sdf model_b.sdf

# assign a secondary InChI string:
./inchi-1 model_b.sdf inchi_b.txt -AuxNone -NoLabels -fixedH -SUCF on
rm model_b.sdf.log
rm model_b.sdf.prb

# show the two InChI strings anticipated to match each other:
cat model_a.txt
cat model_b.txt
