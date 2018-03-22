#!/bin/bash

# names are reused a lot so defined them once for convenience
SYSFN="example-sys.hdf5"
RUNFN="example-run.hdf5"

# create the systemfile
./nurse.py init $SYSFN

# set parameters
./nurse.py setp $SYSFN \
           -i 100 \
           -a 1000000 \
           -z 250000 \
           -k 1.0 \
           -m 500000 \
           -p 1e-7

# add drug parameters
# for brevity, lets look at only two mutants and two drugs
# first, define the interesting mutants
./nurse.py gene $SYSFN \
           -r G250E \
           -r T315A \
           -r T315I

# then, give their IC50s for drugs
# you can define more drugs than you intend to use in a specific protocol.
./nurse.py drug $SYSFN \
           -n Imatinib \
           -x 5000 \
           -s 3000 \
           -w 500 \
           -m G250E 3600 \
           -m T315I 9200
./nurse.py drug $SYSFN \
           -n Dasatinib \
           -x 100 \
           -s 80 \
           -w 2 -m G250E 8 \
           -m T315I 140

# set reproduction rates
./nurse.py rate $SYSFN \
           --uni-deathrate 0.006 \
           --uni-birthrate 0.016 \
           --uni-minbirthrate 0.0

# finally define the starting population
./nurse.py cell $SYSFN \
           -c GGGACT 250000

# that's all for the system file
# now to make a protocol

# make the file
./plan.py init $RUNFN

# set simulation repeats
./plan.py setp $RUNFN \
          -r 100 # simulation repeats (per protocol)

# create a simple control protocol
./plan.py newp $RUNFN \
          -n control \
          -l 500000

# And an Imatinib-Dasatinib rotation
./plan.py newp $RUNFN \
          -n IMDA \
          -l 500000
./plan.py edip $RUNFN \
          -n IMDA \
          -d Imatinib \
          -q squarewave 0 2000 400 2000
./plan.py edip $RUNFN \
          -n IMDA \
          -d Dasatinib \
          -q squarewave 1 400 2000 2000
