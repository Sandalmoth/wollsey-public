#!/bin/bash

# names are reused a lot so defined them once for convenience
SYSFN="example-sys.hdf5"
RUNFN="example-run.hdf5"

# create the systemfile
./nurse.py init $SYSFN

# set parameters
# there are a few different population size parameters here as there are multiple simulation engines
# that handle population size balancing differently.
# these are
# constant_branching_process, which keeps the mean constant by adjusting death rate (described in paper)
# limited_branching_process, which applies the spring outside of a defined interval but allows variation within
# branching_process, which doesn't limit population size at all. Be careful with overflow on this one.
# -i [integer > 0] minimum number of cells for limited_branching_process simulation engine
# -a [integer > 0] maximum number of cells for limited_branching_process simulation engine
# -z [integer > 0] population size for constant_branching_process simulation engine
# -k [float > 0]   spring constant describing how strongly constant branching process tends towards mean
# -m [integer > 0] max iterations in simulation
# -p [float > 0]   mutation probability
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
# the given residues have to be in numerical order
# -r [mutation code on form A123BCDE] defines both the wildtype sequence (before number) and allowed SNVs (after number)
./nurse.py gene $SYSFN \
           -r G250E \
           -r T315AI

# then, give their IC50s for drugs
# you can define more drugs than you intend to use in a specific protocol.
# -n the drug name (used later when defining the protocols)
# -w IC50 for the wildtype (unit is irrelevant, but keep it consistent, at least within one drug)
# -m IC50 for a mutation. Can be a SNV or a compound mutation as shown.
# All SNVs included must have an IC50! It is optional for compound mutations
# (in which case it is assumed to be the maximum IC50 of its constituent SNVs)
./nurse.py drug $SYSFN \
           -n Imatinib \
           -w 500 \
           -m G250E 3600 \
           -m T315A 890 \
           -m T315I 9200 \
           -m 'G250E T315I' 10000
./nurse.py drug $SYSFN \
           -n Dasatinib \
           -w 2 -m G250E 8 \
           -m T315A 110 \
           -m T315I 140 \
           -m 'G250E T315I' 300

# set reproduction rates
# --uni-deathrate sets a deathrate, (which is ignored in the constant_branching_process simulation engine)
# --uni-birthrate sets a birthrate for all cells
# --uni-minbirthrate sets a minimum birthrate. Cells cannot divide slower than this regardless of drug dose
# it is also possible to specify reproduction rates for all SNVs individually.
# In that case, use
# -p, and -q to specify wildtype birthrate and minimum birthrate
# and -b and -m to specify birthrate and minimum birthrate for individual mutations.
# remember to specify all SNVs
# However, individual birthrate simulations are only experimentally implemented, be cautious of results.
./nurse.py rate $SYSFN \
           --uni-deathrate 0.006 \
           --uni-birthrate 0.016 \
           --uni-minbirthrate 0.0

# finally define the starting population
# -c adds a population of starting cells, (can be used multiple times for a mixed starting pouplation)
# -c [genetic code for chosen residues (gene section)] [number]
# in this case, we have the codons for G250 and T315 in order.
./nurse.py cell $SYSFN \
           -c GGGACT 250000

# that's all for the system file
# now to make a protocol

# make the file
./plan.py init $RUNFN

# set simulation repeats
# since it's stochastic with a lot of variation, it is generally necessary to run simulations a
# large number of times (1000+) to gather adequate staitstics.
./plan.py setp $RUNFN \
          -r 100 # simulation repeats (per protocol)

# create a simple control protocol
# -n [Protocol name (string, [a-ZA-Z_0-9])]
# -l [length of protocol. Strongly recommend that this matches max iterations from before]
./plan.py newp $RUNFN \
          -n control \
          -l 500000

# And an Imatinib-Dasatinib rotation
# this is more complicated, most importantly it involves the edip command, which edits a raw blank protocol
# -n [name of protocol to edit]
# -d [drug  to edit]
# -q [edit function, a variety of equations and parameters available]. Chech ./plan.py edip --help for details
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
