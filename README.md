# wollsey-public
Public stable version of the wollsey stochastic evolution simulator.

## Build instructions
Project includes a cmake build file. Simply clone the repo and
```bash
cd wollsey-public
mkdir bin && cd bin
cmake ..
make
```

### Dependencies
Main program compilation
+ [TCLAP](https://tclap.sourceforge.net) 1.2.1
+ HDF5 1.8.16

Python input/output handling
+ matplotlib 2.1.2
+ h5py 2.7.1
+ click 6.7
+ numpy 1.14.0

## Usage
Running the program takes two files, a system file detailing system parameters such as:
+ Starting cell count and genotype
+ Relevant mutations and their IC50s
+ Reproduction rates and mutation rates

as well as a run file which describes the treatment protocols for simulation. Specifically:
+ A protocol name
+ Drugs
+ A dose timeline for each of the selected drugs

Both of these files are HDF5 files composed using the programs nurse.py and plan.py respectively.
A detailed commented shellscript that generates example files is available in the examples folder.
Because it's quite the process to create the inputs, assembling a script file to create them is strongly recommended.
Briefly
```bash
# names are reused a lot so defined them once for convenience
SYSFN="example-sys.hdf5"
RUNFN="example-run.hdf5"

# create the systemfile
./nurse.py init $SYSFN

# set parameters
./nurse.py setp $SYSFN -i 100 -a 1000000 -z 250000 -k 1.0 -m 500000 -p 1e-7

# add drug parameters
# for brevity, lets look at only two mutants and two drugs
# first, define the interesting mutants
./nurse.py gene $SYSFN -r G250E -r T315I

# then, give their IC50s for drugs
./nurse.py drug $SYSFN -n Imatinib -x 5000 -s 3000 -w 500 -m G250E 3600 -m T315I 9200
./nurse.py drug $SYSFN -n Dasatinib -x 100 -s 80 -w 2 -m G250E 8 -m T315I 140

# set reproduction rates
./nurse.py rate $SYSFN --uni-deathrate 0.006 --uni-birthrate 0.016 --uni-minbirthrate 0.0

# finally define the starting population
./nurse.py cell $SYSFN -c GGGACT 250000

# that's all for the system file
# now to make a protocol

# make the file
./plan.py init $RUNFN

# set simulation repeats
./plan.py setp $RUNFN -r 100 # simulation repeats (per protocol)

# create a simple control protocol
./plan.py newp $RUNFN -n control -l 500000

# And an Imatinib-Dasatinib rotation
./plan.py newp $RUNFN -n IMDA -l 500000
./plan.py edip $RUNFN -n IMDA -d Imatinib -q squarewave 0 2000 400 2000
./plan.py edip $RUNFN -n IMDA -d Dasatinib -q squarewave 1 400 2000 2000
```

All program parameters should also be explained by running the programs with --help

With the inputs set up, we can do a simulation. Running with nohup usually necessary if running on a remote server, as the simulation takes some time and disconnection would otherwise stop it. On a local machine, nohup can be omitted.
```bash
nohup ./wollsey -i example-sys.hdf5 -r example-run.hdf5 -o examtle-out.hdf5 -s 100 -j 10 -v -e constant_branching_process -m single &
```

Output should look something like:
```
Reading parameters
Reading protocols
  Reading protocol IMDA
    Reading drug Imatinib
    Reading drug Dasatinib
  Reading protocol control
Reading radix
Reading protein options
Reading drugs
  Reading drug Dasatinib
  Reading drug Imatinib
Dasatinib	wt 2 snps | 0 8 | 0 140 | 
Imatinib	wt 500 snps | 0 3600 | 0 9200 | 

Reading genotype properties
Unified death rate: 0.006
Unified reproduction rate: 0.016
Unified minimum reproduction rate: 0

Reading starting cells
GGGACT
GGG ACT
G T
0
0
G T
G T
Reading simulation parameters
max_iters	500000
mutation_probability	1e-07
max_cells	1000000
min_cells	100
population_size	250000
spring_force	1
Loaded simulation system

Running simulation(s)
Single run mode
  Running protocol 1/2	IMDA
  Running protocol 2/2	control
Total simulation time: 132 s

```

Finally, there are a few analysis functions built into the program expo.py
For now, lets just compare the wthalf of the protocols
```
./expo.py pcomp example-out.hdf5 -m wthalf --all
```
(Without treatment there is no resistance, so control group is not very interesting)

and look at an individual simulation timeline
```
./expo.py draw example-out.hdf5 -i 0 0 0
```

## Output file structure
Since only a few simple analysis functions are included, here is a detailed description of the output file structure for easier design of custom analysis functions. For convenience and safekeeping, the output file also contains a full copy of both inputs.

HDF5 files have a hierarchal tree-like nature which makes them relatively organized.

| *File structure* |                                               |                                                                      |                                         |       |        |                              |
| ---------------- | --------------------------------------------- | -------------------------------------------------------------------- | --------------------------------------- | ----- | ------ | ---------------------------- |
| Cells            | Counts                                        | *Starting number of cells*                                           |                                         |       |        |                              |
|                  | Genomes                                       | *of these genotypes*                                                 |                                         |       |        |                              |
| Drugs            | Drug name                                     | Single                                                               | *IC50s for single mutations*            |       |        |                              |
|                  |                                               | Unique_IC50s                                                         | *IC50s for compound mutations*          |       |        |                              |
|                  |                                               | Unique_genes                                                         | *Specifiers of these copound mutations* |       |        |                              |
|                  | *One group like that per drug*                |                                                                      |                                         |       |        |                              |
| Gene             | Options                                       | *Set of allowed SNVs*                                                |                                         |       |        |                              |
|                  | Positions                                     | *Residue numbers*                                                    |                                         |       |        |                              |
|                  | Radix                                         | *Number of allowed SNVs*                                             |                                         |       |        |                              |
| Parameters       | *Empty dataset with parameters as attributes* |                                                                      |                                         |       |        |                              |
| Protocol         | Protocol name                                 | Index                                                                | *Drug names for timelines*              |       |        |                              |
|                  |                                               | Record                                                               | *Dose timeline*                         |       |        |                              |
|                  | *One group like that per protocol*            |                                                                      |                                         |       |        |                              |
| Rates            | Genes                                         | *If unique rates for phenotypes are used, this specifies the genes*  |                                         |       |        |                              |
|                  | Reproduction                                  | *and this specifies their reproduction rates*                        |                                         |       |        |                              |
| Record           | 0 *generation number*                         | 0 *protocol number (matches lexicographically sorted protocol list)* | 0 *Run number*                          | Cells | Index  | *List of observed cells*     |
|                  |                                               |                                                                      |                                         |       | Record | *Their population over time* |
|                  |                                               |                                                                      |                                         | Drugs | Index  | *List of used drugs*         |
|                  |                                               |                                                                      | *One group like this per repeat*        |       | Record | *Their dose over time*       |
|                  |                                               | *One group like this per protocol*                                   |                                         |       |        |                              |

Do explore a file with hdfview or similar to get a grasp of it. Some datasets also have attributes, which should have self-explanatory naming.

An implementation concept that it might be worth knowing about is phenotype id's. 
To simplify internal program structure, a cell is usually assigned an id that depends on what SNVs it has. These are derived by treating the list of SNVs as a mixed radix number. 
Suffice to say, in many places where cells are named (like the record-0-0-0-cell-index above) these numbers are used. expo.py has a hid function that translates them to a human-readable format. The specific implementation is in gene.h and gene.cpp.

