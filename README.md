Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark

The  structuresel.py is a Python script to: i) generate an ensemble of structure from profasi, ii) select a reference 
structure on the base of the calculated Rg.

It needs to have Marshsequence.py for the Rg calculation in the same folder
Marshsequence.py alculates Rh based on N, charge, his-tag, and prolines from Marsh et al. XXX



### REQUIREMENT BEFORE RUNNING ###

A virtual environment for python might be required
The script needs MD_Analysis version XX installed
and profasi version 1.4.9

the path to profasi installation is hardcoded at lines 93, and 108, they need to be replaced before usage


### FOR CBL GROUP MEMBERS only ####
to activate the virtual environment on bioinfo01: 

. /usr/local/envs/py27/bin/activate

### HOW TO RUN ###

detailed options are available with

python strucuresel2.py -h

examples are reported in the folder case_studies 

a standard running line (in interactive mode) is:

python structuresel2.py -s AAAA -i -NACE -CNH2 -nstruct 10000 

where s is the sequence with one-letter aminoacid code, -NACE and -CNH2 takes care of the capping mode at N- and C-terminal of the ID region, -nstruct is the number of cycles for the run and structures to be generated 

if you need to perform a dry run (i.e. without running profasi) use -n

python structuresel.py -s AEPADMRPEIWIAQELRRIGDEFNAYYA -i -NACE -CNH2 -nstruct 10000 -n


### FOR CBL GROUP MEMBERS ONLY ####
to activate the virtual environment: 

. /usr/local/envs/py27/bin/activate


### NOTES ###

We tested it with profasi version 1.4.9 - more recent profasi versions might require updates
Please get in touch with us if this is the case

The script to use is the newest version structuresel2.py where also the possibility to select the random seed for the run
is included. This is fundamental for reproducing an already existing run.

If we do not define a random seed,  profasi will use MersenneTwister (a default random number generator) and it stores the info in the file in n0/random_number_state and in the log file


