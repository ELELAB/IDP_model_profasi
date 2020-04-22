#!/usr/bin/env python

__author__ = "Mads Nygaard"
__date__ = "20200422" # Edited by Burcu Aykac Fas

import MDAnalysis
import subprocess
from numpy import std
from Marshsequence import rhcalculator, rhToRg, sequencechecker


def rglist(uni, v=False):
    rglist = []
    CAs = uni.select_atoms("name CA")
    if map(int, MDAnalysis.__version__.split(".")) < [0, 16, 0]:
        print("WARNING: MDAnalysis can be slow for v. < 0.16.0") 
    for idx, ts in enumerate(uni.trajectory):
        rglist.append(CAs.atoms.radius_of_gyration())
        if (v is True) and (idx % 500 == 0):
            print("i= %s" % idx)
    return rglist


def printframe(uni, frameidx, name="output.pdb"):
    protein = uni.select_atoms("all")
    with MDAnalysis.Writer(name, protein.n_atoms) as W:           
        uni.trajectory[frameidx]
        W.write(protein)
          

def selclosest(alist, val, k):
    from heapq import nsmallest
    return nsmallest(k, enumerate(alist), key=lambda x: abs(x[1]-val))


def bulkyness(Aaa):
    """Bulkyness of amino acids
    Ref:  J. Theor. Biol. 21:170-201(1968)."""
    bulkdir = {"ALA": 11.500, "ARG": 14.280, "ASN": 12.820, 
               "ASP": 11.680, "CYS": 13.460, "GLN": 14.450, 
               "GLU": 13.570, "GLY":  3.400, "HIS": 13.690, 
               "ILE": 21.400, "LEU": 21.400, "LYS": 15.710, 
               "MET": 16.250, "PHE": 19.800, "PRO": 17.430, 
               "SER":  9.470, "THR": 15.770, "TRP": 21.670, 
               "TYR": 18.030, "VAL": 21.570}
    try:
        return bulkdir[Aaa.upper()]
    except:
        raise Exception("Input %s not long form of AA" % Aaa)


def pepBulk(sequence):
    """Returns total bulkyness number of peptide"""
    bnum = 0
    for A in sequence:
        bnum += bulkyness(aatolong(A))
    return bnum


def aatoshort(Aaa):
    aatable = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", 
               "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I", 
               "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", 
               "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", 
               "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"}
    try:
        return aatable[Aaa.upper()]
    except:
        raise Exception("Input %s not long form of AA" % Aaa)


def aatolong(A):
    aatable = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 
               'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 
               'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 
               'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 
               'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}
    try:
        return aatable[A.upper()]
    except:
        raise Exception("Input %s not short form of AA" % A)


def getseqfrompdb(filename):
    uni = MDAnalysis.Universe(filename)
    sequence = []    
    for idx, n in enumerate(uni.residues):
        sequence.append(aatoshort(uni.residues[idx].name))
    return "".join(sequence)


def basicMCrun(sequence, N_ACE=False, C_NH2=False, temp=300, nstruct=100):     
    path = "/data/tools_scripts_repository/PROFASI_muninn/profasi/app/bin/BasicMCRun"
    Nterm = ''
    Cterm = ''
    if N_ACE is True:
        Nterm = "ACE"
    if C_NH2 is True:
        Cterm = "NH2"
    sequence2 = r'"< ' + Nterm + r' * ' + sequence + r' * ' + Cterm + r' >"' 
    command1 = "".join([path, ' -ac', ' 1 ', str(sequence2),
                        ' -t ', str(temp), ' -ncyc ', str(nstruct), 
                        ' -nrt ', str(10), ' -lcyc ', str(-1), 
                        ' -iconf ', str(1), ' -random_number_seed ', str(args.random_number_seed) ])
    print command1
    with open("logfile.txt", "w") as f:
        subprocess.call(command1, stdout=f, shell=True)
        subprocess.call("/data/tools_scripts_repository/PROFASI_muninn/profasi/app/bin/extract_props -sf -pdb n0/traj", stdout=f, shell=True)


def printreport(rglist, sequence):
    '''Prints report, returns touple of marshRg, avg, stddev '''
    avg = (sum(rglist)/len(rglist))
    rh = rhcalculator(sequence)
    marshRg = rhToRg(rh, len(sequence))
    stddev = std(rglist)
    print "\n%s structures created using PROFASI for a peptide" % len(rglist)
    print "with Marsh et al. predicted Rh = %.5g Rg = %.5g" % (rh, marshRg)
    print "| MIN    | AVG    | MAX    |"  
    print "| %.5g | %.5g | %.5g |" % (min(rglist), avg, max(rglist))
    print "\nStd. dev. = %.5g" % stddev
    if len(sequence) < 50:
        print "\nNOTE: Marsh eq. not accurate, less than 50 res in seq.\n"
    return (marshRg, avg, stddev)


def getchoicefromuser(N):
    is_valid = 0
    while not is_valid:
        choice = raw_input('Enter your choice [1-%s] : ' % N)
        if choice.lower() in ["exit", "q"]:
            exit(0)
        try:
            choice = int(choice)
            if choice in range(1, N+1):
                is_valid = 1
            else:
                print "Not valid number, try again"
        except ValueError, e:
            print ("'%s' is not a valid integer." % e.args[0].split(": ")[1])
    return choice


def getnumfromuser():
    is_valid = 0
    while not is_valid:
        try:
            choice = float(raw_input('Enter Rg: '))
            is_valid = 1
        except ValueError, e:
            print ("'%s' is not a valid number." % e.args[0].split(": ")[1])
    return choice
    

if __name__ == "__main__":
    import argparse
    import sys
    choices = {"Marsh": 1, "AvgStd": 2, "Avg": 3}

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""This script creates 'nstruct' random structures using PROFASI and selects one based on several different parameters and outputs it to output.pdb. PROFASI is run in the location of the script.\nIt is possible to run on a self created pdb-ensemble using -n and -f flags.\nNOTE: The script does not delete or overwrite already run PROFASI ensembles in the folder.\nNOTE: The script uses Marshsequence.py which should also be located in its folder. \n\nExample usage:\n    python structuresel.py -s ACDEFGHILMNOPQRSTVY -i \n""")
    parser.add_argument("-s", required=True, metavar="ACDEFGHIK...", help="Sequence", type=str)     
    parser.add_argument("-i", action="store_true", required=False, help="Interactive mode")    
    parser.add_argument("-sel", required=False, type=str, help="Select the structure based on", choices=choices.keys())
    parser.add_argument("-size", required=False, type=float, help="Rg based structure selection", metavar="14.1")
    parser.add_argument("-n", required=False, action="store_true", help="Run without profasi. Set enseble to select from with -f flag.") 
    parser.add_argument("-f", required=False, type=str, default="frames.pdb", help="Ensemble input filename.")
    parser.add_argument("-o", required=False, type=str, default="output", help="Output file name.")
    parser.add_argument("-nstruct", help="Number of structures generated in PROFASI", metavar="100", required=False, type=int, default=100)
    parser.add_argument("-NACE", help="Add neutral ACE group in the N-term.", required=False, action="store_true")
    parser.add_argument("-CNH2", help="Add neutral NH2 group in the C-term.", required=False, action="store_true")
    parser.add_argument("-N_out", help="Output N structures", metavar="N", type=int, default=1, required=False)
    parser.add_argument("-random_number_seed", help="Random Seed Number", type=float, required=False)
    
    args = parser.parse_args()
        
    sequ = args.s
    sequencechecker(sequ)
    if args.n is False:
        basicMCrun(sequ, nstruct=args.nstruct, N_ACE=args.NACE, C_NH2=args.CNH2)
    try:
        print("Reading frames...")
        uni = MDAnalysis.Universe(args.f)
    except IOError:
        sys.stderr.write("Error: cannot find file %s\n" % args.f)
        exit(1)
    print("Calculating Rg")
    rgvals = rglist(uni, v=True)
    print("Done!")
    marshRg, avg, stddev = printreport(rgvals, sequ)
    choice = 0
    if args.sel:
        choice = choices[args.sel] 
    if args.size:
        choice = 4
    if args.i:
        print "Select a suitable structure based on:"
        print "1. Marsh equation   Rg = %.5g" % (marshRg)
        print "2. Avg + stddev     Rg = %.5g" % (avg+stddev)
        print "3. Avg              Rg = %.5g" % (avg)    
        print "4. Own input        Rg = ????"    
        choice = getchoicefromuser(4)
    if not choice:
        choice is False
    elif choice == 1:
        selval = marshRg
    elif choice == 2:
        selval = avg+stddev
    elif choice == 3:
        selval = avg
    elif choice == 4:
        if args.i:
            selval = getnumfromuser()
        else:
            selval = args.size 
    if choice:
        print "Choosing structure closest to Rg = %.5g" % (selval)
        if args.N_out > len(rgvals):
            print "Number of structures to output exceeds number of generated structures, will output 1"
            args.N_out = 1
        idx_rgval_list = selclosest(rgvals, selval, args.N_out)
        if selval > max(rgvals):
            print "\nNOTE: The choosen Rg is larger than all of the generated structures\n"
        elif selval < min(rgvals):
            print "\nNOTE: The choosen Rg is smaller than all of the generated structures\n"
        for idx1, rgval in idx_rgval_list:
            outname = args.o+str(idx1)+".pdb" 
            printframe(uni=uni, frameidx=idx1, name=outname)
            tempuni = MDAnalysis.Universe(outname)
            selsize = rglist(tempuni)[0]
            print "Frame no. %s selected from %s" % (idx1, args.f)
            print "with an Rg of %.5g" % selsize 
    else:
        print "No selection was made, strange. Did you run with -i or select with -sel?"
    print "\n"
