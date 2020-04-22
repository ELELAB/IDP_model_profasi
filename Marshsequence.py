#!/usr/bin/env python

__author__ = "Mads Nygaard"
__date__ = "20170614"


def rhcalculator(seq):  # Published equation
    """Calculates Rh based on N, charge, his-tag, and prolines from Marsh"""
    sequence = list(seq)
    k_A = 1.24
    k_B = 0.904
    k_C = 0.00759
    k_D = 0.963
    k_S_his = 0.901
    k_R_0 = 2.49
    k_nu = 0.509
    seqlen, nopro, ishis = sequencechecker(sequence)
    if ishis is True:
        k_S_his2 = k_S_his
    else:
        k_S_his2 = 1.0
    P_pro = float(nopro)/float(seqlen) 
    a_Charge = abs(seqcharge(sequence))
    R_h_calc = (k_A * P_pro + k_B) * (k_C * a_Charge + k_D) * k_S_his2 * k_R_0 * seqlen**k_nu
    return R_h_calc


def seqcharge(sequence):  # Type: list -> int 
    """Outputs the net charge of a sequence at pH = 7.4"""
    AApKa = {"D": -1, "E": -1, "H": 1, "K": 1, "R": 1}
    totalcharge = 0
    for key in AApKa:
        if AApKa[key] > 0:
            tempcharge = 1*sequence.count(key)
        elif AApKa[key] < 0:
            tempcharge = -1*sequence.count(key)
        totalcharge += tempcharge
    return float(totalcharge)


def sequencechecker(sequence):  # Takes sequence in list, CAPS. Outputs (seqlength, no of Prolines,histagstatus) 
    histag = ("H", "H", "H", "H", "H", "H")
    RealAA = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
    if set(RealAA).issuperset(sequence):
        prolines = float(sequence.count("P"))
    else:
        raise Exception("Something wrong with seqence")
    if (sequence[:6] == histag) or (sequence[-6:] == histag):
        ishis = True
    else:
        ishis = False

    return float(len(sequence)), prolines, ishis


#def rhToRg(rh_val, seqlen):  # Not published equation
#    """Calculates rg from rh. Takes rh_val:number, seqlen:int"""
#    k_alpha = 1.64548
#    k_nu = 0.685667
#    k_corr = 26.0435
#    k_m = 0.438091
#    Nfactored = ((seqlen+k_corr)**(k_nu))
#    rg = (k_alpha*rh_val*k_m*Nfactored)/(k_alpha*Nfactored-rh_val)
#    return rg


def rhToRg(rh_val, seqlen):  # Published equation
    """Calculates rg from rh. Takes rh_val:number, seqlen:int"""
    alpha_1 = 0.216
    alpha_2 = 4.06
    alpha_3 = 0.821
    N_33 = seqlen**0.33
    N_60 = seqlen**0.60
    rg = (rh_val*(alpha_1*alpha_2*N_33+alpha_3*(N_33-N_60)))/(alpha_1*rh_val-N_60+N_33)
    return rg


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""A collection of functions that calculates the Rh and Rg of a sequence.\nExample command:\n\n        python Marshscript.py -s GSKGYSMVTCIFPPVPGPKIKGFDAHLLEVTP\n\nThe functions are also used in the structure selection script""")
    parser.add_argument("-s", required=False, metavar="ACDEFGHIK...", help="Sequence", type=str)
    args = parser.parse_args()
    rh = rhcalculator(args.s)  
    rg = rhToRg(rh, len(args.s))
    print "Calculation of sequence:"
    print args.s
    print "Length: %s\t\t| Charge: %s" % (len(args.s), seqcharge(args.s))
    print "Rh(Marsh): %.2f\t| Rg(Converted):  %.2f" % (rh, rg)
