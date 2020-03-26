import ROOT
import math
import numpy as np
import os
import sys
from decimal import Decimal

sys.setrecursionlimit(1000000)

debug = True

#INPUT SETTINGS-------------------------------------------------------------
fInput1 = ROOT.TFile("histos/latest_production/MCbackgroundHistogram.root")
fInput2 = ROOT.TFile("histos/latest_production/histos_Signal.root")

#TREE RETRIEVING------------------------------------------------------------
BKG_tree = fInput1.Get("tree_output")
Signal_tree = fInput2.Get("tree_output")

BKG_entries = BKG_tree.GetEntriesFast()
Signal_entries = Signal_tree.GetEntriesFast()
print "BKG samples have ", BKG_entries, " entries"
print "Signal sample has ", Signal_entries, " entries"



#FACTORIAL DEFINITION------------------------------------------
#it seems treefactorial takes less than math.factorial to compute
def range_prod(lo,hi):
    if lo+1 < hi:
        mid = (hi+lo)//2
        return range_prod(lo,mid) * range_prod(mid+1,hi)
    if lo == hi:
        return lo
    return lo*hi

def treefactorial(n):
    if n < 2:
        return 1
    return range_prod(1,n)

#EVENTS COUNTING------------------------------------------------
BKG_events = 0.
Signal_events = 0.
Hmass_min = 123.1
Hmass_max = 127.3

#BKG events loop
for jentry in xrange(BKG_entries):
    ientry = BKG_tree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = BKG_tree.GetEntry(jentry )
    if nb <= 0:
        continue

    if BKG_tree._HiggsMass >= Hmass_min and BKG_tree._HiggsMass <= Hmass_max:
        BKG_events += BKG_tree._nEvents #b   

#SIGNAL events loop
for jentry in xrange(Signal_entries):
    ientry = Signal_tree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = Signal_tree.GetEntry(jentry )
    if nb <= 0:
        continue

    if Signal_tree._HiggsMass >= Hmass_min and Signal_tree._HiggsMass <= Hmass_max:
        Signal_events += Signal_tree._nEvents #s   

print ""
print "BKG = ", BKG_events
print "Signal = ", Signal_events
print ""


#VARIABLES FOR UL CALCULATION-----------------------------------------------
alpha = 0.05
totalEvents = BKG_events +  Signal_events # n = s + b
s_zero = 100*Signal_events

b = int(BKG_events)
n = int(totalEvents)
summation = Decimal(0.)
delta_s = 0.01


#test values
'''
b = 130
Signal_events = 0.2
n = int(b + Signal_events)
s_zero = 100*Signal_events
'''

print "alpha = ",alpha,"  b = ",b,"  n = ",n, "  s0 starting = ",s_zero 

#RECURSIVE FUNCTION FOR UL CALCULATION--------------------------------------
def upper_limit(s):
   
   
    if s <= 0:
        if debug:
            print "NO values fallen in range!!!"     
            print "s0 = ",s + delta_s
            return s + delta_s #returns the former value of s

    summation = Decimal(0.)
#    exponential = np.exp(b+s,dtype=np.float128)
    exponential = Decimal(np.exp(b+s))

    for k in range(0,n):
           
#        powerk = np.power(b+s,k,dtype=np.float128)
        powerk = Decimal(math.pow((b+s),k))
        fac = Decimal(treefactorial(k))
        toAdd = powerk/(exponential*fac)
        summation += toAdd
        
        if debug:            
            print "--------for s0 = ",s," and k = ",k,"-----"
            print "exp = ",exponential
            print "pow = ",powerk
            print "factorial = ", fac
            print "to add = ",toAdd
            print "summation = ",summation
            print "---------------------------------------"
            
    if summation <= alpha + alpha/100 and summation >= alpha - alpha/100:  
        print ""
        print "Done!"
        print "***************************************************************"
        print "summation fallen in range = ",round(summation,5), "for s0 = ",s
        print "***************************************************************"
        return s

    else:
        return upper_limit(s-delta_s)

print  upper_limit(s_zero)    





















    


    
