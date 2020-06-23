import ROOT
import math
import numpy as np
import argparse
import os

debug = True
isDataBlind = False

p = argparse.ArgumentParser(description='Select directory')
p.add_argument('dir_path', help='Type directory path')
args = p.parse_args()

#INPUT SETTINGS-------------------------------------------------------------
fInput1 = ROOT.TFile(args.dir_path + "histos_Data.root")
fInput2 = ROOT.TFile(args.dir_path + "histos_Signal.root")

#OUTPUT SETTINGS-----------------------------------------------------------
fOutput = ROOT.TFile("nEvents_Container.root","RECREATE")
fOutput.cd()
nEvents_BKG = np.zeros(1, dtype=float)
nEvents_Signal = np.zeros(1, dtype=float)
tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('nEvents_BKG',nEvents_BKG,'nEvents_BKG/D')
tree_output.Branch('nEvents_Signal',nEvents_Signal,'nEvents_Signal/D')

#TREE RETRIEVING------------------------------------------------------------
BKG_tree = fInput1.Get("tree_output")
Signal_tree = fInput2.Get("tree_output")

BKG_entries = BKG_tree.GetEntriesFast()
Signal_entries = Signal_tree.GetEntriesFast()

#EVENTS COUNTING------------------------------------------------
BKG_events = 0.
Signal_events = 0.
KKgMass_min = 120.
KKgMass_max = 130.

#BKG events loop
for jentry in xrange(BKG_entries):
    ientry = BKG_tree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = BKG_tree.GetEntry(jentry )
    if nb <= 0:
        continue

    if isDataBlind:
        if BKG_tree.mass_KKg <= KKgMass_min or BKG_tree.mass_KKg >= KKgMass_max:
            BKG_events += BKG_tree._eventWeight #b   
    else:
        BKG_events += BKG_tree._eventWeight #b   

#SIGNAL events loop
for jentry in xrange(Signal_entries):
    ientry = Signal_tree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = Signal_tree.GetEntry(jentry )
    if nb <= 0:
        continue

    if Signal_tree.mass_KKg >= KKgMass_min and Signal_tree.mass_KKg <= KKgMass_max:
        Signal_events += Signal_tree._eventWeight #s   

nEvents_BKG[0] = BKG_events
nEvents_Signal[0] = Signal_events
tree_output.Fill()
tree_output.Write()
tree_output.Scan()
fOutput.Close()

print "-----------------------------------------------------------"
print "BKG samples have ", BKG_entries, " entries (bkg taken from data)" 
print "Signal sample has ", Signal_entries, " entries"
print "-----------------------------------------------------------"
print "Events in sidebands:"
print "nBKG = ", BKG_events
print "nSignal = ", Signal_events
print "***********************************************************"
print "SIGNIFICANCE: Signal/sqrt(BKG) =",Signal_events/(BKG_events)**(0.5)
print "***********************************************************"
