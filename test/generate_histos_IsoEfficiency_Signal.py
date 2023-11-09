import ROOT
import argparse
import math
import numpy as np
import sys
import copy
import tdrstyle, CMS_lumi
from ROOT import gROOT
import os, sys


#bools
debug = False #Bool for verbose
isPtBin = True

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)  

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('channel', help='Type which kind of sample it is')
p.add_argument('variable', help='Type the variable')
p.add_argument('eta', help='Type the eta range')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
mytree = fInput.Get("HPhiGammaAnalysis/mytree")

CHANNEL = str(args.channel)
SAMPLE_NAME = "Signal"
VARIABLE = args.variable
ETA = args.eta

#Output
output_filename = args.outputfile_option
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Create lists of bin min e max values
pT_list = [38.,40.,42.,44.,46.,50.,55.,1000.]

# D i c t i o n a r i e s
#-------------------------------
variablesDict = {}
treeDict      = {}

nVariables = len(pT_list)
for i in range(nVariables): # Use adfor loop to create the sequential variables

    #Fill variables dict -------------------------------
    bin_fail_name = "bin{}failEvents".format(i) # Create the variable names using string formatting with the current index
    bin_pass_name = "bin{}passEvents".format(i)
    
    variablesDict[bin_fail_name] = np.zeros(1, dtype=float) # Assign each variable with a zero-initialized NumPy array
    variablesDict[bin_pass_name] = np.zeros(1, dtype=float)

    #Fill tree dict ------------------------------------
    tree_fail_name = "tree_output_bin{}_fail".format(i) # Create the tree names using string formatting with the current index
    tree_pass_name = "tree_output_bin{}_pass".format(i)
    
    treeDict[tree_fail_name] = ROOT.TTree(tree_fail_name,tree_fail_name)
    treeDict[tree_fail_name].Branch('mesonMass',variablesDict[bin_fail_name],'mesonMass/D')
    treeDict[tree_pass_name] = ROOT.TTree(tree_pass_name,tree_pass_name)
    treeDict[tree_pass_name].Branch('mesonMass',variablesDict[bin_pass_name],'mesonMass/D')

print(variablesDict)
print(treeDict)

#CMS-style plotting 
#-------------------------------
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 


#H I S T O S 
#-------------------------------
histo_map = dict()
list_histos = ["h_mesonMassPass","h_mesonMassFail"]

if CHANNEL == "Phi":
    histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{#phi} in "+SAMPLE_NAME+" if meson passes iso cut",100,1.,1.05)
    histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{#phi} in "+SAMPLE_NAME+" if meson fails iso cut",100,1.,1.05)
if CHANNEL == "Rho":
    histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{#rho} in "+SAMPLE_NAME+" if meson passes iso cut",100,0.5,1.)
    histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{#rho} in "+SAMPLE_NAME+" if meson fails iso cut",100,0.5,1.)
if CHANNEL == "K0s":
    histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{K_{0}^{*}} in "+SAMPLE_NAME+" if meson passes iso cut",100,0.8,1.)
    histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{K_{0}^{*}} in "+SAMPLE_NAME+" if meson fails iso cut",100,0.8,1.)

#histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"Trigger efficiency bin by bin in "+SAMPLE_NAME+" (TwoProngs leg)", 3, 0.,3.)


#E V E N T S   L O O P 
#---------------------------------- 

print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

#nentries = 2000000 #FIXMEEEEEEEE

#Counters
nIsolatedMesons = 0

for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry)
    if nb <= 0:
        continue

    if debug:
        print ""
        print "Processing EVENT n.",jentry+1
        print "--------------------------"

    #Retrieve variables from the tree and fill the histos
    mesonPt   = mytree.bestCouplePt
    mesonEta  = mytree.bestCoupleEta
    iso       = mytree.iso_couple
    isoCh     = mytree.iso_couple_ch
    isoNeu    = mesonPt/(mesonPt + (mytree.Couple_sum_pT_05 - mytree.Couple_sum_pT_05_ch))
    if CHANNEL == "Phi": mesonMass = mytree._PhiMass
    if CHANNEL == "Rho": mesonMass = mytree._RhoMass
    if CHANNEL == "K0s": mesonMass = mytree._K0starMass
        
    #weight
    if SAMPLE_NAME == "Signal":
        weight = mytree.PU_Weight
    else:
        weight  = 1.

    #trigger
    if not mytree.isTwoProngTrigger: continue

    #pT cut
    if mesonPt < 38.: continue

    if abs(mesonEta) > 2.1: continue

    if ETA == "barrel":
        if abs(mesonEta) > 1.444: continue

    if ETA == "endcap":
        if abs(mesonEta) < 1.444: continue

    if not VARIABLE == "IsoCh": 
        if (isoCh < 0.9): continue

    #Fill trees -----------------------------------------
    pT_bin = None

    for i in range(len(pT_list) - 1):
        if pT_list[i] <= mesonPt < pT_list[i + 1]:
            pT_bin = i

            if debug:
                print "mesonPt = ",mesonPt
                print "pT_bin n  = ",pT_bin
                print ""

            break

    # Check if the event falls within a pT bin and if probe mu passes the isolation cut
    if (pT_bin is None): continue 

    if VARIABLE == "IsoCh": 
        if(isoCh > 0.9):

            nIsolatedMesons = nIsolatedMesons + 1

            variablesDict["bin{}passEvents".format(pT_bin)][0] = mesonMass
            treeDict["tree_output_bin{}_pass".format(pT_bin)].Fill()
           
            histo_map["h_mesonMassPass"].Fill(mesonMass, weight)

        else:
            variablesDict["bin{}failEvents".format(pT_bin)][0] = mesonMass
            treeDict["tree_output_bin{}_fail".format(pT_bin)].Fill()  

            histo_map["h_mesonMassFail"].Fill(mesonMass, weight)

    if VARIABLE == "IsoNeu": 
        if(isoNeu > 0.8):

            nIsolatedMesons = nIsolatedMesons + 1

            variablesDict["bin{}passEvents".format(pT_bin)][0] = mesonMass
            treeDict["tree_output_bin{}_pass".format(pT_bin)].Fill()
           
            histo_map["h_mesonMassPass"].Fill(mesonMass, weight)

        else:
            variablesDict["bin{}failEvents".format(pT_bin)][0] = mesonMass
            treeDict["tree_output_bin{}_fail".format(pT_bin)].Fill()  

            histo_map["h_mesonMassFail"].Fill(mesonMass, weight)


# F i n a l   P r i n t s
#---------------------------------- 
print ""
print "Summary"
print "----------------"
print "n mesons passed the isolation cut = ", nIsolatedMesons
print "n mesons failed the isolation cut = ", nentries - nIsolatedMesons



# W r i t e   t r e e s
#---------------------------------- 

for i in range(nVariables):
    treeDict["tree_output_bin{}_fail".format(i)].Write()
    treeDict["tree_output_bin{}_pass".format(i)].Write()

for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()

