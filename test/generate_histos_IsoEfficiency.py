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
p.add_argument('is_Data', help='Type which kind of sample it is')
p.add_argument('variable', help='Type the variable')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
mytree = fInput.Get("HPhiGammaIsoEfficiencyAnalysis/mytree")

if args.is_Data == "Data":
    SAMPLE_NAME = "Data"
else:
    SAMPLE_NAME = "MC"

VARIABLE = args.variable

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
    treeDict[tree_fail_name].Branch('mumuMass',variablesDict[bin_fail_name],'mumuMass/D')
    treeDict[tree_pass_name] = ROOT.TTree(tree_pass_name,tree_pass_name)
    treeDict[tree_pass_name].Branch('mumuMass',variablesDict[bin_pass_name],'mumuMass/D')

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
list_histos = ["h_mumuMassPass","h_mumuMassFail","h_tagMuPt","h_probeMuPt","h_isoCh","h_isoNeu","h_tagMuPhi","h_probeMuPhi","h_tagMuEta","h_probeMuEta"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{#mu#mu} in "+SAMPLE_NAME+" if probe mu passes iso cut",100,60.,120.)
histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{#mu#mu} in "+SAMPLE_NAME+" if probe mu fails iso cut",100,60.,120.)
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"p_{T}^{tag#mu} in "+SAMPLE_NAME, 100, 10.,160.)
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T}^{probe#mu} in "+SAMPLE_NAME, 100, 10.,160.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"iso^{probe#mu}_{ch} in "+SAMPLE_NAME, 100, 0.9,1.)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"iso^{probe#mu}_{neu} in "+SAMPLE_NAME, 100, 0.5,1.)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#phi^{tag#mu} in "+SAMPLE_NAME, 100, -3.14,3.14)
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi^{probe#mu} in "+SAMPLE_NAME, 100, -3.14,3.14)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"#eta^{tag#mu} in "+SAMPLE_NAME, 100, -2.5,2.5)
histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"#eta^{probe#mu} in "+SAMPLE_NAME, 100, -2.5,2.5)

#histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"Trigger efficiency bin by bin in "+SAMPLE_NAME+" (TwoProngs leg)", 3, 0.,3.)


#E V E N T S   L O O P 
#---------------------------------- 

print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

nentries = 500000 #FIXMEEEEEEEE

#Counters
nIsolatedProbes = 0

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
    tagMuPhi   = mytree.tagMuPhi
    tagMuEta   = mytree.tagMuEta
    tagMuPt    = mytree.tagMuPt
    probeMuPhi = mytree.probeMuPhi
    probeMuEta = mytree.probeMuEta
    probeMuPt  = mytree.probeMuPt
    iso        = mytree.iso
    isoCh      = mytree.iso_ch
    isoNeu     = probeMuPt/(probeMuPt + (mytree.sum_pT_05 - mytree.sum_pT_05_ch))
    mumuMass   = mytree.MuMuMass
    
    #weight
    if SAMPLE_NAME == "MC":
        weight = mytree.PU_Weight
    else:
        weight  = 1.

    #pT cut
    if probeMuPt < 38.: continue

    if not VARIABLE == "IsoCh": 
        if (isoCh < 0.9): continue

    #Fill trees -----------------------------------------
    pT_bin = None

    for i in range(len(pT_list) - 1):
        if pT_list[i] <= probeMuPt < pT_list[i + 1]:
            pT_bin = i

            if debug:
                print "probeMuPt = ",probeMuPt
                print "pT_bin n  = ",pT_bin
                print ""

            break

    # Check if the event falls within a pT bin and if probe mu passes the isolation cut
    if (pT_bin is None): continue 

    if VARIABLE == "IsoCh": 
        if(isoCh > 0.9):

            nIsolatedProbes = nIsolatedProbes + 1

            variablesDict["bin{}passEvents".format(pT_bin)][0] = mumuMass
            treeDict["tree_output_bin{}_pass".format(pT_bin)].Fill()
           
            histo_map["h_mumuMassPass"].Fill(mumuMass, weight)
            histo_map["h_tagMuPt"].Fill(tagMuPt, weight)
            histo_map["h_probeMuPt"].Fill(probeMuPt, weight)
            histo_map["h_tagMuPhi"].Fill(tagMuPhi, weight)
            histo_map["h_probeMuPhi"].Fill(probeMuPhi, weight)
            histo_map["h_tagMuEta"].Fill(tagMuEta, weight)
            histo_map["h_probeMuEta"].Fill(probeMuEta, weight)
            histo_map["h_isoCh"].Fill(isoCh, weight)
            histo_map["h_isoNeu"].Fill(isoNeu, weight)

        else:
            variablesDict["bin{}failEvents".format(pT_bin)][0] = mumuMass
            treeDict["tree_output_bin{}_fail".format(pT_bin)].Fill()  

            histo_map["h_mumuMassFail"].Fill(mumuMass, weight)

    if VARIABLE == "IsoNeu": 
        if(isoNeu > 0.7):

            nIsolatedProbes = nIsolatedProbes + 1

            variablesDict["bin{}passEvents".format(pT_bin)][0] = mumuMass
            treeDict["tree_output_bin{}_pass".format(pT_bin)].Fill()
           
            histo_map["h_mumuMassPass"].Fill(mumuMass, weight)
            histo_map["h_tagMuPt"].Fill(tagMuPt, weight)
            histo_map["h_probeMuPt"].Fill(probeMuPt, weight)
            histo_map["h_tagMuPhi"].Fill(tagMuPhi, weight)
            histo_map["h_probeMuPhi"].Fill(probeMuPhi, weight)
            histo_map["h_tagMuEta"].Fill(tagMuEta, weight)
            histo_map["h_probeMuEta"].Fill(probeMuEta, weight)
            histo_map["h_isoCh"].Fill(isoCh, weight)
            histo_map["h_isoNeu"].Fill(isoNeu, weight)

        else:
            variablesDict["bin{}failEvents".format(pT_bin)][0] = mumuMass
            treeDict["tree_output_bin{}_fail".format(pT_bin)].Fill()  

            histo_map["h_mumuMassFail"].Fill(mumuMass, weight)

# F i n a l   P r i n t s
#---------------------------------- 
print ""
print "Summary"
print "----------------"
print "n probes passed the isolation cut = ", nIsolatedProbes
print "n probes failed the isolation cut = ", nentries - nIsolatedProbes



# W r i t e   t r e e s
#---------------------------------- 

for i in range(nVariables):
    treeDict["tree_output_bin{}_fail".format(i)].Write()
    treeDict["tree_output_bin{}_pass".format(i)].Write()

for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()

