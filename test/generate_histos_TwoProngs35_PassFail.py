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
isPhiAnalysis = True   

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)  

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('is_Data', help='Type which kind of sample it is')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
mytree = fInput.Get("HPhiGammaTwoProngsTriggerAnalysis/mytree")

if args.is_Data == "Data":

    SAMPLE_NAME = "Data"
else:
    SAMPLE_NAME = "MC"

#Output
output_filename = args.outputfile_option
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Create lists of bin min e max values
pT_list = [38.,42.,50.,2000.]

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
    treeDict[tree_fail_name].Branch('MesonMass',variablesDict[bin_fail_name],'MesonMass/D')
    treeDict[tree_pass_name] = ROOT.TTree(tree_pass_name,tree_pass_name)
    treeDict[tree_pass_name].Branch('MesonMass',variablesDict[bin_pass_name],'MesonMass/D')

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
list_histos = ["h_mass_mumu","h_TwoProngs_pT","h_nTwoProngsPass_pT","h_nTwoProngsFail_pT","h_triggerEff_TwoProngsPt","h_MesonMass_Pass","h_MesonMass_Fail","h_MesonEta","h_triggerEff_TwoProngsEta","h_nEventsPerBin_Fail","h_nEventsPerBin_Pass","h_efficiency_perBin"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{#mu#mu} in "+SAMPLE_NAME+" (TwoProngs leg)",100,65.,115.)
histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 100, 30.,160.)
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"n. events TwoProngs35 Pass as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"n. events TwoProngs35 Fail as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"Trigger efficiency as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
if isPhiAnalysis:
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+SAMPLE_NAME+" for events Pass", 30, 1.,1.05)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"Inv mass of the pair in "+SAMPLE_NAME+" for events Fail", 30, 1.,1.05)
else:
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+SAMPLE_NAME+" for events Pass", 30, 0.5,1.)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"Inv mass of the pair in "+SAMPLE_NAME+" for events Fail", 30, 0.5,1.)

histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"Best pair #eta in "+SAMPLE_NAME+" (TwoProngs leg)", 100, 0.,2.5)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"Trigger efficiency as function of #eta^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 0.,2.5)
histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"nEvents bin by bin in "+SAMPLE_NAME+" Pass (TwoProngs leg)", 3, 0.,3.)
histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"nEvents bin by bin in "+SAMPLE_NAME+" Pass (TwoProngs leg)", 3, 0.,3.)
histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"Trigger efficiency bin by bin in "+SAMPLE_NAME+" (TwoProngs leg)", 3, 0.,3.)

#histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"Trigger efficiency bin by bin in "+SAMPLE_NAME+" (TwoProngs leg)", 3, 0.,3.)


#E V E N T S   L O O P 
#---------------------------------- 

print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

#Counters
nEventsOverCuts           = 0
nEventsTwoProngsTriggered = 0

for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        continue

    if debug:
        print ""
        print "Processing EVENT n.",jentry+1
        print "--------------------------"

    #Retrieve variables from the tree and fill the histos
    MuMuMass              = mytree.MuMuMass
    TwoProngsPt           = mytree.bestCouplePt
    isIsoMuTrigger        = mytree.isIsoMuTrigger
    isTwoProngsTrigger    = mytree.isTwoProngsTrigger      
    isBestPairFound       = mytree.isBestCoupleOfTheEvent_Found
    isPhi                 = mytree.isPhi
    isRho                 = mytree.isRho
    MesonMass             = mytree.MesonMass
    bestPairIsoCh         = mytree.iso_couple_ch
    bestPairEta           = mytree.bestCoupleEta
    
    #Isolation cut
    if bestPairIsoCh < 0.7: continue

    #channel cut
    if isPhiAnalysis:
        if isRho : 
            continue
    else:        
        if isPhi : 
            continue

    #HLT_IsoMu cut
    if not isIsoMuTrigger: continue

    nEventsOverCuts = nEventsOverCuts + 1

    #Fill trees -----------------------------------------
    pT_bin = None

    for i in range(len(pT_list) - 1):
        if pT_list[i] <= TwoProngsPt < pT_list[i + 1]:
            pT_bin = i

            if debug:
                print "TwoProngsPt = ",probeMuPt
                print "pT_bin n    = ",pT_bin
                print ""

            break

    # Check if the event falls within a pT bin and if the event passes the trigger
    if (pT_bin is None): continue 

    if isTwoProngsTrigger:

        nEventsTwoProngsTriggered = nEventsTwoProngsTriggered + 1

        variablesDict["bin{}passEvents".format(pT_bin)][0] = MesonMass
        treeDict["tree_output_bin{}_pass".format(pT_bin)].Fill()

        histo_map["h_MesonMass_Pass"].Fill(MesonMass)

    else:
        variablesDict["bin{}failEvents".format(pT_bin)][0] = MesonMass
        treeDict["tree_output_bin{}_fail".format(pT_bin)].Fill()  

        histo_map["h_MesonMass_Fail"].Fill(MesonMass)

# F i n a l   P r i n t s
#---------------------------------- 
print ""
print "Summary"
print "----------------"
print "n events fired the trigger     = ", nEventsTwoProngsTriggered
print "n probes not fired the trigger = ", nEventsOverCuts - nEventsTwoProngsTriggered



# W r i t e   t r e e s
#---------------------------------- 

for i in range(nVariables):
    treeDict["tree_output_bin{}_fail".format(i)].Write()
    treeDict["tree_output_bin{}_pass".format(i)].Write()

for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()
