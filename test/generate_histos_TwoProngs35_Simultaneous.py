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
isPhiAnalysis = True   
isPtBin = True

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

samplename =(args.rootfile_name.split("HPhiGammaTwoProngsTriggerAnalysis_output_")[1])[:-5] 
print "samplename =", samplename

#Output
output_filename = args.outputfile_option
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree 
bin0DenEvents = np.zeros(1, dtype=float)
bin0NumEvents = np.zeros(1, dtype=float)
bin1DenEvents = np.zeros(1, dtype=float)
bin1NumEvents = np.zeros(1, dtype=float)
bin2DenEvents = np.zeros(1, dtype=float)
bin2NumEvents = np.zeros(1, dtype=float)
bin3DenEvents = np.zeros(1, dtype=float)
bin3NumEvents = np.zeros(1, dtype=float)

tree_output_bin0_den = ROOT.TTree('tree_output_bin0_den','tree_output_bin0_den')
tree_output_bin0_den.Branch('MesonMass',bin0DenEvents,'MesonMass/D')
tree_output_bin0_num = ROOT.TTree('tree_output_bin0_num','tree_output_bin0_num')
tree_output_bin0_num.Branch('MesonMass',bin0NumEvents,'MesonMass/D')
tree_output_bin1_den = ROOT.TTree('tree_output_bin1_den','tree_output_bin1_den')
tree_output_bin1_den.Branch('MesonMass',bin1DenEvents,'MesonMass/D')
tree_output_bin1_num = ROOT.TTree('tree_output_bin1_num','tree_output_bin1_num')
tree_output_bin1_num.Branch('MesonMass',bin1NumEvents,'MesonMass/D')
tree_output_bin2_den = ROOT.TTree('tree_output_bin2_den','tree_output_bin0_den')
tree_output_bin2_den.Branch('MesonMass',bin2DenEvents,'MesonMass/D')
tree_output_bin2_num = ROOT.TTree('tree_output_bin2_num','tree_output_bin0_num')
tree_output_bin2_num.Branch('MesonMass',bin2NumEvents,'MesonMass/D')
tree_output_bin3_den = ROOT.TTree('tree_output_bin3_den','tree_output_bin1_den')
tree_output_bin3_den.Branch('MesonMass',bin3DenEvents,'MesonMass/D')
tree_output_bin3_num = ROOT.TTree('tree_output_bin3_num','tree_output_bin1_num')
tree_output_bin3_num.Branch('MesonMass',bin3NumEvents,'MesonMass/D')

#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

#Normalization for MC dataset ################################################################################
if SAMPLE_NAME == "MC":
    normalization_InputFile = open("rootfiles/latest_production/MC/normalizations/Normalizations_table_TwoProngs.txt","r")
    norm_map = dict()
    for line in normalization_InputFile:
        data_norm = line.split()
        norm_map[data_norm[0]] = float(data_norm[1])

    normalization_weight = norm_map[samplename]
    print "normalization_weight = ",normalization_weight

#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_mass_mumu","h_TwoProngs_pT","h_nTwoProngsTriggered_pT","h_nTwoProngsOffline_pT","h_triggerEff_TwoProngsPt","h_MesonMass","h_MesonMass_triggered","h_MesonEta","h_triggerEff_TwoProngsEta"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{#mu#mu} in "+SAMPLE_NAME+" (TwoProngs leg)",100,65.,115.)
histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 100, 30.,160.)
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"n. events TwoProngs35 triggered as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"n. events two prongs over offline selection as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"Trigger efficiency as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
if isPhiAnalysis:
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+SAMPLE_NAME+" for events at the denominator (TwoProngs leg)", 30, 1.,1.05)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"Inv mass of the pair in "+SAMPLE_NAME+" for events at the numerator (TwoProngs leg)", 30, 1.,1.05)
else:
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+SAMPLE_NAME+" for events at the denominator (TwoProngs leg)", 30, 0.5,1.)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"Inv mass of the pair in "+SAMPLE_NAME+" for events at the numerator (TwoProngs leg)", 30, 0.5,1.)

histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"Best pair #eta in "+SAMPLE_NAME+" (TwoProngs leg)", 100, 0.,2.5)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"Trigger efficiency as function of #eta^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 0.,2.5)


#    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+SAMPLE_NAME+" (TwoProngs leg)", 100, 0.5,1.)

#VARIABLES INITIALIZATION ##################################################################################################
nEventsIsoMuTrigger       = 0
nEventsTwoProngsTriggered = 0
nEventsOfflineTwoProngs   = 0
triggerFraction           = 0

#EVENTS LOOP ##################################################################################################################### 

#Create lists of bin min e max values
pT_list = [38.,45.,53.,1000.]

#Event loop 
print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        continue

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
    
    #Event weights
    weight = 1.
    if samplename == "DY10to50": weight = 0.187376434292/0.0726364661457 # DY10to50_xsec / DY50_xsec

    #Remove events not passing the offline selection
    if not isBestPairFound: continue 
    
    #Isolation cut
    if bestPairIsoCh < 0.8: continue

    if isPhiAnalysis:
        if isRho : 
            continue
    else:        
        if isPhi : 
            continue

    if debug:
        print ""
        print "Processing EVENT n.",jentry+1,"================"
        print "isIsoMuTrigger = ",isIsoMuTrigger
        print "isTwoProngsTrigger = ",isTwoProngsTrigger
        print "isBestPairFound = ",isBestPairFound
        print "Event weight = ",weight
    
    if isIsoMuTrigger:
        nEventsIsoMuTrigger       = nEventsIsoMuTrigger  + 1
    if isTwoProngsTrigger:
        nEventsTwoProngsTriggered = nEventsTwoProngsTriggered + 1
    if isBestPairFound:
        nEventsOfflineTwoProngs   = nEventsOfflineTwoProngs + 1
    
    if debug:
        print "total nEventsTwoProngsTriggered = ",nEventsTwoProngsTriggered
        print "total nEventsOfflineTwoProngs = ",nEventsOfflineTwoProngs

    #Fill histos -----------------------------------------
    histo_map["h_mass_mumu"].Fill(MuMuMass, weight)        
    histo_map["h_TwoProngs_pT"].Fill(TwoProngsPt, weight)
    
    # create pT bins for the denominator (offline selection only)
    if (TwoProngsPt >= pT_list[0] and TwoProngsPt < pT_list[1]): #bin0
        histo_map["h_MesonMass"].Fill(MesonMass, weight)
        bin0DenEvents[0] = MesonMass
        tree_output_bin0_den.Fill()
    
    if (TwoProngsPt >= pT_list[1] and TwoProngsPt < pT_list[2]): #bin1
        histo_map["h_MesonMass"].Fill(MesonMass, weight)
        bin1DenEvents[0] = MesonMass
        tree_output_bin1_den.Fill()

    if (TwoProngsPt >= pT_list[2] and TwoProngsPt < pT_list[3]): #bin2
        histo_map["h_MesonMass"].Fill(MesonMass, weight)
        bin2DenEvents[0] = MesonMass
        tree_output_bin2_den.Fill()

    #if (TwoProngsPt >= pT_list[3] and TwoProngsPt < pT_list[4]): #bin3
     #   histo_map["h_MesonMass"].Fill(MesonMass, weight)
      #  bin3DenEvents[0] = MesonMass
       # tree_output_bin3_den.Fill()

    # create pT bins for the numerator (offline selection AND HLT_TwoProngs35_Mu24)
    if isTwoProngsTrigger: 
            
            histo_map["h_nTwoProngsTriggered_pT"].Fill(TwoProngsPt, weight)
            histo_map["h_triggerEff_TwoProngsPt"].Fill(TwoProngsPt, weight) #these two histos contains the same entries. h_triggerEff_TwoProngsPt is the numerator of the trigger efficiency ratio so far.
            histo_map["h_triggerEff_TwoProngsEta"].Fill(abs(bestPairEta), weight) #these two histos contains the same entries. h_triggerEff_TwoProngsPt is the numerator of the trigger efficiency ratio so far.
            
            if (TwoProngsPt >= pT_list[0] and TwoProngsPt < pT_list[1]): #bin0
                histo_map["h_MesonMass_triggered"].Fill(MesonMass, weight)
                bin0NumEvents[0] = MesonMass
                tree_output_bin0_num.Fill()

            if (TwoProngsPt >= pT_list[1] and TwoProngsPt < pT_list[2]): #bin1
                histo_map["h_MesonMass_triggered"].Fill(MesonMass, weight)
                bin1NumEvents[0] = MesonMass
                tree_output_bin1_num.Fill()

            if (TwoProngsPt >= pT_list[2] and TwoProngsPt < pT_list[3]): #bin2
                histo_map["h_MesonMass_triggered"].Fill(MesonMass, weight)
                bin2NumEvents[0] = MesonMass
                tree_output_bin2_num.Fill()

          #  if (TwoProngsPt >= pT_list[3] and TwoProngsPt < pT_list[4]): #bin3
           #     histo_map["h_MesonMass_triggered"].Fill(MesonMass, weight)
            #    bin3NumEvents[0] = MesonMass
             #   tree_output_bin3_num.Fill()

    if isBestPairFound: 
        histo_map["h_nTwoProngsOffline_pT"].Fill(TwoProngsPt,weight)
        histo_map["h_MesonEta"].Fill(abs(bestPairEta),weight)

tree_output_bin0_den.Write()
tree_output_bin0_num.Write()
tree_output_bin1_den.Write()
tree_output_bin1_num.Write()
tree_output_bin2_den.Write()
tree_output_bin2_num.Write()
#tree_output_bin3_den.Write()
#tree_output_bin3_num.Write()

fOut.Close()

