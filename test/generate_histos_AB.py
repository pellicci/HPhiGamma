import ROOT
import argparse
import math
import numpy as np
import sys
from functions_smuggler import Simplified_Workflow_Handler

#bools
debug          = False #Bool for verbose

#Following bools are given as input
isDataBlind    = False #Bool for blind analysis
isBDT          = False #BDT bool
isPhiAnalysis  = False # for H -> Phi Gamma
isRhoAnalysis  = False # for H -> Rho Gamma

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('Decay_channel_option', help='Type <<Phi>> for Phi, <<Rho>> for Rho') #flag for bkg estimation
p.add_argument('CR_option', help='Type <<0>> for SR, <<1>> for CR') #flag for bkg estimation
p.add_argument('isBDT_option', help='Type <<preselection>> or <<BDT>>') #flag for loose selection or tight selection (from BDT output)
p.add_argument('isBlindAnalysis', help='Type <<blind>> or <<unblind>>') #flag for loose selection or tight selection (from BDT output)
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("HPhiGammaAnalysis/mytree")
h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")

#Bools ######################################################################################################
print "################################################################################"

#SPLIT: I can split a string of chars before and after a split code (that could be a string or a symbol)
#then I can take the string standing before or after with [0] or [1], respectively. 
#[:-n] is not an emoji, it deletes last n symbols from the string
samplename =(args.rootfile_name.split("HPhiGammaAnalysis_")[1])[:-5] 
print "samplename =", samplename

if args.Decay_channel_option == "Phi":
    isPhiAnalysis = True
    print "H -> PhiGamma analysis"
else: 
    isRhoAnalysis = True
    print "H -> RhoGamma analysis"

if args.isBlindAnalysis == "blind":
    isDataBlind = True
    print "isDataBlind = ",isDataBlind

CRflag = int(args.CR_option)
if CRflag > 0 :
    print "Processing the control region ", CRflag
else :
    print "Processing the signal region" 

if args.isBDT_option == "BDT":
    isBDT = True
    BDT_OUT = 0.267508006455 #take this number running MVA/BDT_significance.py: this is the BDT output value which maximizes the significance
print "BDT = ",isBDT

print "-----------------------------------------------"

################################################################################################################
myWF = Simplified_Workflow_Handler("Signal","Data",isBDT)

#Normalization for MC dataset ################################################################################
if not samplename == "Data":
    normalization_InputFile = open("rootfiles/latest_production/MC/normalizations/Normalizations_table.txt","r")
    norm_map = dict()
    for line in normalization_InputFile:
        data_norm = line.split()
        norm_map[data_norm[0]] = float(data_norm[1])

    normalization_weight = norm_map[samplename]

    #Combine luminosity
    #luminosity2018A = 14.00 #fb^-1
    #luminosity2018B = 3.41 #fb^-1   
    #luminosity2018C = 6.94 #fb^-1
    #luminosity2018D = 31.93 #fb^-1
    luminosity = 39.54 #total lumi delivered during the trigger activity: 39.54 #fb^-1

if not samplename == "Data":
    weightSum = 0.
else: weightSum = 1.


#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_InvMass_TwoTrk_Photon","h_meson_InvMass_TwoTrk","h_firstTrk_pT","h_secondTrk_pT","h_firstTrk_Eta","h_secondTrk_Eta","h_firstTrk_Phi","h_secondTrk_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestCoupleDeltaR","h_bestJetPt","h_bestJetEta","h_firstTrk_Iso","h_firstTrk_Iso_ch","h_secondTrk_Iso","h_secondTrk_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons","h_efficiency","h_photonWP90","h_decayChannel","h_couple_Iso_neutral","h_cutOverflow","h_met_pT"]#,"h_BDT_out"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{H}",140,100.,170.) 
if   isPhiAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 1., 1.05) 
elif isRhoAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.5, 1.) 
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"p_{T} of the 1st track", 100, 20.,60.)
if   isPhiAnalysis: histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 11.,55.)
elif isRhoAnalysis: histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 5.,50.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"#eta of the 1st track", 100, -2.5,2.5)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta of the 2nd track", 100, -2.5,2.5)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#phi of the 1st track", 100, -3.14,3.14)
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi of the 2nd track", 100, -3.14,3.14)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"p_{T} of the meson", 100, 38.,110.)
histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"#eta_{meson}", 100, -2.5,2.5)
if   isPhiAnalysis: histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#Delta R_{meson}", 100, 0.,0.026)
elif isRhoAnalysis: histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#Delta R_{meson}", 100, 0.,0.07)
histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"p_{T} of the jet", 100, 40.,170.)
histo_map[list_histos[12]] = ROOT.TH1F(list_histos[12],"#eta of the jet", 100, -2.5,2.5)
histo_map[list_histos[13]] = ROOT.TH1F(list_histos[13],"Iso of the 1st track", 100, 0.45,1.)
histo_map[list_histos[14]] = ROOT.TH1F(list_histos[14],"Iso_ch of the 1st track", 100, 0.8,1.)
histo_map[list_histos[15]] = ROOT.TH1F(list_histos[15],"Iso of the 2nd track", 100, 0.3,1.)
histo_map[list_histos[16]] = ROOT.TH1F(list_histos[16],"Iso_ch of the 2nd track", 100, 0.75,1.)
histo_map[list_histos[17]] = ROOT.TH1F(list_histos[17],"Iso of the meson", 100, 0.6,1.)
histo_map[list_histos[18]] = ROOT.TH1F(list_histos[18],"Iso_ch of the meson", 100, 0.9,1.)
histo_map[list_histos[19]] = ROOT.TH1F(list_histos[19],"E_{T} of the #gamma", 100, 38.,160.)
histo_map[list_histos[20]] = ROOT.TH1F(list_histos[20],"#eta_{#gamma}", 100, -2.5,2.5)
histo_map[list_histos[21]] = ROOT.TH1F(list_histos[21],"n. of jets over pre-filters",  7, -0.5,6.5)
histo_map[list_histos[22]] = ROOT.TH1F(list_histos[22],"n. of muons", 6, -0.5,5.5)
histo_map[list_histos[23]] = ROOT.TH1F(list_histos[23],"n. of electrons", 5, -0.5,4.5)
histo_map[list_histos[24]] = ROOT.TH1F(list_histos[24],"n. of #gamma", 6, -0.5,5.5)
histo_map[list_histos[25]] = ROOT.TH1F(list_histos[25],"Efficiency steps", 7, 0.,7.)
histo_map[list_histos[26]] = ROOT.TH1F(list_histos[26],"Photon wp90 steps", 2, -0.5,1.5)
histo_map[list_histos[27]] = ROOT.TH1F(list_histos[27],"Phi or Rho channel", 2, -0.5,1.5)
histo_map[list_histos[28]] = ROOT.TH1F(list_histos[28],"Iso_neutral of the meson", 100, 0.6,1.)
histo_map[list_histos[29]] = ROOT.TH1F(list_histos[29],"Data cut overflow", 5, 0.,5.)
histo_map[list_histos[30]] = ROOT.TH1F(list_histos[30],"p_{T} of the met", 200, 0.,120.)

#histo_map[list_histos[29]] = ROOT.TH1F(list_histos[29],"BDT output", 40, -1.,1.)

#CREATE OUTPUT ROOTFILE ##################################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree #############################################################################################
mesonGammaMass       = np.zeros(1, dtype=float)
mesonMass            = np.zeros(1, dtype=float)
_coupleIsoCh         = np.zeros(1, dtype=float)
_coupleIso           = np.zeros(1, dtype=float)
_coupleIso0          = np.zeros(1, dtype=float)
_bestJetPt           = np.zeros(1, dtype=float)
_bestCouplePt        = np.zeros(1, dtype=float)
_firstTrkPt          = np.zeros(1, dtype=float)
_secondTrkPt         = np.zeros(1, dtype=float)
_photonEt            = np.zeros(1, dtype=float)
_firstTrkIso         = np.zeros(1, dtype=float)  
_secondTrkIso        = np.zeros(1, dtype=float)  
_firstTrkIsoCh       = np.zeros(1, dtype=float)  
_secondTrkIsoCh      = np.zeros(1, dtype=float) 
_firstTrkEta         = np.zeros(1, dtype=float)  
_secondTrkEta        = np.zeros(1, dtype=float)  
_bestCoupleEta       = np.zeros(1, dtype=float)  
_bestCoupleDeltaR    = np.zeros(1, dtype=float) 
_photonEta           = np.zeros(1, dtype=float)  
_eventWeight         = np.zeros(1, dtype=float)
_metPt               = np.zeros(1, dtype=float)
_nJets               = np.zeros(1, dtype=float)

tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('mesonGammaMass',mesonGammaMass,'mesonGammaMass/D')
tree_output.Branch('mesonMass',mesonMass,'mesonMass/D')
tree_output.Branch('_coupleIsoCh',_coupleIsoCh,'_coupleIsoCh/D')
tree_output.Branch('_coupleIso',_coupleIso,'_coupleIso/D')
tree_output.Branch('_coupleIso0',_coupleIso0,'_coupleIso0/D')
tree_output.Branch('_bestJetPt',_bestJetPt,'_bestJetPt/D')
tree_output.Branch('_bestCouplePt',_bestCouplePt,'_bestCouplePt/D')
tree_output.Branch('_firstTrkPt',_firstTrkPt,'_firstTrkPt/D')
tree_output.Branch('_secondTrkPt',_secondTrkPt,'_secondTrkPt/D')
tree_output.Branch('_photonEt',_photonEt,'_photonEt/D')
tree_output.Branch('_firstTrkIso',_firstTrkIso,'_firstTrkIso/D')
tree_output.Branch('_firstTrkIsoCh',_firstTrkIsoCh,'_firstTrkIsoCh/D')
tree_output.Branch('_secondTrkIso',_secondTrkIso,'_secondTrkIso/D')
tree_output.Branch('_secondTrkIsoCh',_secondTrkIsoCh,'_secondTrkIsoCh/D')
tree_output.Branch('_firstTrkEta',_firstTrkEta,'_firstTrkEta/D')
tree_output.Branch('_secondTrkEta',_secondTrkEta,'_secondTrkEta/D')
tree_output.Branch('_bestCoupleEta',_bestCoupleEta,'_bestCoupleEta/D')
tree_output.Branch('_bestCoupleDeltaR',_bestCoupleDeltaR,'_bestCoupleDeltaR/D')
tree_output.Branch('_photonEta',_photonEta,'_photonEta/D')
tree_output.Branch('_eventWeight',_eventWeight,'_eventWeight/D')
tree_output.Branch('_metPt',_metPt,'_metPt/D')
tree_output.Branch('_nJets',_nJets,'_nJets/D')
    
print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

#------------- counters -----------------
nEventsMesonAnalysis      = 0
nEventsOverLeptonVeto     = 0
nEventsAfterRegionDefiner = 0
nEventsOutOfHmassRange    = 0
nEventsOverCuts           = 0
nEventsLeftSB             = 0
nEventsRightSB            = 0


#EVENTS LOOP ##################################################################################################################### 
for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        print "nb < 0"
        continue

    
    if debug: 
        print "Processing EVENT n.",jentry+1," ..."

    #Retrieve variables from the tree 
    Hmass          = mytree.Hmass_From2K_Photon
    MesonMass      = mytree.MesonMass
    MesonIsoCh     = mytree.iso_couple_ch
    jetPt          = mytree.bestJet_pT      
    MesonPt        = mytree.bestCouplePt      
    firstTrkPt     = mytree.firstCandPt         
    secondTrkPt    = mytree.secondCandPt        
    photonEt       = mytree.photon_eT        
    firstTrketa    = mytree.firstCandEta
    secondTrketa   = mytree.secondCandEta
    firstTrkphi    = mytree.firstCandPhi
    secondTrkphi   = mytree.secondCandPhi
    jetEta         = mytree.bestJet_eta
    MesonIso       = mytree.iso_couple
    firstTrkiso    = mytree.iso_K1
    firstTrkisoCh  = mytree.iso_K1_ch
    secondTrkiso   = mytree.iso_K2
    secondTrkisoCh = mytree.iso_K2_ch
    photonEta      = mytree.photon_eta
    nJets          = mytree.nJets_25
    nMu            = mytree.nMuons
    nEle           = mytree.nElectrons
    nPhotons       = mytree.nPhotonsOverSelection
    MesonEta       = mytree.bestCoupleEta
    photonWP90     = mytree.is_photon_wp90
    isPhiEvent     = mytree.isPhi
    isRhoEvent     = mytree.isRho
    nElectrons     = mytree.nElectrons
    nMuons         = mytree.nMuons
    MesonIso0      = MesonPt/(MesonPt + (mytree.Couple_sum_pT_05 - mytree.Couple_sum_pT_05_ch))
    metPt          = mytree.met_pT
    
    #If I'm performing a PhiGamma analysis I don't want to choose those events tagged as a RhoGamma events, and viceversa
    if isPhiAnalysis and not isPhiEvent: 
        continue
    if isRhoAnalysis and not isRhoEvent: 
        continue
    nEventsMesonAnalysis+=1


    #Define Control and Signal regions: 
    if isPhiAnalysis: #for Phi meson
        if CRflag == 0 and not (MesonMass > 1.005 and MesonMass < 1.035) :
            continue

        if CRflag == 1 and (MesonMass > 1.005 and MesonMass < 1.035) :
            continue

    if isRhoAnalysis: #for Rho meson
        if CRflag == 0 and not (MesonMass > 0.63 and MesonMass < 0.91) :
            continue

        if CRflag == 1 and (MesonMass > 0.63 and MesonMass < 0.91) :
            continue

################### line used to take simmetrical sidebands ################
    if (isPhiAnalysis and MesonMass > 1.05): continue
    nEventsAfterRegionDefiner+=1
############################################################################
    
    if debug:
        print ""
        print"EVENT FEATURES"
        print"--------------------------------------"
        print "isRhoAnalysis = ",isRhoAnalysis
        print "isRhoEvent    = ",isRhoEvent
        print "CRflag        = ",CRflag
        print "isDataBlind   = ",isDataBlind
        print "MesonMass     = ",MesonMass
        print ""


    #TIGHT SELECTION from BDT output -------------------------------------------------  
    if isBDT:
        BDT_out = myWF.get_BDT_output(firstTrkisoCh,MesonIso,MesonPt,photonEt,metPt,jetPt,Hmass)
        #histo_map["h_BDT_out"].Fill(BDT_out)

        if debug: print "BDT value before selection = ", BDT_out
        if BDT_out < BDT_OUT: #Cut on BDT output
            if debug: print "BDT cut NOT passed"
            continue


    #NORMALIZATION -------------------------------------------------------------------
    #normalization for MC
    if not samplename == "Data":
        
        PUWeight               = mytree.PU_Weight
        weight_sign            = mytree.MC_Weight/abs(mytree.MC_Weight)
        photonSF, photonSF_err = myWF.get_photon_scale(mytree.photon_eT,mytree.photon_eta)
        photonSF              -= photonSF_err
        eventWeight            = weight_sign * luminosity * normalization_weight * PUWeight * photonSF
        
        if debug:
            print "EVENT WEIGHT"
            print "--------------------------------------"
            print "luminosity            = ",luminosity
            print "normalization_weight  = ",normalization_weight
            print "mytree.MC_Weight      = ",mytree.MC_Weight
            print "weight sign           = ",weight_sign
            print "PUWeight              = ",PUWeight
            print "photonSF              = ",photonSF
            print "photonSF_err          = ",photonSF_err
            print "eventWeight           = ",eventWeight
            print ""

    #normalization for DATA 
    if samplename == "Data":
        eventWeight = 1.  
        
    #--------------------------------------------------------------------------------------

    #phi angle folding
    coupleDeltaPhi = math.fabs(mytree.firstCandPhi - mytree.secondCandPhi)
    if coupleDeltaPhi > 3.14:
        coupleDeltaPhi = 6.28 - coupleDeltaPhi
    deltaR = math.sqrt((mytree.firstCandEta - mytree.secondCandEta)**2 + (coupleDeltaPhi)**2)
   

    #-------------- n events in the sidebands -----------------------------
    if (Hmass < 100. or Hmass > 170.): continue
    nEventsOutOfHmassRange+=1

    if not samplename == 'Signal':
         if (CRflag == 0 and Hmass > 100. and Hmass < 115.) : nEventsLeftSB  += 1
         if (CRflag == 0 and Hmass > 135. and Hmass < 170.) : nEventsRightSB += 1


    #Electron and muon veto ------------------------------------------
    histo_map["h_nMuons"].Fill(nMu, eventWeight) # fill this histo to check the number of electrons and muons before the veto
    histo_map["h_nElectrons"].Fill(nEle, eventWeight)
    
    if nElectrons > 0: continue
    if nMuons     > 0: continue
    nEventsOverLeptonVeto += 1


    #FILL HISTOS --------------------------------------------------------------------------
    #if DATA -> Blind Analysis on H inv mass plot
    if samplename == "Data":
        if isDataBlind:
            if Hmass < 115. or Hmass > 135.:
                histo_map["h_InvMass_TwoTrk_Photon"].Fill(Hmass, eventWeight)
        else:
            histo_map["h_InvMass_TwoTrk_Photon"].Fill(Hmass, eventWeight)
    else:
        histo_map["h_InvMass_TwoTrk_Photon"].Fill(Hmass, eventWeight)
            
    histo_map["h_meson_InvMass_TwoTrk"].Fill(MesonMass, eventWeight)       
    histo_map["h_couple_Iso_ch"].Fill(MesonIsoCh, eventWeight)
    histo_map["h_bestJetPt"].Fill(jetPt, eventWeight)
    histo_map["h_bestCouplePt"].Fill(MesonPt, eventWeight)
    histo_map["h_firstTrk_pT"].Fill(firstTrkPt, eventWeight)
    histo_map["h_secondTrk_pT"].Fill(secondTrkPt, eventWeight)
    histo_map["h_photon_energy"].Fill(photonEt, eventWeight)
    histo_map["h_photonWP90"].Fill(photonWP90, eventWeight)
    histo_map["h_firstTrk_Eta"].Fill(firstTrketa, eventWeight)    
    histo_map["h_secondTrk_Eta"].Fill(secondTrketa, eventWeight)   
    histo_map["h_firstTrk_Phi"].Fill(firstTrkphi, eventWeight)    
    histo_map["h_secondTrk_Phi"].Fill(secondTrkphi, eventWeight)   
    histo_map["h_bestCoupleDeltaR"].Fill(deltaR, eventWeight)
    histo_map["h_bestJetEta"].Fill(jetEta, eventWeight)
    histo_map["h_couple_Iso"].Fill(MesonIso, eventWeight)
    histo_map["h_firstTrk_Iso"].Fill(firstTrkiso, eventWeight)
    histo_map["h_firstTrk_Iso_ch"].Fill(firstTrkisoCh, eventWeight)
    histo_map["h_secondTrk_Iso"].Fill(secondTrkiso, eventWeight)
    histo_map["h_secondTrk_Iso_ch"].Fill(secondTrkisoCh, eventWeight)
    histo_map["h_photon_eta"].Fill(photonEta, eventWeight)
    histo_map["h_nJets_25"].Fill(nJets, eventWeight)
    histo_map["h_nPhotons"].Fill(nPhotons, eventWeight)
    histo_map["h_bestCoupleEta"].Fill(MesonEta, eventWeight)
    histo_map["h_decayChannel"].Fill(isPhiEvent, eventWeight)
    histo_map["h_couple_Iso_neutral"].Fill(MesonIso0, eventWeight)
    histo_map["h_met_pT"].Fill(metPt, eventWeight)



    #------------------------------------------------------------------------------------
    
    #FILL TREE
    mesonGammaMass[0]     = Hmass
    mesonMass[0]          = MesonMass
    _coupleIsoCh[0]       = MesonIsoCh
    _coupleIso[0]         = MesonIso
    _coupleIso0[0]        = MesonIso0
    _bestJetPt[0]         = jetPt        
    _bestCouplePt[0]      = MesonPt        
    _firstTrkPt[0]        = firstTrkPt        
    _secondTrkPt[0]       = secondTrkPt           
    _photonEt[0]          = photonEt
    _firstTrkIso[0]       = firstTrkiso        
    _secondTrkIso[0]      = secondTrkiso        
    _firstTrkIsoCh[0]     = firstTrkisoCh       
    _secondTrkIsoCh[0]    = secondTrkisoCh
    _firstTrkEta[0]       = firstTrketa        
    _secondTrkEta[0]      = secondTrketa        
    _bestCoupleEta[0]     = MesonEta
    _bestCoupleDeltaR[0]  = deltaR
    _photonEta[0]         = photonEta  
    _eventWeight[0]       = eventWeight
    _metPt[0]             = metPt
    _nJets[0]             = nJets
    
    tree_output.Fill()

    if debug:
        print "***********************************"
        print "*** EVENT RECORDED: tree filled ***"
        print "***********************************"
        print ""

    #counters
    nEventsOverCuts += 1

    if not samplename == 'Data' :
        weightSum += _eventWeight


#HISTO LABELS
histo_map["h_InvMass_TwoTrk_Photon"].GetXaxis().SetTitle("m_{trk^{+}trk^{-}#gamma} [GeV/c^2]")
histo_map["h_InvMass_TwoTrk_Photon"].SetTitle("Tracks+Photon invariant mass (Cut on phi inv. mass)")
histo_map["h_meson_InvMass_TwoTrk"].GetXaxis().SetTitle("m_{trk^{+}trk^{-}} [GeV/c^2]")
histo_map["h_meson_InvMass_TwoTrk"].SetTitle("Tracks invariant mass")
histo_map["h_firstTrk_pT"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_firstTrk_pT"].SetTitle("Transverse momentum of the first charged particle (pT_{max} of the couple)")
histo_map["h_secondTrk_pT"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_secondTrk_pT"].SetTitle("Transverse momentum of the second charged particle (pT_{min} of the couple)")
histo_map["h_firstTrk_Eta"].GetXaxis().SetTitle("#eta")
histo_map["h_firstTrk_Eta"].SetTitle("Pseudorapidity of the first charged particle (pT_{max} of the couple)")
histo_map["h_secondTrk_Eta"].GetXaxis().SetTitle("#eta")
histo_map["h_secondTrk_Eta"].SetTitle("Pseudorapidity of the second charged particle (pT_{max} of the couple)")
histo_map["h_firstTrk_Phi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_firstTrk_Phi"].SetTitle("Azimuthal angle of the first charged particle (pT_{max} of the couple)")
histo_map["h_secondTrk_Phi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_secondTrk_Phi"].SetTitle("Azimuthal angle of the second charged particle (pT_{max} of the couple)")
histo_map["h_bestCouplePt"].GetXaxis().SetTitle("pT_{trk^{+}trk^{-}} [GeV]")
histo_map["h_bestCouplePt"].SetTitle("Transverse momentum of the couple")
histo_map["h_bestCoupleEta"].GetXaxis().SetTitle("#eta")
histo_map["h_bestCoupleEta"].SetTitle("Pseudorapidity of the couple")
histo_map["h_bestJetPt"].GetXaxis().SetTitle("pT_{jet} [GeV]")
histo_map["h_bestJetPt"].SetTitle("Transverse momentum of the jet")
histo_map["h_bestJetEta"].GetXaxis().SetTitle("#eta_{jet}")
histo_map["h_bestJetEta"].SetTitle("Pseudorapidity of the jet")
histo_map["h_firstTrk_Iso"].GetXaxis().SetTitle("sum pT_{trk_{1}}/pT_{trk_{1}}")
histo_map["h_firstTrk_Iso"].SetTitle("Isolation of the trk_{1} candidate")
histo_map["h_firstTrk_Iso_ch"].GetXaxis().SetTitle("sum pT_{trk_{1}}/pT_{trk_{1}}")
histo_map["h_firstTrk_Iso_ch"].SetTitle("Isolation of the trk_{1} candidate")
histo_map["h_secondTrk_Iso"].GetXaxis().SetTitle("sum pT_{trk_{2}}/pT_{trk_{2}}")
histo_map["h_secondTrk_Iso"].SetTitle("Isolation of the trk_{2} candidate")
histo_map["h_secondTrk_Iso_ch"].GetXaxis().SetTitle("sum pT_{trk_{2}}/pT_{trk_{2}}")
histo_map["h_secondTrk_Iso_ch"].SetTitle("Isolation of the trk_{1} candidate")
histo_map["h_couple_Iso"].GetXaxis().SetTitle("sum pT_{2trk}/pT_{2trk}")
histo_map["h_couple_Iso"].SetTitle("Isolation of the couple candidate")
histo_map["h_couple_Iso_ch"].GetXaxis().SetTitle("sum pT_{2trk}/pT_{2trk}")
histo_map["h_couple_Iso_ch"].SetTitle("Isolation of the couple candidate")
histo_map["h_bestCoupleDeltaR"].GetXaxis().SetTitle("#DeltaR_{trk^{+}trk^{-}}")
histo_map["h_bestCoupleDeltaR"].SetTitle("Delta R of the couple")
histo_map["h_photon_energy"].GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
histo_map["h_photon_energy"].SetTitle("Energy of the photon")
histo_map["h_photon_eta"].GetXaxis().SetTitle("#eta_{#gamma}")
histo_map["h_photon_eta"].SetTitle("Eta of the photon")
histo_map["h_nJets_25"].GetXaxis().SetTitle("nJets over 25 GeV")
histo_map["h_nJets_25"].SetTitle("# of jets with pT > 25 GeV")
histo_map["h_nPhotons"].GetXaxis().SetTitle("n.#gamma")
histo_map["h_nPhotons"].SetTitle("n.#gamma over selections")
histo_map["h_efficiency"].GetXaxis().SetTitle("")
histo_map["h_efficiency"].GetYaxis().SetTitle("#epsilon (%)")
histo_map["h_photonWP90"].GetXaxis().SetTitle("#gamma wp90 bool")
histo_map["h_photonWP90"].GetYaxis().SetTitle("")
histo_map["h_decayChannel"].SetTitle("Decay channel")
histo_map["h_decayChannel"].GetXaxis().SetBinLabel(1,"Rho events")
histo_map["h_decayChannel"].GetXaxis().SetBinLabel(2,"Phi events")
histo_map["h_couple_Iso_neutral"].GetXaxis().SetTitle("sum pT_{2trk}/pT_{2trk}")
histo_map["h_couple_Iso_neutral"].SetTitle("Isolation of the couple candidate from neutral particles")
histo_map["h_met_pT"].GetXaxis().SetTitle("pT_{met} [GeV/c]")
histo_map["h_met_pT"].SetTitle("Transverse momentum of the missing energy")


#########################################################################
#Tree writing
tree_output.Write()
#tree_output.Scan()

#EFFICIENCY STEP PLOT MANAGEMENT
if not samplename == "Data":

    bin1content  = h_Events.GetBinContent(1)
    bin2content  = h_Events.GetBinContent(2)
    bin3content  = h_Events.GetBinContent(3)
    bin4content  = h_Events.GetBinContent(4)
    bin5content  = h_Events.GetBinContent(5)
    bin6content  = h_Events.GetBinContent(6)
    bin7content  = nEventsOverCuts
    nSignal      = bin1content
    scale_factor = 100/nSignal

    histo_map["h_efficiency"].Fill(0.5,bin1content*scale_factor)
    histo_map["h_efficiency"].Fill(1.5,bin2content*scale_factor)
    histo_map["h_efficiency"].Fill(2.5,bin3content*scale_factor)
    histo_map["h_efficiency"].Fill(3.5,bin4content*scale_factor)
    histo_map["h_efficiency"].Fill(4.5,bin5content*scale_factor)
    histo_map["h_efficiency"].Fill(5.5,bin6content*scale_factor)
    histo_map["h_efficiency"].Fill(6.5,bin7content*scale_factor)
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(1,"Events processed")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(2,"Events triggered")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(3,"Photon requested")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(4,"Best couple found")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(5,"trk-cand pT selection")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(6,"trktrk iso selection")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(7,"Tight selections")

    c11 = ROOT.TCanvas()
    c11.cd()
    histo_map["h_efficiency"].SetFillColor(1) 
    histo_map["h_efficiency"].SetFillStyle(3003)
    ROOT.gStyle.SetPaintTextFormat("4.2f %")
    ROOT.gStyle.SetOptStat(0)
    histo_map["h_efficiency"].SetMarkerSize(1.4)
    histo_map["h_efficiency"].GetXaxis().SetRangeUser(0.,7.1)
    histo_map["h_efficiency"].GetYaxis().SetRangeUser(0.,30.)
    #histo_map["h_efficiency"].SetMaximum(max(histo_map["h_efficiency"].GetHistogram().GetMaximum(),30.))
    histo_map["h_efficiency"].Draw("HIST TEXT0")
    c11.SaveAs("plots/h_efficiency.pdf")
    c11.SaveAs("plots/h_efficiency.png")    

else: #Only if data
    #Variables for cut overflow
    nEventsProcessed   = h_Events.GetBinContent(1)
    nEventsTriggered   = h_Events.GetBinContent(2)
    nEventsPhoton      = h_Events.GetBinContent(3)
    nEventsBestPair    = h_Events.GetBinContent(6)
    nEventsMesonMassSR = nEventsRightSB + nEventsLeftSB

    histo_map["h_cutOverflow"].Fill(0.5,nEventsProcessed)
    histo_map["h_cutOverflow"].Fill(1.5,nEventsTriggered)
    histo_map["h_cutOverflow"].Fill(2.5,nEventsPhoton)
    histo_map["h_cutOverflow"].Fill(3.5,nEventsBestPair)
    histo_map["h_cutOverflow"].Fill(4.5,nEventsMesonMassSR)
    histo_map["h_cutOverflow"].GetXaxis().SetBinLabel(1,"Events processed")
    histo_map["h_cutOverflow"].GetXaxis().SetBinLabel(2,"Events triggered")
    histo_map["h_cutOverflow"].GetXaxis().SetBinLabel(3,"Photon requested")
    histo_map["h_cutOverflow"].GetXaxis().SetBinLabel(4,"Best pair found")
    histo_map["h_cutOverflow"].GetXaxis().SetBinLabel(5,"Mmass SR and MGammaMass SB")

    print "CUT OVERFLOW on DATA"
    print "---------------------------------------"
    print "CRflag                    = ",CRflag," (SR if 0, CR if 1)"
    print "nEventsProcessed          = ",nEventsProcessed," (n. events in dataset)"
    print "nEventsTriggered          = ",nEventsTriggered," (n. events over HLT)"
    print "nEventsPhoton             = ",nEventsPhoton," (n. events with a best photon found)"
    print "nEventsBestPair           = ",nEventsBestPair," (n. events with a best pair found)"
    print "nEventsMesonAnalysis      = ",nEventsMesonAnalysis," (split in PhiGamma or RhoGamma analysis)"
    print "nEventsAfterRegionDefiner = ",nEventsAfterRegionDefiner," (split in SR or CR of the ditrack inv mass)"
    print "nEventsOutOfHmassRange    = ",nEventsOutOfHmassRange," (n. events in 100 < Hmass < 170 GeV)"
    print "nEventsMesonMassSR        = ",nEventsMesonMassSR," (Events in SR of MesonMass and in SBs of Hmass)"
    print "nEventsOverLeptonVeto     = ",nEventsOverLeptonVeto," (n. events without any electron or muon)"
    print "---------------------------------------"

    c11 = ROOT.TCanvas()
    c11.cd()
    histo_map["h_cutOverflow"].SetFillColor(1) 
    histo_map["h_cutOverflow"].SetFillStyle(3003)
    ROOT.gStyle.SetPaintTextFormat("4.2f")
    ROOT.gStyle.SetOptStat(0)
    histo_map["h_cutOverflow"].SetMarkerSize(1.4)
    histo_map["h_cutOverflow"].GetXaxis().SetRangeUser(1.,4.1)
#    histo_map["h_cutOverflow"].GetYaxis().SetRangeUser(0.,30.)
    #histo_map["h_cutOverflow"].SetMaximum(max(histo_map["h_cutOverflow"].GetHistogram().GetMaximum(),30.))
    histo_map["h_cutOverflow"].Draw("HIST TEXT0")
    c11.SaveAs("plots/h_cutOverflow.pdf")
    c11.SaveAs("plots/h_cutOverflow.png")

#FINAL PRINTS ###########################################################
print ""
print ""
print "########### SUMMARY ####################"
if isBDT:
    print "BDT output used = ",BDT_OUT
if samplename == "Signal":
    print "Signal MC sample (~ 2M events)"
if isBDT:
    print "n. events after preselection = ",jentry
print "n. events after cuts = " , nEventsOverCuts
if not samplename == 'Signal' and CRflag == 0:
    print "n. events in the left sideband counted = ",nEventsLeftSB
    print "n. events in the right sideband counted = ",nEventsRightSB 
    print "Total events in the sidebands = ", nEventsRightSB + nEventsLeftSB

if samplename == "Signal":
    print "Signal weight sum   = ",float(weightSum)
    print "Total signal events = ",bin1content
    print "Signal efficiency   = ",nEventsOverCuts/bin1content
print "########################################"
print ""
print ""

#HISTOS WRITING
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()
fOut.Close()
