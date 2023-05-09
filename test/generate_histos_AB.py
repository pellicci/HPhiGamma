import ROOT
import argparse
import math
import numpy as np
import sys
from array import array
from functions_smuggler import Simplified_Workflow_Handler

#bools
debug          = False #Bool for verbose
#isSys          = True

#Following bools are given as input
isDataBlind    = False #Bool for blind analysis
isBDT          = False #BDT bool
isPhiAnalysis  = False # for H -> Phi Gamma
isRhoAnalysis  = False # for H -> Rho Gamma
isPhotonEtaCat = False # for barrel and endcap categories

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('Decay_channel_option', help='Type <<Phi>> for Phi, <<Rho>> for Rho') #flag for bkg estimation
p.add_argument('CR_option', help='Type <<0>> for SR, <<1>> for CR') #flag for bkg estimation
p.add_argument('isBDT_option', help='Type <<preselection>> or <<BDT>>') #flag for loose selection or tight selection (from BDT output)
p.add_argument('isBlindAnalysis', help='Type <<blind>> or <<unblind>>') #flag for loose selection or tight selection (from BDT output)
p.add_argument("isEtaCat",help='Type the photon eta category',default='noPhotonEtaCat')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("HPhiGammaAnalysis/mytree")
h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")

#Trigger eff scale factors
fInput_TwoProngsTriggerSF = ROOT.TFile("scale_factors/TwoProngsTriggerSF.root") 
fInput_PhotonTriggerSF    = ROOT.TFile("scale_factors/Photon35TriggerSF.root") #Photon35 is the new histo containing SFs calculated changing the Photon eT threshold to 35 GeV at HLT level
h_TwoProngsTriggerSF = fInput_TwoProngsTriggerSF.Get("h_trigger_efficiency")
h_PhotonTriggerSF    = fInput_PhotonTriggerSF.Get("h_triggerEff_eT")

#Track iso syst scale factors
fInput_isoNeutralSF = ROOT.TFile("scale_factors/isoTrackSyst.root")
h_isoNeutralSF = fInput_isoNeutralSF.Get("h_iso_couple_all_Data")

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

if (args.isBDT_option == "BDT0" or args.isBDT_option == "BDT1"):
    isBDT = True
    fInputBDT = ROOT.TFile("MVA/BDToutput.root","READ")
    BDTtree = fInputBDT.Get("BDTtree")
    BDTtree.GetEntry(0)
    BDT_OUT = BDTtree._BDT_output
    #BDT_OUT = 0.
    print "BDT_OUT = ",BDT_OUT
    fInputBDT.Close()

if (args.isEtaCat == "EB" or args.isEtaCat == "EE"):
    isPhotonEtaCat = True
    PHOTONCAT = args.isEtaCat
else:
    PHOTONCAT = "No photon cat"


print "Category = ",args.isBDT_option, ", ",PHOTONCAT
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
list_histos = ["h_InvMass_TwoTrk_Photon","h_meson_InvMass_TwoTrk","h_firstTrk_pT","h_secondTrk_pT","h_firstTrk_Eta","h_secondTrk_Eta","h_firstTrk_Phi","h_secondTrk_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestCoupleDeltaR","h_bestJetPt","h_bestJetEta","h_firstTrk_Iso","h_firstTrk_Iso_ch","h_secondTrk_Iso","h_secondTrk_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons38WP80","h_efficiency","h_decayChannel","h_couple_Iso_neutral","h_cutOverflow","h_met_pT","h_dPhiGammaTrk","h_nPhotons20WP90","h_pTOverHmass","h_eTOverHmass","h_JetChargedEmEnergy","h_JetNeutralEmEnergy","h_JetChargedHadEnergy","h_JetNeutralHadEnergy","h_massResolution","h_genPhotonEt","h_genMesonPt","h_RecoVsGenPhotonPtRel","h_RecoVsGenMesonPtRel","h_MrecoMinusMgen","h_trkPlus_dxy","h_trkPlus_dz","h_trkPlus_dxyErr","h_trkPlus_dzErr","h_trkMinus_dxy","h_trkMinus_dz","h_trkMinus_dxyErr","h_trkMinus_dzErr","h_firstTrkCharge","h_secondTrkCharge","h_deltaRKpm"]#,"h_GenVsReco_DeltaR_trkPlus","h_GenVsReco_DeltaR_trkMinus"]#,"h_BDT_out"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{H}",300,70.,170.) 
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
histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"p_{T} of the jet", 100, 30.,170.)
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
histo_map[list_histos[24]] = ROOT.TH1F(list_histos[24],"n. of #gamma 38wp80", 6, -0.5,5.5)
histo_map[list_histos[25]] = ROOT.TH1F(list_histos[25],"Efficiency steps", 7, 0.,7.)
histo_map[list_histos[26]] = ROOT.TH1F(list_histos[26],"Phi or Rho channel", 2, -0.5,1.5)
histo_map[list_histos[27]] = ROOT.TH1F(list_histos[27],"Iso_neutral of the meson", 100, 0.6,1.)
histo_map[list_histos[28]] = ROOT.TH1F(list_histos[28],"Data cut overflow", 5, 0.,5.)
histo_map[list_histos[29]] = ROOT.TH1F(list_histos[29],"p_{T} of the met", 200, 0.,120.)
histo_map[list_histos[30]] = ROOT.TH1F(list_histos[30],"#Delta#phi between leading track and photon", 100, 0.,3.14)
histo_map[list_histos[31]] = ROOT.TH1F(list_histos[31],"n. of #gamma 20wp90", 6, -0.5,5.5)
histo_map[list_histos[32]] = ROOT.TH1F(list_histos[32],"p_{T} of the meson divided by m_{ditrk#gamma}", 100, 0.2,0.75)
histo_map[list_histos[33]] = ROOT.TH1F(list_histos[33],"E_{T} of the #gamma divided by m_{ditrk#gamma}", 100, 0.2,1.1)
histo_map[list_histos[34]] = ROOT.TH1F(list_histos[34],"Charged em energy of the jet", 100, 0.,0.001)
histo_map[list_histos[35]] = ROOT.TH1F(list_histos[35],"Neutral em energy of the jet", 100, 0.,0.7)
histo_map[list_histos[36]] = ROOT.TH1F(list_histos[36],"Charged had energy of the jet", 100, 0.,1.)
histo_map[list_histos[37]] = ROOT.TH1F(list_histos[37],"Neutral had energy of the jet", 100, 0.,0.4)
histo_map[list_histos[38]] = ROOT.TH1F(list_histos[38],"Mass resolution", 100, 0.,0.1)
histo_map[list_histos[39]] = ROOT.TH1F(list_histos[39],"genPhoton eT", 100, 38.,160.)
histo_map[list_histos[40]] = ROOT.TH1F(list_histos[40],"genMeson pT", 100, 38.,110.)
histo_map[list_histos[41]] = ROOT.TH2F(list_histos[41],"genPhotonPt vs PhotonPt/genPhotonPt;genPhotonPt;PhotonPt/genPhotonPt", 100, 38.,160.,100,0.88,1.12)
histo_map[list_histos[42]] = ROOT.TH2F(list_histos[42],"genMesonPt vs MesonPt/genMesonPt;genMesonPt;MesonPt/genMesonPt", 100, 38., 160, 100, 0.88, 1.12)
histo_map[list_histos[43]] = ROOT.TH1F(list_histos[43],"mMesonReco - mMesonGen", 100, -0.1,0.1)
histo_map[list_histos[44]] = ROOT.TH1F(list_histos[44],"trk+ dxy", 100, -0.01,0.01)
histo_map[list_histos[45]] = ROOT.TH1F(list_histos[45],"trk+ dz", 100, -0.02,0.02)
histo_map[list_histos[46]] = ROOT.TH1F(list_histos[46],"trk+ dxy error", 100, 0.,0.004)
histo_map[list_histos[47]] = ROOT.TH1F(list_histos[47],"trk+ dz error", 100, -0.,15.)#0.012)
histo_map[list_histos[48]] = ROOT.TH1F(list_histos[48],"trk- dxy", 100, -0.01,0.01)
histo_map[list_histos[49]] = ROOT.TH1F(list_histos[49],"trk- dz", 100, -0.02,0.02)
histo_map[list_histos[50]] = ROOT.TH1F(list_histos[50],"trk- dxy error", 100, 0.,0.004)
histo_map[list_histos[51]] = ROOT.TH1F(list_histos[51],"trk- dz error", 100, 0.,15.)#0.012)
histo_map[list_histos[52]] = ROOT.TH1F(list_histos[52],"trk 1 charge", 3, -1.5,1.5)
histo_map[list_histos[53]] = ROOT.TH1F(list_histos[53],"trk 2 charge", 3, -1.5,1.5)
histo_map[list_histos[54]] = ROOT.TH1F(list_histos[54],"#Delta R_{K^{#pm}}", 100, 0.,0.026)
#histo_map[list_histos[29]] = ROOT.TH1F(list_histos[29],"BDT output", 40, -1.,1.)

#CREATE OUTPUT ROOTFILE ##################################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree #############################################################################################
mesonGammaMass       = np.zeros(1, dtype=float)   
massResolution       = np.zeros(1, dtype=float)
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
_dPhiGammaTrk        = np.zeros(1, dtype=float)
_runNumber           = np.zeros(1, dtype=float)
_eventNumber         = np.zeros(1, dtype=float)
_JetChargedEmEnergy  = np.zeros(1, dtype=float)
_JetNeutralEmEnergy  = np.zeros(1, dtype=float)
_JetChargedHadEnergy = np.zeros(1, dtype=float)
_JetNeutralHadEnergy = np.zeros(1, dtype=float)

tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('mesonGammaMass',mesonGammaMass,'mesonGammaMass/D')
tree_output.Branch('massResolution',massResolution,'massResolution/D')
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
tree_output.Branch('_dPhiGammaTrk',_dPhiGammaTrk,'_dPhiGammaTrk/D')
tree_output.Branch('_runNumber',_runNumber,'_runNumber/D')
tree_output.Branch('_eventNumber',_eventNumber,'_eventNumber/D')
tree_output.Branch('_JetChargedEmEnergy',_JetChargedEmEnergy,'_JetChargedEmEnergy/D')
tree_output.Branch('_JetNeutralEmEnergy',_JetNeutralEmEnergy,'_JetNeutralEmEnergy/D')
tree_output.Branch('_JetChargedHadEnergy',_JetChargedHadEnergy,'_JetChargedHadEnergy/D')
tree_output.Branch('_JetNeutralHadEnergy',_JetNeutralHadEnergy,'_JetNeutralHadEnergy/D')

print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

#------------- counters -----------------
nEventsMesonAnalysis      = 0
nEventsOverLeptonVeto     = 0
nEventsOverDiPhotonVeto   = 0
nEventsOverVBFOrt         = 0
nEventsAfterRegionDefiner = 0
nEventsInHmassRange       = 0
nEventsOverCuts           = 0
nEventsLeftSB             = 0
nEventsRightSB            = 0
nTrkPlusMatched           = 0
nTrkMinusMatched          = 0
nPhotonMatched            = 0
nHmatched                 = 0
nMesonMatched             = 0
nHiggsFound               = 0
nMesonPtMatched           = 0
nMesonPtNotMatched        = 0
nPhotonEtMatched          = 0
nPhotonEtNotMatched       = 0
nLeadingTrkMatched        = 0
nLeadingTrkNotMatched     = 0
nSubleadingTrkMatched     = 0
nSubleadingTrkNotMatched  = 0
nOneTrkMatched            = 0

#EVENTS LOOP ##################################################################################################################### 
for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry)
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
    massResolution = mytree.photonRegressionError/photonEt #calculated as photonErr / photonEt 
    phEt_sigmaDW   = mytree.photon_eT_sigmaDW
    phEt_sigmaUP   = mytree.photon_eT_sigmaUP
    phET_scaleDW   = mytree.photon_eT_scaleDW
    phET_scaleUP   = mytree.photon_eT_scaleUP
    firstTrkCharge = mytree.firstCandCharge
    secondTrkCharge= mytree.secondCandCharge
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
    nJets          = mytree.nJets25
    nPhotons38WP80 = mytree.nPhotons38WP80
    nPhotons20WP90 = mytree.nPhotons20WP90
    MesonEta       = mytree.bestCoupleEta
    isPhiEvent     = mytree.isPhi
    isRhoEvent     = mytree.isRho
    nElectrons     = mytree.nElectrons10
    nMuons         = mytree.nMuons10
    MesonIso0      = MesonPt/(MesonPt + (mytree.Couple_sum_pT_05 - mytree.Couple_sum_pT_05_ch))
    metPt          = mytree.met_pT
    trk1dxy        = mytree.firstCand_dxy
    trk1dz         = mytree.firstCand_dz
    trk1dxyErr     = mytree.firstCand_dxyErr
    trk1dzErr      = mytree.firstCand_dzErr
    trk2dxy        = mytree.secondCand_dxy
    trk2dz         = mytree.secondCand_dz
    trk2dxyErr     = mytree.secondCand_dxyErr
    trk2dzErr      = mytree.secondCand_dzErr

    dPhiGammaTrk   = math.fabs(firstTrkphi - mytree.photon_phi)
    if dPhiGammaTrk > 3.14:
        dPhiGammaTrk = 6.28 - dPhiGammaTrk

    JetChargedEmEn = mytree.bestJet_chargedEmEnergyFraction
    JetNeutralEmEn = mytree.bestJet_neutralEmEnergyFraction
    JetChargedHadEn= (mytree.bestJet_chargedHadEnergy - MesonPt)/jetPt
    JetNeutralHadEn= mytree.bestJet_neutralHadEnergyFraction

    if not samplename == "Data" and isPhiAnalysis:
        deltaR_trkPlus = mytree.deltaR_Kplus
        deltaR_trkMinus= mytree.deltaR_Kminus
        KPlusPt      = mytree.KplusPt
        KMinusPt     = mytree.KminusPt
        KPlusEta     = mytree.Kplus_eta
        KMinusEta    = mytree.Kminus_eta
        KPlusPhi     = mytree.Kplus_phi
        KMinusPhi    = mytree.Kminus_phi

    if not samplename == "Data" and not isPhiAnalysis:
        deltaR_trkPlus = mytree.deltaR_Piplus
        deltaR_trkMinus= mytree.deltaR_Piminus

    if samplename == "Data":
        eventNumber    = mytree.event_number
        runNumber      = mytree.run_number
    else:
        isHiggsMatched  = mytree.isHiggsMatched
        isPhotonMatched = mytree.isPhotonMatched
        genPhotonEt     = mytree.genPhoton_eT
        genMesonPt      = mytree.genMeson_pT
        genMesonMass    = mytree.genMeson_m

    #If I'm performing a PhiGamma analysis I don't want to choose those events tagged as a RhoGamma events, and viceversa
    if isPhiAnalysis and not isPhiEvent: 
        continue
    if isRhoAnalysis and not isRhoEvent: 
        continue
    nEventsMesonAnalysis+=1

    #MC Truth Study -------------------------------------------
    #if not samplename == "Data":
     #   if not isPhotonMatched: continue #FIXMEEEE

    #VBF ORTHOGONALITY ------------------------------------------
    #if nJets > 1 : continue
    nEventsOverVBFOrt += 1


    #Define Control and Signal regions: ------------------------------------------
    if isPhiAnalysis: #for Phi meson
        if CRflag == 0 and not (MesonMass > 1.008 and MesonMass < 1.032) :
            continue

        if CRflag == 1 and (MesonMass > 1.008 and MesonMass < 1.032) :
            continue

        if CRflag == 3 and (MesonMass > 1.008): #left sideband tests
            continue

        if CRflag == 4 and (MesonMass < 1.032): #right sideband tests
            continue

    if isRhoAnalysis: #for Rho meson
        if CRflag == 0 and not (MesonMass > 0.62 and MesonMass < 0.92) :
            continue

        if CRflag == 1 and (MesonMass > 0.62 and MesonMass < 0.92) :
            continue

################### line used to take simmetrical sidebands ################
    if (isPhiAnalysis and MesonMass > 1.04): continue
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

    #--------------------------------------------------------------------------------------

    #phi angle folding
    coupleDeltaPhi = math.fabs(mytree.firstCandPhi - mytree.secondCandPhi)
    if coupleDeltaPhi > 3.14:
        coupleDeltaPhi = 6.28 - coupleDeltaPhi
    deltaR = math.sqrt((mytree.firstCandEta - mytree.secondCandEta)**2 + (coupleDeltaPhi)**2)

    #---------------- FIXMEEEEEE
    #if firstTrkphi > 0.65 and firstTrkphi < 0.8: continue
    #if abs(photonEta) > 2.1: continue
    #if abs(MesonEta) > 2.1: continue
    #if (Hmass > 114. and Hmass < 118.): continue #NoTail
    #if (Hmass < 114. or  Hmass > 118.): continue #Tail

    #--------------------------------------------------------------------------
    central_photonEt = photonEt
    #photonEt = phEt_sigmaUP

    #NORMALIZATION -------------------------------------------------------------------
    #normalization for MC
    if not samplename == "Data":
        
        PUWeight               = mytree.PU_Weight
        weight_sign            = mytree.MC_Weight/abs(mytree.MC_Weight) #just take the sign of the MC gen weight
        photonSF, photonSF_err = myWF.get_photon_scale(photonEt,mytree.photon_eta)
     #   photonSF              += photonSF_err
        
        #Trigger eff scale factors
        if photonEt < 80.:
            PhotonTriggerSF        = h_PhotonTriggerSF.GetBinContent(h_PhotonTriggerSF.GetXaxis().FindBin(photonEt))
            PhotonTriggerSF_err    = h_PhotonTriggerSF.GetBinError(h_PhotonTriggerSF.GetXaxis().FindBin(photonEt))
        else: PhotonTriggerSF,PhotonTriggerSF_err = h_PhotonTriggerSF.GetBinContent(9),h_PhotonTriggerSF.GetBinError(9)
        
        PhotonTriggerSF += PhotonTriggerSF_err     

        if MesonPt < 100.:
            TwoProngsTriggerSF     = h_TwoProngsTriggerSF.GetBinContent(h_TwoProngsTriggerSF.GetXaxis().FindBin(MesonPt))
            TwoProngsTriggerSF_err = h_TwoProngsTriggerSF.GetBinError(h_TwoProngsTriggerSF.GetXaxis().FindBin(MesonPt))
        else:
            TwoProngsTriggerSF,TwoProngsTriggerSF_err = h_TwoProngsTriggerSF.GetBinContent(4),h_TwoProngsTriggerSF.GetBinError(4)

#        TwoProngsTriggerSF -= TwoProngsTriggerSF_err
        

        #iso track neutral 
        isoNeutralSF     = h_isoNeutralSF.GetBinContent(h_isoNeutralSF.GetXaxis().FindBin(MesonIso0))
        isoNeutralSF_err = h_isoNeutralSF.GetBinError(h_isoNeutralSF.GetXaxis().FindBin(MesonIso0))

        #isoNeutralSF += isoNeutralSF_err

        eventWeight =  luminosity * normalization_weight * weight_sign * PUWeight * photonSF * PhotonTriggerSF * TwoProngsTriggerSF #* isoNeutralSF


        if debug:
            print "EVENT WEIGHT"
            print "--------------------------------------"
            print "phi angle trk 1        = ",firstTrkphi
            print "MesonPt                = ",MesonPt
            print "luminosity             = ",luminosity
            print "normalization_weight   = ",normalization_weight
            print "mytree.MC_Weight       = ",mytree.MC_Weight, " (this is not considered)"
            print "weight sign            = ",weight_sign
            print "PUWeight               = ",PUWeight
            print "photonSF               = ",photonSF
            print "photonSF_err           = ",photonSF_err
            print "TwoProngsTriggerSF     = ",TwoProngsTriggerSF
            print "TwoProngsTriggerSF_err = ",TwoProngsTriggerSF_err
            print "Photon Et              = ",central_photonEt
            print "Photon Et smeared down = ",phEt_sigmaDW
            print "Photon Et smeared up   = ",phEt_sigmaUP
            print "Photon Et scaled down  = ",phET_scaleDW
            print "Photon Et scaled up    = ",phET_scaleUP
            print "PhotonTriggerSF        = ",PhotonTriggerSF
            print "Meson iso0             = ",MesonIso0
            print "isoNeutralSF           = ",isoNeutralSF
    #        print "PhotonTriggerSF_err    = ",PhotonTriggerSF_err
            print "Final eventWeight **** = ",eventWeight
            print "isHiggsMatched         = ",isHiggsMatched
            print "isPhotonMatched        = ",isPhotonMatched
            print ""

    #normalization for DATA 
    if samplename == "Data":
        eventWeight = 1.  

    #Electron and muon veto ------------------------------------------
    histo_map["h_nMuons"].Fill(nMuons, eventWeight) # fill this histo to check the number of electrons and muons before the veto
    histo_map["h_nElectrons"].Fill(nElectrons, eventWeight)

    #Lepton veto
    if nElectrons > 0: continue
    if nMuons     > 0: continue
    nEventsOverLeptonVeto += 1


    #DiPhoton veto ------------------------------------------
    if nPhotons20WP90 > 1: continue
    nEventsOverDiPhotonVeto += 1

    #-------------- n events in the sidebands -----------------------------
    if (Hmass < 100. or Hmass > 170.): continue

    #TIGHT SELECTION from BDT output -------------------------------------------------  
    if isBDT:
        BDT_out = myWF.get_BDT_output(firstTrkisoCh,MesonIso0,MesonPt,photonEt,Hmass)#,JetNeutralEmEn,JetChargedHadEn,JetNeutralHadEn)
        #histo_map["h_BDT_out"].Fill(BDT_out)

        if debug: print "BDT value before selection = ", BDT_out
        if args.isBDT_option == "BDT0":
            if BDT_out < BDT_OUT: #Cut on BDT output
                if debug: print "BDT cut NOT passed"
                continue
        
        if args.isBDT_option == "BDT1":
            if (BDT_out < 0. or BDT_out > BDT_OUT):
                if debug: print "BDT cut NOT passed"
                continue

    '''
    #Photon category -------------------------------------------------
    if isPhotonEtaCat:
        if PHOTONCAT == "EB":
            if (abs(photonEta) < 0. or abs(photonEta) > 1.444): continue
        
        elif PHOTONCAT == "EE":
            if (abs(photonEta) < 1.566 or abs(photonEta) > 2.5): continue
    
    if (abs(photonEta) < 1.566 or abs(photonEta) > 2.5): continue
    if (abs(photonEta) < 0. or abs(photonEta) > 1.444): continue
    if (abs(MesonEta) < 1.566 or abs(MesonEta) > 2.5): continue
    if (abs(MesonEta) < 0. or abs(MesonEta) > 1.444): continue
    '''

    nEventsInHmassRange+=1

    if samplename == 'Data':
         if (CRflag == 0 and Hmass >  70. and Hmass < 115.) : nEventsLeftSB  += 1
         if (CRflag == 0 and Hmass > 135. and Hmass < 170.) : nEventsRightSB += 1
    
    if (abs(photonEta) > 1.444 and abs(photonEta) < 1.566): continue

    '''
    #------- MC Truth -------------------
    if not samplename == "Data":
        if deltaR_trkPlus < 0. or deltaR_trkMinus < 0.: continue
        if deltaR_trkPlus < 0.001: nTrkPlusMatched += 1
        if deltaR_trkMinus < 0.001: nTrkMinusMatched += 1
        
        if deltaR_trkPlus < 0.9*deltaR and deltaR_trkMinus < 0.9*deltaR: 
            nMesonMatched +=1
        else:
            continue

        if isPhotonMatched: nPhotonMatched += 1
        if deltaR_trkPlus < 0.002 and deltaR_trkMinus < 0.002 and isPhotonMatched: nHmatched +=1

        if isHiggsMatched and isPhotonMatched: nHiggsFound += 1
    
    '''
    #if massResolution > 0.025: continue #FIXMEEEE
    #if deltaR > 0.01: continue

    '''
    print ""
    print "genMesonPt  = ",genMesonPt
    print "MesonPt     = ",MesonPt
    print "genPhotonEt = ",genPhotonEt
    print "photonEt    = ",photonEt
    print ""
    '''
    
    if not samplename == "Data":

        #bools
        leadingTrkMatched    = True 
        subleadingTrkMatched = True 
        mesonMatched         = False
        trkPlusMatched       = True
        trkMinusMatched      = True

        #CANDIDATES CHARGE
        if firstTrkCharge > 0.: 
            trkPlusPt      = firstTrkPt
            trkPlus_dxy    = trk1dxy
            trkPlus_dz     = trk1dz
            trkPlus_dxyErr = trk1dxyErr
            trkPlus_dzErr  = trk1dzErr

            trkMinusPt      = secondTrkPt
            trkMinus_dxy    = trk2dxy
            trkMinus_dz     = trk2dxy
            trkMinus_dxyErr = trk2dxyErr
            trkMinus_dzErr  = trk2dxyErr

        else:
            trkPlusPt      = secondTrkPt
            trkPlus_dxy    = trk2dxy
            trkPlus_dz     = trk2dz
            trkPlus_dxyErr = trk2dxyErr
            trkPlus_dzErr  = trk2dzErr

            trkMinusPt      = firstTrkPt
            trkMinus_dxy    = trk1dxy
            trkMinus_dz     = trk1dxy
            trkMinus_dxyErr = trk1dxyErr
            trkMinus_dzErr  = trk1dxyErr

        if isPhiAnalysis: #to date, only PhiGamma channel has genTrk info

            #CANDIDATES SORTING
            leadingKpT    = max(KPlusPt,KMinusPt)
            subleadingKpT = min(KPlusPt,KMinusPt)


            if (firstTrkPt < 0.95*leadingKpT or firstTrkPt > 1.05*leadingKpT): leadingTrkMatched = False #range within 5 GeV
            if (secondTrkPt < 0.8*subleadingKpT or secondTrkPt > 1.2*subleadingKpT): subleadingTrkMatched = False #range within 5 GeV

            if (trkPlusPt < 0.8*KPlusPt or trkPlusPt > 1.2*KPlusPt): trkPlusMatched = False #range within 5 GeV
            if (trkMinusPt < 0.8*KMinusPt or trkMinusPt > 1.2*KMinusPt ): trkMinusMatched = False #range within 5 GeV  

            if leadingTrkMatched:
                nLeadingTrkMatched += 1
            else:
                nLeadingTrkNotMatched += 1

            if subleadingTrkMatched:
                nSubleadingTrkMatched += 1
            else:
                nSubleadingTrkNotMatched += 1

            if (leadingTrkMatched and subleadingTrkMatched): 
                nMesonPtMatched += 1 
                mesonMatched = True
            #    continue
            #else:
             #   continue

            if not leadingTrkMatched or not subleadingTrkMatched: 
                nMesonPtNotMatched += 1
                #continue
            if (leadingTrkMatched == True and subleadingTrkMatched == False) or (subleadingTrkMatched == True and leadingTrkMatched == False): 
                nOneTrkMatched +=1

            #if (leadingTrkMatched == False and subleadingTrkMatched == False): continue
            #if (leadingTrkMatched == False and subleadingTrkMatched == True): continue

            #phi angle folding
            deltaPhi_Kpm = math.fabs(KPlusPhi - KMinusPhi)
            if deltaPhi_Kpm > 3.14: deltaPhi_Kpm = 6.28 - deltaPhi_Kpm

            deltaR_Kpm = math.sqrt((KPlusEta - KMinusEta) * (KPlusEta - KMinusEta) + deltaPhi_Kpm * deltaPhi_Kpm)
        
            if (trkPlusMatched and trkMinusMatched): continue
            if (trkPlusMatched == False and trkMinusMatched == False): continue
            if (trkPlusMatched == True and trkMinusMatched == False): continue

            print"-------- Event n.",jentry+1,"----------"
            print"firstTrk charge      = ",firstTrkCharge
            print"firstTrkPt           = ",firstTrkPt
            print"leadingKpT           = ",leadingKpT
            print"leadingTrkMatched    = ",str(leadingTrkMatched)
            print ""
            print"secondTrk charge     = ",secondTrkCharge,""
            print"secondTrkPt          = ",secondTrkPt
            print"subleadingKpT        = ",subleadingKpT
            print"subleadingTrkMatched = ",str(subleadingTrkMatched)
            print ""
            print"MesonPt              = ",MesonPt
            print"genMesonPt           = ",genMesonPt
            print"MesonMatched         = ",str(mesonMatched)
            print""


        '''
        if (MesonPt < 0.95 * genMesonPt or MesonPt > 1.05 * genMesonPt): 
            nMesonPtNotMatched += 1
            continue
        else:
            nMesonPtMatched += 1
            #continue

        if (photonEt < 0.95 * genPhotonEt or photonEt > 1.05 * genPhotonEt): 
            nPhotonEtNotMatched += 1
            continue
        else:
            nPhotonEtMatched += 1
        '''

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
    histo_map["h_firstTrk_Eta"].Fill(firstTrketa, eventWeight)    
    histo_map["h_secondTrk_Eta"].Fill(secondTrketa, eventWeight)   
    histo_map["h_firstTrk_Phi"].Fill(firstTrkphi, eventWeight)    
    histo_map["h_secondTrk_Phi"].Fill(secondTrkphi, eventWeight)   
    histo_map["h_bestCoupleDeltaR"].Fill(deltaR/MesonMass, eventWeight)
    histo_map["h_bestJetEta"].Fill(jetEta, eventWeight)
    histo_map["h_couple_Iso"].Fill(MesonIso, eventWeight)
    histo_map["h_firstTrk_Iso"].Fill(firstTrkiso, eventWeight)
    histo_map["h_firstTrk_Iso_ch"].Fill(firstTrkisoCh, eventWeight)
    histo_map["h_secondTrk_Iso"].Fill(secondTrkiso, eventWeight)
    histo_map["h_secondTrk_Iso_ch"].Fill(secondTrkisoCh, eventWeight)
    histo_map["h_photon_eta"].Fill(photonEta, eventWeight)
    histo_map["h_nJets_25"].Fill(nJets, eventWeight)
    histo_map["h_nPhotons38WP80"].Fill(nPhotons38WP80, eventWeight)
    histo_map["h_nPhotons20WP90"].Fill(nPhotons20WP90, eventWeight)
    histo_map["h_bestCoupleEta"].Fill(MesonEta, eventWeight)
    histo_map["h_decayChannel"].Fill(isPhiEvent, eventWeight)
    histo_map["h_couple_Iso_neutral"].Fill(MesonIso0, eventWeight)
    histo_map["h_met_pT"].Fill(metPt, eventWeight)
    histo_map["h_dPhiGammaTrk"].Fill(dPhiGammaTrk,eventWeight)
    histo_map["h_pTOverHmass"].Fill(MesonPt/Hmass, eventWeight)
    histo_map["h_eTOverHmass"].Fill(photonEt/Hmass, eventWeight)
    histo_map["h_JetChargedEmEnergy"].Fill(JetChargedEmEn, eventWeight)
    histo_map["h_JetNeutralEmEnergy"].Fill(JetNeutralEmEn, eventWeight)
    histo_map["h_JetChargedHadEnergy"].Fill(JetChargedHadEn, eventWeight)
    histo_map["h_JetNeutralHadEnergy"].Fill(JetNeutralHadEn, eventWeight)
    histo_map["h_JetNeutralHadEnergy"].Fill(JetNeutralHadEn, eventWeight)
    histo_map["h_massResolution"].Fill(massResolution, eventWeight)
    histo_map["h_trkPlus_dxy"].Fill(trkPlus_dxy,eventWeight)
    histo_map["h_trkPlus_dz"].Fill(trkPlus_dz,eventWeight)
    histo_map["h_trkPlus_dxyErr"].Fill(trkPlus_dxyErr,eventWeight)
    histo_map["h_trkPlus_dzErr"].Fill(trkPlus_dzErr/trkPlus_dz,eventWeight)
    histo_map["h_trkMinus_dxy"].Fill(trkMinus_dxy,eventWeight)
    histo_map["h_trkMinus_dz"].Fill(trkMinus_dz,eventWeight)
    histo_map["h_trkMinus_dxyErr"].Fill(trkMinus_dxyErr,eventWeight)
    histo_map["h_trkMinus_dzErr"].Fill(trkMinus_dzErr/trkMinus_dz,eventWeight)
    histo_map["h_firstTrkCharge"].Fill(firstTrkCharge,eventWeight)
    histo_map["h_secondTrkCharge"].Fill(secondTrkCharge,eventWeight)

    if not samplename == "Data":
        histo_map["h_genPhotonEt"].Fill(genPhotonEt, eventWeight)
        histo_map["h_genMesonPt"].Fill(genMesonPt, eventWeight)
        histo_map["h_RecoVsGenPhotonPtRel"].Fill(genPhotonEt, (photonEt/genPhotonEt)) #TH2F
        histo_map["h_RecoVsGenMesonPtRel"].Fill(genMesonPt, (MesonPt/genMesonPt)) #TH2F
        histo_map["h_MrecoMinusMgen"].Fill(MesonMass - genMesonMass, eventWeight)
        histo_map["h_deltaRKpm"].Fill(deltaR_Kpm, eventWeight)

    #------------------------------------------------------------------------------------
    #If preselection, reweight events as function of the mass resolution. Just for the BDT.
    if not (isBDT and samplename == "Data"): eventWeight = eventWeight/massResolution

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
    _dPhiGammaTrk[0]      = dPhiGammaTrk
    _JetChargedEmEnergy[0]  = JetChargedEmEn
    _JetNeutralEmEnergy[0]  = JetNeutralEmEn
    _JetChargedHadEnergy[0] = JetChargedHadEn
    _JetNeutralHadEnergy[0] = JetNeutralHadEn
    if samplename == "Data":
        _eventNumber[0]       = eventNumber
        _runNumber[0]         = runNumber
    
    if not isBDT: #this is done to train the BDT with events under the peak
        if (Hmass > 115. and Hmass < 135.): tree_output.Fill()
    else:
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
histo_map["h_nPhotons38WP80"].GetXaxis().SetTitle("n.#gamma_{38wp80}")
histo_map["h_nPhotons38WP80"].SetTitle("n.#gamma over selections")
histo_map["h_nPhotons20WP90"].GetXaxis().SetTitle("n.#gamma_{20wp90}")
histo_map["h_nPhotons20WP90"].SetTitle("n.#gamma over selections")
histo_map["h_efficiency"].GetXaxis().SetTitle("")
histo_map["h_efficiency"].GetYaxis().SetTitle("#epsilon (%)")
histo_map["h_decayChannel"].SetTitle("Decay channel")
histo_map["h_decayChannel"].GetXaxis().SetBinLabel(1,"Rho events")
histo_map["h_decayChannel"].GetXaxis().SetBinLabel(2,"Phi events")
histo_map["h_couple_Iso_neutral"].GetXaxis().SetTitle("sum pT_{2trk}/pT_{2trk}")
histo_map["h_couple_Iso_neutral"].SetTitle("Isolation of the couple candidate from neutral particles")
histo_map["h_met_pT"].GetXaxis().SetTitle("pT_{met} [GeV/c]")
histo_map["h_met_pT"].SetTitle("Transverse momentum of the missing energy")
histo_map["h_dPhiGammaTrk"].GetXaxis().SetTitle("#Delta#Phi_{#gamma,Trk_1} [rad]")
histo_map["h_dPhiGammaTrk"].SetTitle("#Delta#Phi_{#gamma,Trk_1} [rad]")
histo_map["h_pTOverHmass"].GetXaxis().SetTitle("p_T^{trk^{+}trk^{-}}/m_{ditrk#gamma}")
histo_map["h_pTOverHmass"].SetTitle("pT of the meson divided by Hmass")
histo_map["h_eTOverHmass"].GetXaxis().SetTitle("eT_{#gamma}/m_{ditrk#gamma}")
histo_map["h_eTOverHmass"].SetTitle("eT of the photon divided by Hmass")
histo_map["h_JetChargedEmEnergy"].GetXaxis().SetTitle("E_{jet}^{ch,em}")
histo_map["h_JetNeutralEmEnergy"].GetXaxis().SetTitle("E_{jet}^{0,em}")
histo_map["h_JetChargedHadEnergy"].GetXaxis().SetTitle("E_{jet}^{ch,had}")
histo_map["h_JetNeutralHadEnergy"].GetXaxis().SetTitle("E_{jet}^{0,had}")
histo_map["h_massResolution"].GetXaxis().SetTitle("#sigmaM/M")
histo_map["h_trkPlus_dxy"].GetXaxis().SetTitle("dxy_{trk^{+}}")
histo_map["h_trkPlus_dz"].GetXaxis().SetTitle("dz_{trk^{+}}")
histo_map["h_trkPlus_dxyErr"].GetXaxis().SetTitle("#sigma_dxy^{trk^{+}}")
histo_map["h_trkPlus_dzErr"].GetXaxis().SetTitle("#sigma_dz^{trk^{+}}")
histo_map["h_trkMinus_dxy"].GetXaxis().SetTitle("dxy_{trk^{-}}")
histo_map["h_trkMinus_dz"].GetXaxis().SetTitle("dz_{trk^{-}}")
histo_map["h_trkMinus_dxyErr"].GetXaxis().SetTitle("#sigma_dxy^{trk^{-}}")
histo_map["h_trkMinus_dzErr"].GetXaxis().SetTitle("#sigma_dz^{trk^{-}}")
histo_map["h_firstTrkCharge"].GetXaxis().SetTitle("q_{trk_{1}}")
histo_map["h_secondTrkCharge"].GetXaxis().SetTitle("q_{trk_{2}}")

if not samplename == "Data":
    histo_map["h_genPhotonEt"].GetXaxis().SetTitle("E_{T}^{gen#gamma}[GeV]")
    histo_map["h_genPhotonEt"].SetTitle("eT of the genPhoton")
    histo_map["h_genMesonPt"].GetXaxis().SetTitle("pT_{genMeson} [GeV]")
    histo_map["h_genMesonPt"].SetTitle("Transverse momentum of the genMeson")
    histo_map["h_RecoVsGenPhotonPtRel"].GetXaxis().SetTitle("E_{T}^{gen#gamma}[GeV]")
    histo_map["h_RecoVsGenPhotonPtRel"].GetYaxis().SetTitle("E_{T}^{reco#gamma}/E_{T}^{gen#gamma}")
    histo_map["h_RecoVsGenMesonPtRel"].GetXaxis().SetTitle("p_{T}^{genMeson}[GeV]")
    histo_map["h_RecoVsGenMesonPtRel"].GetYaxis().SetTitle("p_{T}^{recoMeson}/p_{T}^{genMeson}")
    histo_map["h_MrecoMinusMgen"].GetXaxis().SetTitle("m_{recoMeson}-m_{genMeson} [GeV]")
    histo_map["h_deltaRKpm"].GetXaxis().SetTitle("#DeltaR_{K^{#pm}}")

 #   histo_map["h_GenVsReco_DeltaR_trkPlus"].GetXaxis().SetTitle("#DeltaR_{GenVsReco}^{trk^{+}}")
  #  histo_map["h_GenVsReco_DeltaR_trkPlus"].SetTitle("#DeltaR_{GenVsReco}^{trk^{+}}")
  #  histo_map["h_GenVsReco_DeltaR_trkMinus"].GetXaxis().SetTitle("#DeltaR_{GenVsReco}^{trk^{-}}")
  #  histo_map["h_GenVsReco_DeltaR_trkMinus"].SetTitle("#DeltaR_{GenVsReco}^{trk^{-}}")

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
    bin7content  = h_Events.GetBinContent(7)
    bin8content  = nEventsOverCuts
    nSignal      = bin1content
    scale_factor = 100/nSignal

    histo_map["h_efficiency"].Fill(0.5,bin1content*scale_factor)
    histo_map["h_efficiency"].Fill(1.5,bin2content*scale_factor)
    histo_map["h_efficiency"].Fill(2.5,bin3content*scale_factor)
    histo_map["h_efficiency"].Fill(3.5,bin4content*scale_factor)
    histo_map["h_efficiency"].Fill(4.5,bin5content*scale_factor)
    histo_map["h_efficiency"].Fill(5.5,bin6content*scale_factor)
    histo_map["h_efficiency"].Fill(6.5,bin7content*scale_factor)
    histo_map["h_efficiency"].Fill(7.5,bin8content*scale_factor)
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(1,"Events processed")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(2,"Events triggered")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(3,"Photon requested")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(4,"Best couple found")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(5,"trk-cand pT selection")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(6,"trktrk iso selection")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(7,"VBF veto")
    histo_map["h_efficiency"].GetXaxis().SetBinLabel(8,"Tight selections")

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
    print "nEventsOverVBFOrt         = ",nEventsOverVBFOrt," (n. events over the VBF orthogonality: events with more than one jet are discarded)"
    print "nEventsAfterRegionDefiner = ",nEventsAfterRegionDefiner," (split in SR or CR of the ditrack inv mass)"
    print "nEventsOverLeptonVeto     = ",nEventsOverLeptonVeto," (n. events without any electron or muon)"
    print "nEventsOverDiPhotonVeto   = ",nEventsOverDiPhotonVeto," (n. events with just one photon)"
    print "nEventsInHmassRange       = ",nEventsInHmassRange," (n. events in 100 < Hmass < 170 GeV)"
    print "nEventsMesonMassSR        = ",nEventsMesonMassSR," (Events in SR of MesonMass and in SBs of Hmass)"
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
print "----------- SUMMARY -----------------------"
if isBDT:
    print "BDT output used = ",BDT_OUT
if not samplename == "Data":
    print "Signal MC sample (~ 2M events)"
if isBDT:
    print "n. events after preselection = ",jentry
print "n. events after cuts = " , nEventsOverCuts

if samplename == 'Data' and CRflag == 0:
    print "n. events in the left sideband counted = ",nEventsLeftSB
    print "n. events in the right sideband counted = ",nEventsRightSB 
    print "Total events in the sidebands = ", nEventsRightSB + nEventsLeftSB

if not samplename == "Data":
    print "n Photon matched = ",nPhotonMatched
    print "n trk +  matched = ",nTrkPlusMatched
    print "n trk -  matched = ",nTrkMinusMatched
    print "n meson  matched = ",nMesonMatched
    print "n Higgs  matched = ",nHmatched
    print "Signal weight sum   = ",float(weightSum)
    print "Signal integral     = ",histo_map["h_InvMass_TwoTrk_Photon"].Integral()
    print "Total signal events = ",bin1content
    print "Signal efficiency   = ",nEventsOverCuts/bin1content
    print "nMesonPtMatched           = ",nMesonPtMatched
    print "nMesonPtNotMatched        = ",nMesonPtNotMatched," (", 100*nMesonPtNotMatched/(nMesonPtMatched + nMesonPtNotMatched)," percento of total)"
    print "nLeadingTrkMatched        = ",nLeadingTrkMatched
    print "nLeadingTrkNotMatched     = ",nLeadingTrkNotMatched
    print "nSubleadingTrkMatched     = ",nSubleadingTrkMatched
    print "nSubleadingTrkNotMatched  = ",nSubleadingTrkNotMatched
    print "nOneTrkMatched            = ",nOneTrkMatched

    #print "nPhotonEtMatched    = ",nPhotonEtMatched
    #print "nPhotonEtNotMatched = ",nPhotonEtNotMatched," (", 100*nPhotonEtNotMatched/(nPhotonEtMatched + nPhotonEtNotMatched)," percento of total)"
    
print "-------------------------------------------"
print ""
print ""

#HISTOS WRITING
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()
fOut.Close()
