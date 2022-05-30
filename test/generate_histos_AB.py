import ROOT
import argparse
import math
import numpy as np
import sys
from functions_smuggler import Simplified_Workflow_Handler

#bools
debug          = False #Bool for verbose
doSelection    = False #if you turn this to False disable select_all_but_one

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

#Bools ######################################################################################################
print "########################################"

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
    BDT_OUT = 0.0775711324343 #take this number running MVA/BDT_significance.py: this is the BDT output value which maximizes the significance
print "BDT = ",isBDT


#SPLIT: I can split a string of chars before and after a split code (that could be a string or a symbol)
#then I can take the string standing before or after with [0] or [1], respectively. 
#[:-n] is not an emoji, it deletes last n symbols from the string
samplename =(args.rootfile_name.split("HPhiGammaAnalysis_")[1])[:-5] 
print "samplename =", samplename

print "######################################"

################################################################################################################
myWF = Simplified_Workflow_Handler("Signal","Data",isBDT)

#Normalization for MC dataset ################################################################################
if not samplename == "Data":
    normalization_InputFile = open("rootfiles/latest_production/MC/normalizations/Normalizations_table.txt","r")
    h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")
    norm_map = dict()
    for line in normalization_InputFile:
        data_norm = line.split()
        norm_map[data_norm[0]] = float(data_norm[1])

    normalization_weight = norm_map[samplename]

    #Combine luminosity
    luminosity2018A = 14.00 #fb^-1
    luminosity2018B = 3.41 #fb^-1   
    luminosity2018C = 6.94 #fb^-1
    luminosity2018D = 31.93 #fb^-1
    luminosity      = 39.54 #total lumi delivered during the trigger activity: 39.54 #fb^-1

if samplename == "Signal":
    sigWeightSum = 0.

#MCfromDATA correction for photons ###############################################################################
ph_ID_scale_name_2018  = "scale_factors/2018_PhotonsMVAwp90.root"
ph_ID_scale_file_2018  = ROOT.TFile(ph_ID_scale_name_2018)
ph_ID_scale_histo_2018 = ROOT.TH2F()
ph_ID_scale_histo_2018 = ph_ID_scale_file_2018.Get("EGamma_SF2D")

ph_pixVeto_scale_name_2018  = "scale_factors/HasPix_2018.root"
ph_pixVeto_scale_file_2018  = ROOT.TFile(ph_pixVeto_scale_name_2018)
ph_pixVeto_scale_histo_2018 = ROOT.TH1F()
ph_pixVeto_scale_histo_2018 = ph_pixVeto_scale_file_2018.Get("eleVeto_SF")

#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_InvMass_TwoTrk_Photon","h_InvMass_TwoTrk_Photon_NoPhiMassCut","h_meson_InvMass_TwoTrk","h_firstTrk_pT","h_secondTrk_pT","h_firstTrk_Eta","h_secondTrk_Eta","h_firstTrk_Phi","h_secondTrk_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestCoupleDeltaR","h_bestJetPt","h_bestJetEta","h_firstTrk_Iso","h_firstTrk_Iso_ch","h_secondTrk_Iso","h_secondTrk_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons","h_efficiency","h_photonWP90","h_decayChannel"]#,"h_BDT_out"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{H}",100,100.,150.) 
histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{H} (no cut on phi mass)",100,100.,150.) 
if   isPhiAnalysis: histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"M_{meson}", 100, 1., 1.04) 
elif isRhoAnalysis: histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"M_{meson}", 100, 0.5, 1.) 
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 1st track", 100, 20.,60.)
if   isPhiAnalysis: histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"p_{T} of the 2nd track", 100, 11.,55.)
elif isRhoAnalysis: histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"p_{T} of the 2nd track", 100, 5.,50.)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta of the 1st track", 100, -2.5,2.5)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#eta of the 2nd track", 100, -2.5,2.5)
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi of the 1st track", 100, -3.14,3.14)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"#phi of the 2nd track", 100, -3.14,3.14)
histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"p_{T} of the meson", 100, 38.,110.)
histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#eta_{meson}", 100, -2.5,2.5)
if   isPhiAnalysis: histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"#Delta R_{meson}", 100, 0.,0.026)
elif isRhoAnalysis: histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"#Delta R_{meson}", 100, 0.,0.07)
histo_map[list_histos[12]] = ROOT.TH1F(list_histos[12],"p_{T} of the jet", 100, 40.,170.)
histo_map[list_histos[13]] = ROOT.TH1F(list_histos[13],"#eta of the jet", 100, -2.5,2.5)
histo_map[list_histos[14]] = ROOT.TH1F(list_histos[14],"Iso of the 1st track", 100, 0.45,1.)
histo_map[list_histos[15]] = ROOT.TH1F(list_histos[15],"Iso_ch of the 1st track", 100, 0.8,1.)
histo_map[list_histos[16]] = ROOT.TH1F(list_histos[16],"Iso of the 2nd track", 100, 0.3,1.)
histo_map[list_histos[17]] = ROOT.TH1F(list_histos[17],"Iso_ch of the 2nd track", 100, 0.75,1.)
histo_map[list_histos[18]] = ROOT.TH1F(list_histos[18],"Iso of the meson", 100, 0.6,1.)
histo_map[list_histos[19]] = ROOT.TH1F(list_histos[19],"Iso_ch of the meson", 100, 0.9,1.)
histo_map[list_histos[20]] = ROOT.TH1F(list_histos[20],"E_{T} of the #gamma", 100, 38.,160.)
histo_map[list_histos[21]] = ROOT.TH1F(list_histos[21],"#eta_{#gamma}", 100, -2.5,2.5)
histo_map[list_histos[22]] = ROOT.TH1F(list_histos[22],"n. of jets over pre-filters",  7, -0.5,6.5)
histo_map[list_histos[23]] = ROOT.TH1F(list_histos[23],"n. of muons", 6, -0.5,5.5)
histo_map[list_histos[24]] = ROOT.TH1F(list_histos[24],"n. of electrons", 5, -0.5,4.5)
histo_map[list_histos[25]] = ROOT.TH1F(list_histos[25],"n. of #gamma", 6, -0.5,5.5)
histo_map[list_histos[26]] = ROOT.TH1F(list_histos[26],"Efficiency steps", 7, 0.,7.)
histo_map[list_histos[27]] = ROOT.TH1F(list_histos[27],"Photon wp90 steps", 2, -0.5,1.5)
histo_map[list_histos[28]] = ROOT.TH1F(list_histos[28],"Phi or Rho channel", 2, -0.5,1.5)

#histo_map[list_histos[29]] = ROOT.TH1F(list_histos[29],"BDT output", 40, -1.,1.)



#CREATE OUTPUT ROOTFILE ##################################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#CREATE OUTPUT ROOTFILE
#fOutput = ROOT.TFile("histos/latest_production/histos_"+samplename+".root","RECREATE")
#fOutput.cd()

#Variables to go in the output tree #############################################################################################
mesonGammaMass       = np.zeros(1, dtype=float)
mesonMass            = np.zeros(1, dtype=float)
_coupleIsoCh         = np.zeros(1, dtype=float)
_coupleIso           = np.zeros(1, dtype=float)
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

tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('mesonGammaMass',mesonGammaMass,'mesonGammaMass/D')
tree_output.Branch('mesonMass',mesonMass,'mesonMass/D')
tree_output.Branch('_coupleIsoCh',_coupleIsoCh,'_coupleIsoCh/D')
tree_output.Branch('_coupleIso',_coupleIso,'_coupleIso/D')
tree_output.Branch('_bestJetPt',_bestJetPt,'_bestJetPt/D')
tree_output.Branch('_bestCouplePt',_bestCouplePt,'_bestCouplePt/D')
tree_output.Branch('_firstTrkPt',_firstTrkPt,'_firstTrkPt/D')
tree_output.Branch('_secondTrkPt',_secondTrkPt,'_secondTrkPt/D')
tree_output.Branch('_photonEt',_photonEt,'_photonEt/D')
tree_output.Branch('_firstTrkIso',_firstTrkIso,'_firstTrkIso/D')
tree_output.Branch('_firstTrkIsoCh',_firstTrkIso,'_firstTrkIsoCh/D')
tree_output.Branch('_secondTrkIso',_secondTrkIso,'_secondTrkIso/D')
tree_output.Branch('_secondTrkIsoCh',_secondTrkIso,'_secondTrkIsoCh/D')
tree_output.Branch('_firstTrkEta',_firstTrkEta,'_firstTrkEta/D')
tree_output.Branch('_secondTrkEta',_secondTrkEta,'_secondTrkEta/D')
tree_output.Branch('_bestCoupleEta',_bestCoupleEta,'_bestCoupleEta/D')
tree_output.Branch('_bestCoupleDeltaR',_bestCoupleDeltaR,'_bestCoupleDeltaR/D')
tree_output.Branch('_photonEta',_photonEta,'_photonEta/D')
tree_output.Branch('_eventWeight',_eventWeight,'_eventWeight/D')
        

#LOOSE SELECTION ########################################################################################################

#Higgs mass
higgs_min_invMass = 100.
higgs_max_invMass = 150.    
#Phi iso
iso_coupleCh_max  = 1.  
#Phi pT
meson_pT_min      = 30. 
#firstTrk pT
firstTrk_pT_min   = 20.
#secondTrk pT
secondTrk_pT_min  = 5.
#photon eT
photon_eT_min     = 35.
#best jet pT
bestJet_pT_max    = 150.
#photon WP90
photonWP90bool    = 1.

    
#counters
nEventsOverCuts = 0

#SELECT ALL BUT ONE FUNCTION ########################################################################################################
#dictionary: allows to apply all cuts except one related to the plotted variable
#if it returns FALSE it throws that event away
def select_all_but_one(h_string):

    if not doSelection:
        return True

    selection_bools = dict()
    selection_bools["h_InvMass_TwoTrk_Photon"] = mytree.Hmass_From2K_Photon >= higgs_min_invMass and mytree.Hmass_From2K_Photon <= higgs_max_invMass 
    selection_bools["h_couple_Iso_ch"]         = mytree.iso_couple_ch <= iso_coupleCh_max
    selection_bools["h_bestCouplePt"]          = mytree.bestCouplePt >= meson_pT_min
    selection_bools["h_firstTrk_pT"]           = mytree.firstCandPt >= firstTrk_pT_min
    selection_bools["h_secondTrk_pT"]          = mytree.secondCandPt >= secondTrk_pT_min
    selection_bools["h_photon_energy"]         = mytree.photon_eT >= photon_eT_min
    selection_bools["h_bestJetPt"]             = mytree.bestJet_pT <= bestJet_pT_max
    selection_bools["h_photonWP90"]            = mytree.is_photon_wp90 == photonWP90bool

    result = True

    for hname in selection_bools:
        if h_string == hname:
            continue
        else:
            result = result and selection_bools[hname]
               
    return result


#Photon scaling function ###########################################################################################################
def get_photon_scale(ph_pt, ph_eta):

    local_ph_pt = ph_pt
    if local_ph_pt > 499.: # This is because corrections are up to 499 GeV
        local_ph_pt = 499.
    
    local_ph_eta = ph_eta
    if local_ph_eta >= 2.5:
        local_ph_eta = 2.49
    if local_ph_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
        local_ph_eta = -2.5

    scale_factor_ID      = ph_ID_scale_histo_2018.GetBinContent( ph_ID_scale_histo_2018.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2018.GetYaxis().FindBin(local_ph_pt) )
    scale_factor_pixVeto = ph_pixVeto_scale_histo_2018.GetBinContent( ph_pixVeto_scale_histo_2018.GetXaxis().FindBin(local_ph_pt), ph_pixVeto_scale_histo_2018.GetYaxis().FindBin(math.fabs(local_ph_eta)) )
    scale_factor         = scale_factor_ID * scale_factor_pixVeto

    #if debug:
     #   print "local_ph_pt          = ",local_ph_pt
      #  print "local_ph_eta         = ",local_ph_eta
       # print "scale_factor_ID      = ",scale_factor_ID
       # print "scale_factor_pixVeto = ",scale_factor_pixVeto

    return scale_factor
    
print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

photonId_true  = 0
photonId_false = 0


#EVENTS LOOP ##################################################################################################################### 
for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        continue

    
    print "Processing EVENT n.",jentry+1,"================"

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


    if photonWP90:
        photonId_true += 1
    else:
        photonId_false += 1

    #If I'm performing a PhiGamma analysis I don't want to choose those events tagged as a RhoGamma events, and viceversa
    if isPhiAnalysis and not isPhiEvent: 
        continue
    if isRhoAnalysis and not isRhoEvent: 
        continue


    #Define Control and Signal regions: 
    if isPhiAnalysis: #for Phi meson
        if CRflag == 0 and not (MesonMass > 1.008 and MesonMass < 1.032) :
            continue

        if CRflag == 1 and (MesonMass > 1.008 and MesonMass < 1.032) :
            continue

    if isRhoAnalysis: #for Rho meson
        if CRflag == 0 and not (MesonMass > 0.63 and MesonMass < 0.91) :
            continue

        if CRflag == 1 and (MesonMass > 0.63 and MesonMass < 0.91) :
            continue

################### line used to take simmetrical sidebands ################
    if (isPhiAnalysis and MesonMass > 1.04): continue
############################################################################

    print"#################"
    print "isRhoAnalysis = ",isRhoAnalysis
    print "isRhoEvent = ",isRhoEvent
    print "CRflag = ",CRflag
    print "MesonMass = ",MesonMass
    print "isDataBlind = ",isDataBlind
    print"#################"


    #TIGHT SELECTION from BDT output -------------------------------------------------  
    if isBDT:
        BDT_out = myWF.get_BDT_output(firstTrkiso,MesonPt,photonEt,Hmass)
        #histo_map["h_BDT_out"].Fill(BDT_out)

        print "BDT value before selection = ", BDT_out
        if BDT_out < BDT_OUT: #Cut on BDT output
            print "BDT cut NOT passed"
            continue


    #NORMALIZATION -------------------------------------------------------------------
    #normalization for MC
    if not samplename == "Data":
        
        PUWeight          = mytree.PU_Weight
        weight_sign       = mytree.MC_Weight/abs(mytree.MC_Weight)
        photonScaleFactor = get_photon_scale(mytree.photon_eT, mytree.photon_eta)
        eventWeight       = weight_sign * luminosity * normalization_weight * PUWeight * photonScaleFactor
        
        if debug:
            print "luminosity = ",luminosity
            print "normalization_weight  = ",normalization_weight
            print "mytree.MC_Weight = ",mytree.MC_Weight
            print "weight sign = ",weight_sign
            print "PUWeight = ",PUWeight
            print "photonScaleFactor = ",photonScaleFactor
            print "eventWeight = ",eventWeight
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
   

    #FILL HISTOS --------------------------------------------------------------------------
    #if DATA -> Blind Analysis on H inv mass plot
    if select_all_but_one("h_InvMass_TwoTrk_Photon"):
        if samplename == "Data":
            if isDataBlind:
                if Hmass < 120. or Hmass > 130.:
                    histo_map["h_InvMass_TwoTrk_Photon"].Fill(Hmass, eventWeight)
            else:
                histo_map["h_InvMass_TwoTrk_Photon"].Fill(Hmass, eventWeight)
        else:
            histo_map["h_InvMass_TwoTrk_Photon"].Fill(Hmass, eventWeight)
 
    if select_all_but_one(""): #H inv mass plot without the phi mass cut    
        if samplename == "Data":
            if isDataBlind:
                if Hmass < 120. or Hmass > 130.:
                    histo_map["h_InvMass_TwoTrk_Photon_NoPhiMassCut"].Fill(Hmass, eventWeight)
            else:
                    histo_map["h_InvMass_TwoTrk_Photon_NoPhiMassCut"].Fill(Hmass, eventWeight)
        else:
            histo_map["h_InvMass_TwoTrk_Photon_NoPhiMassCut"].Fill(Hmass, eventWeight)
            
    if select_all_but_one("h_meson_InvMass_TwoTrk"): 
        histo_map["h_meson_InvMass_TwoTrk"].Fill(MesonMass, eventWeight)       
    if select_all_but_one("h_couple_Iso_ch"):    
        histo_map["h_couple_Iso_ch"].Fill(MesonIsoCh, eventWeight)
    if select_all_but_one("h_bestJetPt"):    
        histo_map["h_bestJetPt"].Fill(jetPt, eventWeight)
    if select_all_but_one("h_bestCouplePt"):    
        histo_map["h_bestCouplePt"].Fill(MesonPt, eventWeight)
    if select_all_but_one("h_firstTrk_pT"):    
        histo_map["h_firstTrk_pT"].Fill(firstTrkPt, eventWeight)
    if select_all_but_one("h_secondTrk_pT"):    
        histo_map["h_secondTrk_pT"].Fill(secondTrkPt, eventWeight)
    if select_all_but_one("h_photon_energy"): 
        histo_map["h_photon_energy"].Fill(photonEt, eventWeight)
    if select_all_but_one("h_photonWP90"): 
        histo_map["h_photonWP90"].Fill(photonWP90, eventWeight)

      
        
    if select_all_but_one(""):    
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
        histo_map["h_nMuons"].Fill(nMu, eventWeight)
        histo_map["h_nElectrons"].Fill(nEle, eventWeight)
        histo_map["h_nPhotons"].Fill(nPhotons, eventWeight)
        histo_map["h_bestCoupleEta"].Fill(MesonEta, eventWeight)
        histo_map["h_decayChannel"].Fill(isPhiEvent, eventWeight)


    #------------------------------------------------------------------------------------
    
    #FILL TREE
    if select_all_but_one(""):    
        mesonGammaMass[0]     = Hmass
        mesonMass[0]          = MesonMass
        _coupleIsoCh[0]       = MesonIsoCh
        _coupleIso[0]         = MesonIso
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
        tree_output.Fill()
    
        if debug:
            print "tree filled"


    #counters
    if select_all_but_one(""):
        nEventsOverCuts += 1

        if samplename == "Signal":
            sigWeightSum += _eventWeight


if debug:        
    print "#############################"
    print "wp90 photons = ",photonId_true
    print "NOT wp90 photons = ",photonId_false
    print "#############################"


#HISTO LABELS
histo_map["h_InvMass_TwoTrk_Photon"].GetXaxis().SetTitle("m_{trk^{+}trk^{-}#gamma} [GeV/c^2]")
histo_map["h_InvMass_TwoTrk_Photon"].SetTitle("Tracks+Photon invariant mass (Cut on phi inv. mass)")
histo_map["h_InvMass_TwoTrk_Photon_NoPhiMassCut"].GetXaxis().SetTitle("m_{trk^{+}trk^{-}}#gamma} [GeV/c^2]")
histo_map["h_InvMass_TwoTrk_Photon_NoPhiMassCut"].SetTitle("Tracks+Photon invariant mass (No cuts)")
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

#FINAL PRINTS ###########################################################
print ""
print "########### SUMMARY #############"
if isBDT:
    print "BDT output used = ",BDT_OUT
if samplename == "Signal":
    print "Total n. of MC signal events ~ 50000"
if isBDT:
    print "n. events after preselection = ",jentry
print "n. events after cuts = " , nEventsOverCuts
if samplename == "Signal":
    print "Signal weight sum = ",float(sigWeightSum)
    if isPhiAnalysis:
        print "Signal efficiency = ",nEventsOverCuts/49750.
    if isRhoAnalysis:
        print "Signal efficiency = ",nEventsOverCuts/50000.
print "#################################"
print ""
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
    histo_map["h_efficiency"].GetXaxis().SetRangeUser(1.,7.1)
    histo_map["h_efficiency"].GetYaxis().SetRangeUser(0.,30.)
    #histo_map["h_efficiency"].SetMaximum(max(histo_map["h_efficiency"].GetHistogram().GetMaximum(),30.))
    histo_map["h_efficiency"].Draw("HIST TEXT0")
    c11.SaveAs("plots/h_efficiency.pdf")
    c11.SaveAs("plots/h_efficiency.png")    


#HISTOS WRITING
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()
fOut.Close()
