import ROOT
import argparse
import math
import numpy as np

pdf_flag=False
debug=False
doSelection=True
tightSelection = True
isDataBlind = False

p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('rootfile_name', help='Type rootfile name')
args = p.parse_args()

fInput = ROOT.TFile(args.rootfile_name)
h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")

#SPLIT: I can split a string of chars before and after a split code (that could be a string or a symbol)
#then I can take the string stands before or after with [0] or [1], respectively. 
#[:-n] is not an emoji, it deletes last n symbols from the string
samplename =(args.rootfile_name.split("HPhiGammaAnalysis_")[1])[:-5] 
print "samplename =", samplename

if not samplename == "Data":
    normalization_InputFile = open("rootfiles/latest_production/MC/normalizations/Normalizations_table.txt","r")
#    debug = False

mytree = fInput.Get("HPhiGammaAnalysis/mytree")

if not samplename == "Data":
    h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")

#normalization for MC dataset
if not samplename == "Data":
    norm_map = dict()
    for line in normalization_InputFile:
        data_norm = line.split()
        norm_map[data_norm[0]] = float(data_norm[1])

    normalization_weight = norm_map[samplename]

    luminosity2018A = 14.00 #fb^-1
    luminosity2018B = 3.41 #fb^-1
    luminosity2018C = 6.94 #fb^-1
    luminosity2018D = 31.93 #fb^-1

    #luminosity = luminosity2018B + luminosity2018C + luminosity2018D
    luminosity = 39.54 #fb^-1

#MCfromDATA correction for photons 
ph_ID_scale_name_2018  = "scale_factors/2018_PhotonsMVAwp90.root"
ph_ID_scale_file_2018  = ROOT.TFile(ph_ID_scale_name_2018)
ph_ID_scale_histo_2018 = ROOT.TH2F()
ph_ID_scale_histo_2018 = ph_ID_scale_file_2018.Get("EGamma_SF2D")

ph_pixVeto_scale_name_2018  = "scale_factors/HasPix_2018.root"
ph_pixVeto_scale_file_2018  = ROOT.TFile(ph_pixVeto_scale_name_2018)
ph_pixVeto_scale_histo_2018 = ROOT.TH1F()
ph_pixVeto_scale_histo_2018 = ph_pixVeto_scale_file_2018.Get("eleVeto_SF")

histo_map = dict()

list_histos = ["h_InvMass_TwoTrk_Photon","h_InvMass_TwoTrk_Photon_NoPhiMassCut","h_phi_InvMass_TwoTrk","h_firstKCand_pT","h_secondKCand_pT","h_firstKCand_Eta","h_secondKCand_Eta","h_firstKCand_Phi","h_secondKCand_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestCoupleDeltaR","h_bestJetPt","h_bestJetEta","h_K1_Iso","h_K1_Iso_ch","h_K2_Iso","h_K2_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons","h_efficiency"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{H}",90,80.,170.) 
histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{H} (no cut on phi mass)",90,80.,170.) 
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"M_{#phi}", 200, 1., 1.05) 
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 1st K", 100, 0.,150.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"p_{T} of the 2nd K", 100, 0.,150.)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta of the 1st K", 100, -2.5,2.5)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#eta of the 2nd K", 100, -2.5,2.5)
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi of the 1st K", 100, -3.14,3.14)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"#phi of the 2nd K", 100, -3.14,3.14)
histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"p_{T} of the #phi", 100, 0.,150.)
histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#eta_{#phi}", 100, -2.5,2.5)
histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"#Delta R_{#phi}", 100, 0.,0.02)
histo_map[list_histos[12]] = ROOT.TH1F(list_histos[12],"p_{T} of the jet", 100, 0.,400.)
histo_map[list_histos[13]] = ROOT.TH1F(list_histos[13],"#eta of the jet", 100, -2.5,2.5)
histo_map[list_histos[14]] = ROOT.TH1F(list_histos[14],"Iso of the 1st K", 100, 0.,1.)
histo_map[list_histos[15]] = ROOT.TH1F(list_histos[15],"Iso_ch of the 1st K", 100, 0.,1.)
histo_map[list_histos[16]] = ROOT.TH1F(list_histos[16],"Iso of the 2nd K", 100, 0.,1.)
histo_map[list_histos[17]] = ROOT.TH1F(list_histos[17],"Iso_ch of the 2nd K", 100, 0.,1.)
histo_map[list_histos[18]] = ROOT.TH1F(list_histos[18],"Iso of the #phi", 100, 0.,1.)
histo_map[list_histos[19]] = ROOT.TH1F(list_histos[19],"Iso_ch of the #phi", 100, 0.,1.)
histo_map[list_histos[20]] = ROOT.TH1F(list_histos[20],"E_{T} of the #gamma", 100, 0.,300.)
histo_map[list_histos[21]] = ROOT.TH1F(list_histos[21],"#eta_{#gamma}", 100, -2.5,2.5)
histo_map[list_histos[22]] = ROOT.TH1F(list_histos[22],"n. of jets over pre-filters", 13, -0.5,12.5)
histo_map[list_histos[23]] = ROOT.TH1F(list_histos[23],"n. of muons", 11, -0.5,10.5)
histo_map[list_histos[24]] = ROOT.TH1F(list_histos[24],"n. of electrons", 6, -0.5,5.5)
histo_map[list_histos[25]] = ROOT.TH1F(list_histos[25],"n. of #gamma", 6, -0.5,5.5)
histo_map[list_histos[26]] = ROOT.TH1F(list_histos[26],"Efficiency steps", 7, 0.,7.)

fOutput = ROOT.TFile("histos/latest_production/histos_"+samplename+".root","RECREATE")
fOutput.cd()


#CREATE OUTPUT ROOTFILE
#Variables to go in the tree
mass_KKg = np.zeros(1, dtype=float)
mass_KK = np.zeros(1, dtype=float)
_coupleIsoCh = np.zeros(1, dtype=float)
_bestJetPt = np.zeros(1, dtype=float)
_bestCouplePt = np.zeros(1, dtype=float)
_firstCandPt = np.zeros(1, dtype=float)
_secondCandPt = np.zeros(1, dtype=float)
_photonEt = np.zeros(1, dtype=float)
_eventWeight = np.zeros(1, dtype=float)

tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('mass_KKg',mass_KKg,'mass_KKg/D')
tree_output.Branch('mass_KK',mass_KK,'mass_KK/D')
tree_output.Branch('_coupleIsoCh',_coupleIsoCh,'_coupleIsoCh/D')
tree_output.Branch('_bestJetPt',_bestJetPt,'_bestJetPt/D')
tree_output.Branch('_bestCouplePt',_bestCouplePt,'_bestCouplePt/D')
tree_output.Branch('_firstCandPt',_firstCandPt,'_firstCandt/D')
tree_output.Branch('_secondCandPt',_secondCandPt,'_secondCandt/D')
tree_output.Branch('_photonEt',_photonEt,'_photonEt/D')
tree_output.Branch('_eventWeight',_eventWeight,'_eventWeight/D')

if samplename == "Data":
    print "Data tree output set!"  #FIXME


#CUTS
if tightSelection:
    #Higgs mass
    higgs_min_invMass = 80.
    higgs_max_invMass = 170.
    #Phi mass
    phi_min_invMass = 1.015 #1.015
    phi_max_invMass = 1.027 #1.027
    #Phi iso
    iso_coupleCh_max = 0.2  #0.2
    #Phi pT
    phi_pT_min = 42. #42. 
    #K1 pT
    firstCand_pT_min = 27.
    #K2 pT
    secondCand_pT_min = 20.
    #photon eT
    photon_eT_min = 35.
    #best jet pT
    bestJet_pT_max = 120. 
    
#counters
nEventsOverCuts = 0

#dictionary: allows to apply all cuts except one related to the plotted variable
#if it returns FALSE it throws that event away
def select_all_but_one(h_string):

    if not doSelection:
        return True

    selection_bools = dict()
    selection_bools["h_phi_InvMass_TwoTrk"] = mytree.Phimass >= phi_min_invMass and mytree.Phimass <= phi_max_invMass 
    selection_bools["h_InvMass_TwoTrk_Photon"] = mytree.Hmass_From2K_Photon >= higgs_min_invMass and mytree.Hmass_From2K_Photon <= higgs_max_invMass 
    selection_bools["h_couple_Iso_ch"] = mytree.iso_couple_ch <= iso_coupleCh_max
    selection_bools["h_bestCouplePt"] = mytree.bestCouplePt >= phi_pT_min
    selection_bools["h_firstKCand_pT"] = mytree.firstCandPt >= firstCand_pT_min
    selection_bools["h_secondKCand_pT"] = mytree.secondCandPt >= secondCand_pT_min
    selection_bools["h_photon_energy"] = mytree.photon_eT >= photon_eT_min
    selection_bools["h_bestJetPt"] = mytree.bestJet_pT <= bestJet_pT_max

    result = True

    for hname in selection_bools:
        if h_string == hname:
            continue
        else:
            result = result and selection_bools[hname]
               
    return result


#photon scaling
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
    scale_factor = scale_factor_ID * scale_factor_pixVeto

    if debug:
        print "local_ph_pt = ",local_ph_pt
        print "local_ph_eta = ",local_ph_eta
        print "scale_factor_ID = ",scale_factor_ID
        print "scale_factor_pixVeto = ",scale_factor_pixVeto

    return scale_factor
    
print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()


#EVENTS LOOP
for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        continue

    #normalization calculation for MC dataset
    if not samplename == "Data":

        if debug:
            print "EVENT n.",jentry+1,"================"

        PUWeight = mytree.PU_Weight

        weight_sign = mytree.MC_Weight/abs(mytree.MC_Weight)
            
        photonScaleFactor = get_photon_scale(mytree.photon_eT, mytree.photon_eta)

        eventWeight = weight_sign * luminosity * normalization_weight * PUWeight * photonScaleFactor
        
        if debug:
            print "luminosity = ",luminosity
            print "normalization_weight  = ",normalization_weight
            print "mytree.MC_Weight = ",mytree.MC_Weight
            print "weight sign = ",weight_sign
            print "PUWeight = ",PUWeight
            print "photonScaleFactor = ",photonScaleFactor
            print "eventWeight = ",eventWeight
            print ""

    #retrieving from the tree
    Hmass = mytree.Hmass_From2K_Photon
    PhiMass = mytree.Phimass
    PhiIsoCh = mytree.iso_couple_ch
    jetPt = mytree.bestJet_pT      
    PhiPt = mytree.bestCouplePt      
    firstKpT = mytree.firstCandPt       
    secondKpT = mytree.secondCandPt        
    photonEt = mytree.photon_eT        
    firstKeta = mytree.firstCandEta
    secondKeta = mytree.secondCandEta
    firstKphi = mytree.firstCandPhi
    secondKphi = mytree.secondCandPhi
    jetEta = mytree.bestJet_eta
    PhiIso = mytree.iso_couple
    firstKiso = mytree.iso_K1
    firstKisoCh = mytree.iso_K1_ch
    secondKiso = mytree.iso_K2
    secondKisoCh = mytree.iso_K2_ch
    photonEta  =mytree.photon_eta
    nJets = mytree.nJets_25
    nMu = mytree.nMuons
    nEle = mytree.nElectrons
    nPhotons = mytree.nPhotonsOverSelection
    PhiEta = mytree.bestCoupleEta


    #phi angle folding
    coupleDeltaPhi = math.fabs(mytree.firstCandPhi - mytree.secondCandPhi)
    if coupleDeltaPhi > 3.14:
        coupleDeltaPhi = 6.28 - coupleDeltaPhi
    deltaR = math.sqrt((mytree.firstCandEta - mytree.secondCandEta)**2 + (coupleDeltaPhi)**2)

    #NO normalization for DATA 
    if samplename == "Data":
        eventWeight = 1.   

    
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
            
    if select_all_but_one("h_phi_InvMass_TwoTrk"):
        histo_map["h_phi_InvMass_TwoTrk"].Fill(PhiMass, eventWeight)
        
    if select_all_but_one("h_couple_Iso_ch"):    
        histo_map["h_couple_Iso_ch"].Fill(PhiIsoCh, eventWeight)
        
    if select_all_but_one("h_bestJetPt"):    
        histo_map["h_bestJetPt"].Fill(jetPt, eventWeight)

    if select_all_but_one("h_bestCouplePt"):    
        histo_map["h_bestCouplePt"].Fill(PhiPt, eventWeight)

    if select_all_but_one("h_firstKCand_pT"):    
        histo_map["h_firstKCand_pT"].Fill(firstKpT, eventWeight)

    if select_all_but_one("h_secondKCand_pT"):    
        histo_map["h_secondKCand_pT"].Fill(secondKpT, eventWeight)

    if select_all_but_one("h_photon_energy"):    
        histo_map["h_photon_energy"].Fill(photonEt, eventWeight)

        
        
    if select_all_but_one(""):    
        histo_map["h_firstKCand_Eta"].Fill(firstKeta, eventWeight)    
    if select_all_but_one(""):    
        histo_map["h_secondKCand_Eta"].Fill(secondKeta, eventWeight)   
    if select_all_but_one(""):    
        histo_map["h_firstKCand_Phi"].Fill(firstKphi, eventWeight)    
    if select_all_but_one(""):    
        histo_map["h_secondKCand_Phi"].Fill(secondKphi, eventWeight)   
    if select_all_but_one(""):    
        histo_map["h_bestCoupleDeltaR"].Fill(deltaR, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_bestJetEta"].Fill(jetEta, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_couple_Iso"].Fill(PhiIso, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_K1_Iso"].Fill(firstKiso, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_K1_Iso_ch"].Fill(firstKisoCh, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_K2_Iso"].Fill(secondKiso, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_K2_Iso_ch"].Fill(secondKisoCh, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_photon_eta"].Fill(photonEta, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_nJets_25"].Fill(nJets, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_nMuons"].Fill(nMu, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_nElectrons"].Fill(nEle, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_nPhotons"].Fill(nPhotons, eventWeight)
    if select_all_but_one(""):    
        histo_map["h_bestCoupleEta"].Fill(PhiEta, eventWeight)


    if select_all_but_one(""):    
        mass_KKg[0] = Hmass
        mass_KK[0] = PhiMass
        _coupleIsoCh[0] = PhiIsoCh
        _bestJetPt[0] = jetPt        
        _bestCouplePt[0] = PhiPt        
        _firstCandPt[0] = firstKpT        
        _secondCandPt[0] = secondKpT           
        _photonEt[0] = photonEt        
        _eventWeight[0] = eventWeight
        tree_output.Fill()
        
        if samplename == "Data":
            print "event n.",jentry  #FIXME
            print "tree filled"
        
    #counters
    if select_all_but_one(""):
        nEventsOverCuts += 1
        
histo_map["h_InvMass_TwoTrk_Photon"].GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^2]")
histo_map["h_InvMass_TwoTrk_Photon"].SetTitle("Tracks+Photon invariant mass (Cut on phi inv. mass)")

histo_map["h_InvMass_TwoTrk_Photon_NoPhiMassCut"].GetXaxis().SetTitle("m_{K^{+}K^{-}}#gamma} [GeV/c^2]")
histo_map["h_InvMass_TwoTrk_Photon_NoPhiMassCut"].SetTitle("Tracks+Photon invariant mass (No cuts)")

histo_map["h_phi_InvMass_TwoTrk"].GetXaxis().SetTitle("m_{K^{+}K^{-}} [GeV/c^2]")
histo_map["h_phi_InvMass_TwoTrk"].SetTitle("Tracks invariant mass")

histo_map["h_firstKCand_pT"].GetXaxis().SetTitle("pT_{K} [GeV/c]")
histo_map["h_firstKCand_pT"].SetTitle("Transverse momentum of the first charged particle (pT_{max} of the couple)")

histo_map["h_secondKCand_pT"].GetXaxis().SetTitle("pT_{K} [GeV/c]")
histo_map["h_secondKCand_pT"].SetTitle("Transverse momentum of the second charged particle (pT_{min} of the couple)")

histo_map["h_firstKCand_Eta"].GetXaxis().SetTitle("#eta")
histo_map["h_firstKCand_Eta"].SetTitle("Pseudorapidity of the first charged particle (pT_{max} of the couple)")

histo_map["h_secondKCand_Eta"].GetXaxis().SetTitle("#eta")
histo_map["h_secondKCand_Eta"].SetTitle("Pseudorapidity of the second charged particle (pT_{max} of the couple)")

histo_map["h_firstKCand_Phi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_firstKCand_Phi"].SetTitle("Azimuthal angle of the first charged particle (pT_{max} of the couple)")

histo_map["h_secondKCand_Phi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_secondKCand_Phi"].SetTitle("Azimuthal angle of the second charged particle (pT_{max} of the couple)")

histo_map["h_bestCouplePt"].GetXaxis().SetTitle("pT_{K^{+}K^{-}} [GeV]")
histo_map["h_bestCouplePt"].SetTitle("Transverse momentum of the couple")

histo_map["h_bestCoupleEta"].GetXaxis().SetTitle("#eta")
histo_map["h_bestCoupleEta"].SetTitle("Pseudorapidity of the couple")

histo_map["h_bestJetPt"].GetXaxis().SetTitle("pT_{jet} [GeV]")
histo_map["h_bestJetPt"].SetTitle("Transverse momentum of the jet")

histo_map["h_bestJetEta"].GetXaxis().SetTitle("#eta_{jet}")
histo_map["h_bestJetEta"].SetTitle("Pseudorapidity of the jet")

histo_map["h_K1_Iso"].GetXaxis().SetTitle("sum pT_{K_{1}}/pT_{K_{1}}")
histo_map["h_K1_Iso"].SetTitle("Isolation of the K_{1} candidate")

histo_map["h_K1_Iso_ch"].GetXaxis().SetTitle("sum pT_{K_{1}}/pT_{K_{1}}")
histo_map["h_K1_Iso_ch"].SetTitle("Isolation of the K_{1} candidate")

histo_map["h_K2_Iso"].GetXaxis().SetTitle("sum pT_{K_{2}}/pT_{K_{2}}")
histo_map["h_K2_Iso"].SetTitle("Isolation of the K_{2} candidate")

histo_map["h_K2_Iso_ch"].GetXaxis().SetTitle("sum pT_{K_{2}}/pT_{K_{2}}")
histo_map["h_K2_Iso_ch"].SetTitle("Isolation of the K_{1} candidate")

histo_map["h_couple_Iso"].GetXaxis().SetTitle("sum pT_{2K}/pT_{2K}")
histo_map["h_couple_Iso"].SetTitle("Isolation of the couple candidate")

histo_map["h_couple_Iso_ch"].GetXaxis().SetTitle("sum pT_{2K}/pT_{2K}")
histo_map["h_couple_Iso_ch"].SetTitle("Isolation of the couple candidate")

histo_map["h_bestCoupleDeltaR"].GetXaxis().SetTitle("#DeltaR_{K^{+}K^{-}}")
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


#if samplename == "Signal":
print "n. events after cuts: " , nEventsOverCuts

print "entries in histos= ", histo_map["h_nPhotons"].GetEntries()
print "entries in histos= ", histo_map["h_K2_Iso"].GetEntries()

print ""

tree_output.Write()
#tree_output.Scan()

#EFFICIENCY STEP PLOT MANAGEMENT
bin1content = h_Events.GetBinContent(1)
bin2content = h_Events.GetBinContent(2)
bin3content = h_Events.GetBinContent(3)
bin4content = h_Events.GetBinContent(4)
bin5content = h_Events.GetBinContent(5)
bin6content = h_Events.GetBinContent(6)
bin7content = nEventsOverCuts
nSignal = bin1content
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
histo_map["h_efficiency"].GetXaxis().SetBinLabel(5,"K-cand pT selection")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(6,"KK iso selection")
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


#HISTOS WRITING
for histo in histo_map:
    histo_map[histo].Write(histo)
    if pdf_flag: 
        if samplename == "Signal" or samplename == "Data":
            canvas = ROOT.TCanvas()
            canvas.cd()
            histo_map[histo].Draw("E1")
            canvas.SaveAs("plots/" + histo +".pdf")

fOutput.Close()
