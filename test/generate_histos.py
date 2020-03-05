import ROOT
import argparse
import math
#import numpy as np

pdf_flag=False
debug=False
doSelection=True

p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('rootfile_name', help='Type rootfile name')
args = p.parse_args()

#print "input = ",args.rootfile_name

fInput = ROOT.TFile(args.rootfile_name)

#SPLIT: I can split a string of chars before and after a split code (that could be a string or a symbol)
#then I can take the string stands before or after with [0] or [1], respectively. 
#[:-n] is not an emoji, it deletes last n symbols from the string
samplename =(args.rootfile_name.split("HPhiGammaAnalysis_")[1])[:-5] 
print "samplename =", samplename

if not samplename == "Data":
    normalization_InputFile = open("rootfiles/latest_production/MC/normalizations/Normalizations_table.txt","r")

mytree = fInput.Get("HPhiGammaAnalysis/mytree")

if not samplename == "Data":
    h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")

#normalization for MC dataset
if not samplename == "Data":
    norm_map = dict()
    for line in normalization_InputFile:
        data_norm = line.split()
        norm_map[data_norm[0]] = float(data_norm[1])

    MC_Weight = norm_map[samplename]

    if debug:
        print "MC_Weight = ",MC_Weight

    luminosity2018A = 14.00 #fb^-1
    luminosity2018B = 3.41 #fb^-1
    luminosity2018C = 6.94 #fb^-1
    luminosity2018D = 31.93 #fb^-1

    luminosity = luminosity2018B + luminosity2018D
    MC_Weight = MC_Weight * luminosity

    if debug:
        print "MC_Weight * lumi = ",MC_Weight

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

#plot: H invariant mass from 2K and photon
h_InvMass_TwoTrk_Photon = ROOT.TH1F("h_InvMass_TwoTrk_Photon","h_InvMass_TwoTrk_Photon",50,100.,150.)
histo_map["h_InvMass_TwoTrk_Photon"] = h_InvMass_TwoTrk_Photon

#plot: H invariant mass from 2K and photon with a cut on phi invariant mass
h_InvMass_TwoTrk_Photon_NoPhiMassCut = ROOT.TH1F("h_InvMass_TwoTrk_Photon_NoPhiMassCut","h_InvMass_TwoTrk_Photon_NoPhiMassCut",50,100.,150.)
histo_map["h_InvMass_TwoTrk_Photon_NoPhiMassCut"] =  h_InvMass_TwoTrk_Photon_NoPhiMassCut

#plot: PHI invariant mass from 2K 
h_phi_InvMass_TwoTrk = ROOT.TH1F("phi_InvMass_TwoTrk","phi_InvMass_TwoTrk", 200, 1., 1.05)
histo_map["h_phi_InvMass_TwoTrk"] = h_phi_InvMass_TwoTrk

#plot: first K-candidate pT (pT max of the couple)
h_firstKCand_pT = ROOT.TH1F("firstKCand_pT","firstKCand_pT", 100, 0.,150.)
histo_map["h_firstKCand_pT"] = h_firstKCand_pT

#plot: second K-candidate pT (pT min of the couple)
h_secondKCand_pT = ROOT.TH1F("secondKCand_pT","secondKCand_pT", 100, 0.,150.)
histo_map["h_secondKCand_pT"] = h_secondKCand_pT

#plot: first K-candidate eta
h_firstKCand_Eta = ROOT.TH1F("firstKCand_Eta","firstKCand_Eta", 100, -2.5,2.5)
histo_map["h_firstKCand_Eta"] = h_firstKCand_Eta

#plot: second K-candidate eta
h_secondKCand_Eta = ROOT.TH1F("secondKCand_Eta","secondKCand_Eta", 100, -2.5,2.5)
histo_map["h_secondKCand_Eta"] = h_secondKCand_Eta

#plot: first K-candidate phi
h_firstKCand_Phi = ROOT.TH1F("firstKCand_Phi","firstKCand_Phi", 100, -3.14,3.14)
histo_map["h_firstKCand_Phi"] = h_firstKCand_Phi

#plot: second K-candidate phi
h_secondKCand_Phi = ROOT.TH1F("secondKCand_Phi","secondKCand_Phi", 100, -3.14,3.14)
histo_map["h_secondKCand_Phi"] = h_secondKCand_Phi

#plot: best couple pt
h_bestCouplePt = ROOT.TH1F("bestCouplePt","bestCouplePt", 100, 0.,150.)
histo_map["h_bestCouplePt"] = h_bestCouplePt

#plot: best couple Eta
h_bestCoupleEta = ROOT.TH1F("bestCoupleEta","bestCoupleEta", 100, -2.5,2.5)
histo_map["h_bestCoupleEta"] = h_bestCoupleEta

#plot: best couple deltaR
h_bestCoupleDeltaR = ROOT.TH1F("bestCoupleDeltaR","bestCoupleDeltaR", 100, 0.,0.02)
histo_map["h_bestCoupleDeltaR"] = h_bestCoupleDeltaR

#plot: best jet pt
h_bestJetPt = ROOT.TH1F("bestJetPt","bestJetPt", 100, 0.,400.)
histo_map["h_bestJetPt"] = h_bestJetPt

#plot: best jet Eta
h_bestJetEta = ROOT.TH1F("bestJetEta","bestJetEta", 100, -2.5,2.5)
histo_map["h_bestJetEta"] = h_bestJetEta

#plot: K1 isolation
h_K1_Iso = ROOT.TH1F("K1_Iso","K1_Iso", 100, 0.,1.)
histo_map["h_K1_Iso"] = h_K1_Iso

#plot: K1 isolation ch
h_K1_Iso_ch = ROOT.TH1F("h_K1_Iso_ch","h_K1_Iso_ch", 100, 0.,1.)
histo_map["h_K1_Iso_ch"] = h_K1_Iso_ch

#plot: K2 isolation
h_K2_Iso = ROOT.TH1F("K2_Iso","K2_Iso", 100, 0.,1.)
histo_map["h_K2_Iso"] = h_K2_Iso

#plot: K2 isolation ch
h_K2_Iso_ch = ROOT.TH1F("h_K2_Iso_ch","h_K2_Iso_ch", 100, 0.,1.)
histo_map["h_K2_Iso_ch"] = h_K2_Iso_ch

#plot: couple isolation
h_couple_Iso = ROOT.TH1F("h_couple_Iso","h_couple_Iso", 100, 0.,1.)
histo_map["h_couple_Iso"] = h_couple_Iso

#plot: couple isolation ch
h_couple_Iso_ch = ROOT.TH1F("h_couple_Iso_ch","h_couple_Iso_ch", 100, 0.,1.)
histo_map["h_couple_Iso_ch"] = h_couple_Iso_ch

#plot: photon energy
h_photon_energy = ROOT.TH1F("h_photon_energy","h_photon_energy", 100, 0.,300.)
histo_map["h_photon_energy"] = h_photon_energy

#plot: photon eta
h_photon_eta = ROOT.TH1F("h_photon_eta","h_photon_eta", 100, -2.5,2.5)
histo_map["h_photon_eta"] = h_photon_eta

#plot: n. of jets over 25 GeV
h_nJets_25 = ROOT.TH1F("h_nJets_25","h_nJets_25", 13, -0.5,12.5)
histo_map["h_nJets_25"] = h_nJets_25

#plot: n. of muons for each event
h_nMuons = ROOT.TH1F("h_nMuons","h_nMuons", 11, -0.5,10.5)
histo_map["h_nMuons"] = h_nMuons

#plot: n. of electrons for each event
h_nElectrons = ROOT.TH1F("h_nElectrons","h_nElectrons", 6, -0.5,5.5)
histo_map["h_nElectrons"] = h_nElectrons

#plot: n. of photons for each event
h_nPhotons = ROOT.TH1F("h_nPhotons","h_nPhotons", 6, -0.5,5.5)
histo_map["h_nPhotons"] = h_nPhotons

"""
#CREATE OUTPUT ROOTFILE
#Variables to go in the tree
_HiggsMass = np.zeros(1, dtype=float)
_PhiMass = np.zeros(1, dtype=float)
_coupleIsoCh = np.zeros(1, dtype=float)
_bestJetPt = np.zeros(1, dtype=float)
_bestCouplePt = np.zeros(1, dtype=float)
_firstCandPt = np.zeros(1, dtype=float)
_secondCandPt = np.zeros(1, dtype=float)

tree_output = ROOT.TTree('Output Tree','tree with branches')
tree_output.Branch('Higgs Mass',_HiggsMass,'HiggsMass/D')
tree_output.Branch('Phi Mass',_PhiMass,'PhiMass/D')
tree_output.Branch('Couple Iso Ch',_coupleIsoCh,'CoupleIsoCh/D')
tree_output.Branch('Best Jet Pt',_bestJetPt,'bestJetPt/D')
tree_output.Branch('Best Couple Pt',_bestCouplePt,'bestCouplePt/D')
tree_output.Branch('First Candidate Pt',_firstCandPt,'PfirstCandt/D')
tree_output.Branch('Second Candidate Pt',_secondCandPt,'PsecondCandt/D')
"""

#CUTS
phi_min_invMass = 1.01
phi_max_invMass = 1.03
higgs_min_invMass = 100.
higgs_max_invMass = 150.

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
    selection_bools["h_couple_Iso_ch"] = mytree.iso_couple_ch <= 0.3
    selection_bools["h_bestJetPt"] = mytree.bestJet_pT <= 150.
    selection_bools["h_bestCouplePt"] = mytree.bestCouplePt >= 40.
    selection_bools["h_firstKCand_pT"] = mytree.firstCandPt >= 25.
    selection_bools["h_secondKCand_pT"] = mytree.secondCandPt >= 15.
    selection_bools["h_photon_energy"] = mytree.photon_eT >= 50.

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
        PUWeight = mytree.PU_Weight

        if debug:
            print "event n.",jentry+1
            print "PUWeight = ",PUWeight

        MC_Weight = abs(MC_Weight) #re-inizialize MC_Weight as positive
        MC_Weight *= mytree.MC_Weight/abs(mytree.MC_Weight)

        if debug:
            print "mytree.MC_Weight = ",mytree.MC_Weight
            print "mytree.MC_Weight/abs(mytree.MC_Weight) = ",mytree.MC_Weight/abs(mytree.MC_Weight)
            print "MC_Weight *= mytree.MC_Weight/abs(mytree.MC_Weight) = ",MC_Weight

        photonScaleFactor = get_photon_scale(mytree.photon_eT, mytree.photon_eta)
        eventWeight = PUWeight*MC_Weight*photonScaleFactor
        
        if debug:
            print "photonScaleFactor = ",photonScaleFactor
            print "PUWeight*MC_Weight*photonScaleFactor = eventWeight = ",eventWeight
            print ""

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
            if mytree.Hmass_From2K_Photon < 120. or mytree.Hmass_From2K_Photon > 130.:
                h_InvMass_TwoTrk_Photon.Fill(mytree.Hmass_From2K_Photon, eventWeight)
        else:
            h_InvMass_TwoTrk_Photon.Fill(mytree.Hmass_From2K_Photon, eventWeight)

    if select_all_but_one(""): #H inv mass plot without the phi mass cut    
        if samplename == "Data":
            if mytree.Hmass_From2K_Photon < 120. or mytree.Hmass_From2K_Photon > 130.:
                h_InvMass_TwoTrk_Photon_NoPhiMassCut.Fill(mytree.Hmass_From2K_Photon, eventWeight)
        else:
            h_InvMass_TwoTrk_Photon_NoPhiMassCut.Fill(mytree.Hmass_From2K_Photon, eventWeight)

    if select_all_but_one("h_phi_InvMass_TwoTrk"):
        h_phi_InvMass_TwoTrk.Fill(mytree.Phimass, eventWeight)
    
    if select_all_but_one("h_couple_Iso_ch"):    
        h_couple_Iso_ch.Fill(mytree.iso_couple_ch, eventWeight)

    if select_all_but_one("h_bestJetPt"):    
        h_bestJetPt.Fill(mytree.bestJet_pT, eventWeight)

    if select_all_but_one("h_bestCouplePt"):    
        h_bestCouplePt.Fill(mytree.bestCouplePt, eventWeight)

    if select_all_but_one("h_firstKCand_pT"):    
        h_firstKCand_pT.Fill(mytree.firstCandPt, eventWeight)

    if select_all_but_one("h_secondKCand_pT"):    
        h_secondKCand_pT.Fill(mytree.secondCandPt, eventWeight)
   
    if select_all_but_one("h_photon_energy"):    
        h_photon_energy.Fill(mytree.photon_eT, eventWeight)

    if select_all_but_one(""):    
        h_firstKCand_Eta.Fill(mytree.firstCandEta, eventWeight)    
    if select_all_but_one(""):    
        h_secondKCand_Eta.Fill(mytree.secondCandEta, eventWeight)   
    if select_all_but_one(""):    
        h_firstKCand_Phi.Fill(mytree.firstCandPhi, eventWeight)    
    if select_all_but_one(""):    
        h_secondKCand_Phi.Fill(mytree.secondCandPhi, eventWeight)   
    if select_all_but_one(""):    
        h_bestCoupleDeltaR.Fill(deltaR, eventWeight)
    if select_all_but_one(""):    
        h_bestJetEta.Fill(mytree.bestJet_eta, eventWeight)
    if select_all_but_one(""):    
        h_couple_Iso.Fill(mytree.iso_couple, eventWeight)
    if select_all_but_one(""):    
        h_K1_Iso.Fill(mytree.iso_K1, eventWeight)
    if select_all_but_one(""):    
        h_K1_Iso_ch.Fill(mytree.iso_K1_ch, eventWeight)
    if select_all_but_one(""):    
        h_K2_Iso.Fill(mytree.iso_K2, eventWeight)
    if select_all_but_one(""):    
        h_K2_Iso_ch.Fill(mytree.iso_K2_ch, eventWeight)
    if select_all_but_one(""):    
        h_photon_eta.Fill(mytree.photon_eta, eventWeight)
    if select_all_but_one(""):    
        h_nJets_25.Fill(mytree.nJets_25, eventWeight)
    if select_all_but_one(""):    
        h_nMuons.Fill(mytree.nMuons, eventWeight)
    if select_all_but_one(""):    
        h_nElectrons.Fill(mytree.nElectrons, eventWeight)
    if select_all_but_one(""):    
        h_nPhotons.Fill(mytree.nPhotonsOverSelection, eventWeight)
    if select_all_but_one(""):    
        h_bestCoupleEta.Fill(mytree.bestCoupleEta, eventWeight)



    #counters
    if select_all_but_one(""):
        nEventsOverCuts += 1

h_InvMass_TwoTrk_Photon.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^2]")
h_InvMass_TwoTrk_Photon.SetTitle("Tracks+Photon invariant mass (Cut on phi inv. mass)")

h_InvMass_TwoTrk_Photon_NoPhiMassCut.GetXaxis().SetTitle("m_{K^{+}K^{-}}#gamma} [GeV/c^2]")
h_InvMass_TwoTrk_Photon_NoPhiMassCut.SetTitle("Tracks+Photon invariant mass (No cuts)")

h_phi_InvMass_TwoTrk.GetXaxis().SetTitle("m_{K^{+}K^{-}} [GeV/c^2]")
h_phi_InvMass_TwoTrk.SetTitle("Tracks invariant mass")

h_firstKCand_pT.GetXaxis().SetTitle("pT_{K} [GeV/c]")
h_firstKCand_pT.SetTitle("Transverse momentum of the first charged particle (pT_{max} of the couple)")

h_secondKCand_pT.GetXaxis().SetTitle("pT_{K} [GeV/c]")
h_secondKCand_pT.SetTitle("Transverse momentum of the second charged particle (pT_{min} of the couple)")

h_firstKCand_Eta.GetXaxis().SetTitle("#eta")
h_firstKCand_Eta.SetTitle("Pseudorapidity of the first charged particle (pT_{max} of the couple)")

h_secondKCand_Eta.GetXaxis().SetTitle("#eta")
h_secondKCand_Eta.SetTitle("Pseudorapidity of the second charged particle (pT_{max} of the couple)")

h_firstKCand_Phi.GetXaxis().SetTitle("#phi [rad]")
h_firstKCand_Phi.SetTitle("Azimuthal angle of the first charged particle (pT_{max} of the couple)")

h_secondKCand_Phi.GetXaxis().SetTitle("#phi [rad]")
h_secondKCand_Phi.SetTitle("Azimuthal angle of the second charged particle (pT_{max} of the couple)")

h_bestCouplePt.GetXaxis().SetTitle("pT_{K^{+}K^{-}} [GeV]")
h_bestCouplePt.SetTitle("Transverse momentum of the couple")

h_bestCoupleEta.GetXaxis().SetTitle("#eta")
h_bestCoupleEta.SetTitle("Pseudorapidity of the couple")

h_bestJetPt.GetXaxis().SetTitle("pT_{jet} [GeV]")
h_bestJetPt.SetTitle("Transverse momentum of the jet")

h_bestJetEta.GetXaxis().SetTitle("#eta_{jet}")
h_bestJetEta.SetTitle("Pseudorapidity of the jet")

h_K1_Iso.GetXaxis().SetTitle("sum pT_{K_{1}}/pT_{K_{1}}")
h_K1_Iso.SetTitle("Isolation of the K_{1} candidate")

h_K1_Iso_ch.GetXaxis().SetTitle("sum pT_{K_{1}}/pT_{K_{1}}")
h_K1_Iso_ch.SetTitle("Isolation of the K_{1} candidate")

h_K2_Iso.GetXaxis().SetTitle("sum pT_{K_{2}}/pT_{K_{2}}")
h_K2_Iso.SetTitle("Isolation of the K_{2} candidate")

h_K2_Iso_ch.GetXaxis().SetTitle("sum pT_{K_{2}}/pT_{K_{2}}")
h_K2_Iso_ch.SetTitle("Isolation of the K_{1} candidate")

h_couple_Iso.GetXaxis().SetTitle("sum pT_{2K}/pT_{2K}")
h_couple_Iso.SetTitle("Isolation of the couple candidate")

h_couple_Iso_ch.GetXaxis().SetTitle("sum pT_{2K}/pT_{2K}")
h_couple_Iso_ch.SetTitle("Isolation of the couple candidate")

h_bestCoupleDeltaR.GetXaxis().SetTitle("#DeltaR_{K^{+}K^{-}}")
h_bestCoupleDeltaR.SetTitle("Delta R of the couple")

h_photon_energy.GetXaxis().SetTitle("E_{#gamma}")
h_photon_energy.SetTitle("Energy of the photon")

h_photon_eta.GetXaxis().SetTitle("#eta_{#gamma}")
h_photon_eta.SetTitle("Eta of the photon")

h_nJets_25.GetXaxis().SetTitle("nJets over 25 GeV")
h_nJets_25.SetTitle("# of jets with pT > 25 GeV")

h_nPhotons.GetXaxis().SetTitle("n.#gamma")
h_nPhotons.SetTitle("n.#gamma over selections")


fOutput = ROOT.TFile("histos/latest_production/histos_"+samplename+".root","RECREATE")

if samplename == "Signal":
    print "n. events after cuts: " , nEventsOverCuts

print "entries in histos= ",h_phi_InvMass_TwoTrk.GetEntries()
print ""

for histo in histo_map:
    histo_map[histo].Write(histo)
    if pdf_flag: 
        if samplename == "Signal" or samplename == "Data":
            canvas = ROOT.TCanvas()
            canvas.cd()
            histo_map[histo].Draw("E1")
            canvas.SaveAs("plots/" + histo +".pdf")

fOutput.Close()
