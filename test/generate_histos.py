import ROOT
import argparse
import math

pdf_flag=False

p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('rootfile_name', help='Type rootfile name')
args = p.parse_args()

print "input = ",args.rootfile_name

fInput = ROOT.TFile(args.rootfile_name)
normalization_InputFile = open("rootfiles/latest_production/MC/normalizations/Normalizations_table.txt","r")
mytree = fInput.Get("HPhiGammaAnalysis/mytree")
h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")

#SPLIT: I can split a string of chars before and after a split code (that could be a string or a symbol)
#then I can take the string stands before or after with [0] or [1], respectively. 
#[:-n] is not an emoji, it deletes last n symbols from the string
samplename =(args.rootfile_name.split("HPhiGammaAnalysis_")[1])[:-5] 
print "samplename =", samplename

norm_map = dict()
for line in normalization_InputFile:
    data_norm = line.split()
    norm_map[data_norm[0]] = float(data_norm[1])

MC_Weight = norm_map[samplename]

luminosity = 59.76 #fb^-1 (in crate_normalization_table.py there's a factor 1000)
MC_Weight = MC_Weight * luminosity


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

#plot: photon energy
h_nJets_25 = ROOT.TH1F("h_nJets_25","h_nJets_25", 15, -0.5,15.)
histo_map["h_nJets_25"] = h_nJets_25

#cuts
phi_min_invMass = 1.
phi_max_invMass = 1.05
higgs_min_invMass = 100.
higgs_max_invMass = 150.

#dictionary: allows to apply all cuts except one related to the plotted variable
def select_all_but_one(h_string):

    selection_bools = dict()
    selection_bools["h_phi_InvMass_TwoTrk"] = mytree.Phimass >= phi_min_invMass and mytree.Phimass <= phi_max_invMass 
    selection_bools["h_InvMass_TwoTrk_Photon"] = mytree.Hmass_From2K_Photon >= higgs_min_invMass and mytree.Hmass_From2K_Photon <= higgs_max_invMass 
    selection_bools["h_couple_Iso_ch"] = mytree.iso_couple_ch <= 0.3
    selection_bools["h_bestJetPt"] = mytree.bestJet_pT <= 200.
    selection_bools["h_bestCouplePt"] = mytree.bestCouplePt >= 40.
    selection_bools["h_secondKCand_pT"] = mytree.secondCandPt >= 15.
    selection_bools["h_photon_energy"] = mytree.photon_energy >= 50.

    result = True

    for hname in selection_bools:
        if h_string == hname:
            continue
        else:
            result = result and selection_bools[hname]
        return result


print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        continue

    PUWeight = mytree.PU_Weight
    MC_Weight *= mytree.MC_Weight/abs(mytree.MC_Weight)
    eventWeight = PUWeight*MC_Weight
    coupleDeltaPhi = math.fabs(mytree.firstCandPhi - mytree.secondCandPhi)
    if coupleDeltaPhi > 3.14:
        coupleDeltaPhi = 6.28 - coupleDeltaPhi
    deltaR = math.sqrt((mytree.firstCandEta - mytree.secondCandEta)**2 + (coupleDeltaPhi)**2)

    if select_all_but_one(h_InvMass_TwoTrk_Photon.GetName()):
        h_InvMass_TwoTrk_Photon.Fill(mytree.Hmass_From2K_Photon, eventWeight)
        
    if select_all_but_one(h_phi_InvMass_TwoTrk.GetName()):
        h_phi_InvMass_TwoTrk.Fill(mytree.Phimass, eventWeight)

    if select_all_but_one(h_couple_Iso_ch.GetName()):    
        h_couple_Iso_ch.Fill(mytree.iso_couple_ch, eventWeight)

   
    if select_all_but_one(h_bestJetPt.GetName()):    
        h_bestJetPt.Fill(mytree.bestJet_pT, eventWeight)

    if select_all_but_one(h_bestCouplePt.GetName()):    
        h_bestCouplePt.Fill(mytree.bestCouplePt, eventWeight)

    if select_all_but_one(h_secondKCand_pT.GetName()):    
        h_secondKCand_pT.Fill(mytree.secondCandPt, eventWeight)
   
    if select_all_but_one(h_photon_energy.GetName()):    
        h_photon_energy.Fill(mytree.photon_energy, eventWeight)

    if select_all_but_one(""):    
        h_InvMass_TwoTrk_Photon_NoPhiMassCut.Fill(mytree.Hmass_From2K_Photon, eventWeight)
    if select_all_but_one(""):    
        h_firstKCand_pT.Fill(mytree.firstCandPt, eventWeight)    
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
        h_nJets_25.Fill(mytree.nJets_25, eventWeight)

    
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

h_nJets_25.GetXaxis().SetTitle("nJets over 25 GeV")
h_nJets_25.SetTitle("# of jets with pT > 25 GeV")


fOutput = ROOT.TFile("histos/latest_production/histos_"+samplename+".root","RECREATE")


for histo in histo_map:
    histo_map[histo].Write(histo)
    if pdf_flag and samplename == "Signal":
        canvas = ROOT.TCanvas()
        canvas.cd()
        histo_map[histo].Draw("E1")
        canvas.SaveAs("plots/" + histo +".pdf")

#h_InvMass_TwoTrk_Photon.Write("h_InvMass_TwoTrk_Photon_PhiMassCut")
#h_phi_InvMass_TwoTrk.Write("h_phi_InvMass_TwoTrk")
#h_InvMass_TwoTrk_Photon_NoPhiMassCut.Write("h_InvMass_TwoTrk_Photon_NoPhiMassCut")
#h_firstKCand_pT.Write("h_firstKCand_pT")
#h_secondKCand_pT.Write("h_secondKCand_pT")
#h_firstKCand_Eta.Write("h_firstKCand_Eta")
#h_secondKCand_Eta.Write("h_secondKCand_Eta")
#h_firstKCand_Phi.Write("h_firstKCand_Phi")
#h_secondKCand_Phi.Write("h_secondKCand_Phi")
#h_bestCouplePt.Write("h_bestCouplePt")
#h_bestJetPt.Write("h_bestJetPt")
#h_bestJetEta.Write("h_bestJetEta")
#h_K1_Iso.Write("h_K1_Iso")
#h_K1_Iso_ch.Write("h_K1_Iso_ch")
#h_K2_Iso.Write("h_K2_Iso")
#h_K2_Iso_ch.Write("h_K2_Iso_ch")
#h_couple_Iso.Write("h_couple_Iso")
#h_couple_Iso_ch.Write("h_couple_Iso_ch")
#h_Events.Write("h_Events")

fOutput.Close()


