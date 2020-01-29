import ROOT

#fInput = ROOT.TFile("HPhiGammaAnalysis_output.root")
#fInput = ROOT.TFile("rootfiles/HPhiGamma_Signal.root")
fInput = ROOT.TFile("rootfiles/latest_production/HPhiGammaAnalysis_output_all.root")
mytree = fInput.Get("HPhiGammaAnalysis/mytree")

cross_section = 48.58*10**(-12)
branching_ratio = 10.**(-5)
luminosity = 59.76*10**(15)
nEvents = 27654.
normalizationWeight = cross_section*branching_ratio*luminosity*(1/nEvents) 


print "normalization weight = ", normalizationWeight


#plot: H invariant mass from 2K and photon
h_InvMass_TwoTrk_Photon = ROOT.TH1F("h_InvMass_TwoTrk_Photon","h_InvMass_TwoTrk_Photon",50,100.,150.)

#plot: H invariant mass from 2K and photon with a cut on phi invariant mass
h_InvMass_TwoTrk_Photon_NoPhiMassCut = ROOT.TH1F("h_InvMass_TwoTrk_Photon_NoPhiMassCut","h_InvMass_TwoTrk_Photon_NoPhiMassCut",50,100.,150.)

#plot: PHI invariant mass from 2K 
h_phi_InvMass_TwoTrk = ROOT.TH1F("phi_InvMass_TwoTrk","phi_InvMass_TwoTrk", 200, 0., 1.3)

#plot: first K-candidate pT (pT max of the couple)
h_firstKCand_pT = ROOT.TH1F("firstKCand_pT","firstKCand_pT", 100, 0.,150.)

#plot: second K-candidate pT (pT min of the couple)
h_secondKCand_pT = ROOT.TH1F("secondKCand_pT","secondKCand_pT", 100, 0.,150.)

#plot: first K-candidate eta
h_firstKCand_Eta = ROOT.TH1F("firstKCand_Eta","firstKCand_Eta", 100, -2.5,2.5)

#plot: second K-candidate eta
h_secondKCand_Eta = ROOT.TH1F("secondKCand_Eta","secondKCand_Eta", 100, -2.5,2.5)

#plot: first K-candidate phi
h_firstKCand_Phi = ROOT.TH1F("firstKCand_Phi","firstKCand_Phi", 100, -3.14,3.14)

#plot: second K-candidate phi
h_secondKCand_Phi = ROOT.TH1F("secondKCand_Phi","secondKCand_Phi", 100, -3.14,3.14)

#plot: best couple pt
h_bestCouplePt = ROOT.TH1F("bestCouplePt","bestCouplePt", 100, 0.,150.)


#cuts
phi_min_invMass = 0.3
phi_max_invMass = 1.2
higgs_min_invMass = 100.
higgs_max_invMass = 150.

#dictionary: allows to apply all cuts except one related to the plotted variable
def select_all_but_one(h_string):

    selection_bools = dict()
    selection_bools["h_phi_InvMass_TwoTrk"] = mytree.Phimass >= phi_min_invMass and mytree.Phimass <= phi_max_invMass 
    selection_bools["h_InvMass_TwoTrk_Photon"] = mytree.Hmass_From2K_Photon >= higgs_min_invMass and mytree.Hmass_From2K_Photon <= higgs_max_invMass 
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
    MCWeight = mytree.MC_Weight/abs(mytree.MC_Weight)
    eventWeight = normalizationWeight*PUWeight*MCWeight

    print "Event n.", jentry, "  normWeight = ", round(normalizationWeight,7), "  PUWeight = ",PUWeight, "  MCWeight = ",MCWeight 

    if select_all_but_one(h_InvMass_TwoTrk_Photon.GetName()):
        h_InvMass_TwoTrk_Photon.Fill(mytree.Hmass_From2K_Photon, eventWeight)
        
    if select_all_but_one(h_phi_InvMass_TwoTrk.GetName()):
        h_phi_InvMass_TwoTrk.Fill(mytree.Phimass, eventWeight)
    
    h_InvMass_TwoTrk_Photon_NoPhiMassCut.Fill(mytree.Hmass_From2K_Photon, eventWeight)
    h_firstKCand_pT.Fill(mytree.firstCandPt, eventWeight)    
    h_secondKCand_pT.Fill(mytree.secondCandPt, eventWeight)   
    h_firstKCand_Eta.Fill(mytree.firstCandEta, eventWeight)    
    h_secondKCand_Eta.Fill(mytree.secondCandEta, eventWeight)   
    h_firstKCand_Phi.Fill(mytree.firstCandPhi, eventWeight)    
    h_secondKCand_Phi.Fill(mytree.secondCandPhi, eventWeight)   
    h_bestCouplePt.Fill(mytree.bestCouplePt, eventWeight)


h_InvMass_TwoTrk_Photon.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^2]")
h_InvMass_TwoTrk_Photon.SetTitle("Tracks+Photon invariant mass (Cut on phi inv. mass)")

h_InvMass_TwoTrk_Photon_NoPhiMassCut.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^2]")
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

h_bestCouplePt.GetXaxis().SetTitle("#pT_{K^{+}K^{-} [GeV]")
h_bestCouplePt.SetTitle("Transverse momentum of the couple")


c1 = ROOT.TCanvas()
c1.cd()
h_InvMass_TwoTrk_Photon.Draw("E1")
c1.SaveAs("plots/h_InvMass_TwoTrk_Photon.pdf")

c2 = ROOT.TCanvas()
c2.cd()
h_phi_InvMass_TwoTrk.Draw("E1")
c2.SaveAs("plots/h_phi_InvMass_TwoTrk.pdf")

c3 = ROOT.TCanvas()
c3.cd()
h_firstKCand_pT.Draw("E1")
c3.SaveAs("plots/h_firstKCand_pT.pdf")

c4 = ROOT.TCanvas()
c4.cd()
h_secondKCand_pT.Draw("E1")
c4.SaveAs("plots/h_secondKCand_pT.pdf")

c5 = ROOT.TCanvas()
c5.cd()
h_firstKCand_Eta.Draw("E1")
c5.SaveAs("plots/h_firstKCand_Eta.pdf")

c6 = ROOT.TCanvas()
c6.cd()
h_secondKCand_Eta.Draw("E1")
c6.SaveAs("plots/h_secondKCand_Eta.pdf")

c7 = ROOT.TCanvas()
c7.cd()
h_firstKCand_Phi.Draw("E1")
c7.SaveAs("plots/h_firstKCand_Phi.pdf")

c8 = ROOT.TCanvas()
c8.cd()
h_secondKCand_Phi.Draw("E1")
c8.SaveAs("plots/h_secondKCand_Phi.pdf")

c9 = ROOT.TCanvas()
c9.cd()
h_InvMass_TwoTrk_Photon_NoPhiMassCut.Draw("E1")
c9.SaveAs("plots/h_InvMass_TwoTrk_Photon_NoPhiMassCut.pdf")

c10 = ROOT.TCanvas()
c10.cd()
h_bestCouplePt.Draw("E1")
c10.SaveAs("plots/h_bestCouplePt.pdf")
