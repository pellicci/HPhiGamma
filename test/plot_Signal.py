import ROOT
import math

#fInput = ROOT.TFile("HPhiGammaAnalysis_output.root")
#fInput = ROOT.TFile("rootfiles/HPhiGamma_Signal.root")
fInput = ROOT.TFile("rootfiles/200130/secondRun/HPhiGammaAnalysis_Signal.root")
mytree = fInput.Get("HPhiGammaAnalysis/mytree")
h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")

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
h_phi_InvMass_TwoTrk = ROOT.TH1F("phi_InvMass_TwoTrk","phi_InvMass_TwoTrk", 200, 0.9, 1.3)

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

#plot: best couple deltaR
h_bestCoupleDeltaR = ROOT.TH1F("bestCoupleDeltaR","bestCoupleDeltaR", 100, 0.,0.22)

#plot: best jet pt
h_bestJetPt = ROOT.TH1F("bestJetPt","bestJetPt", 100, 0.,400.)

#plot: best jet Eta
h_bestJetEta = ROOT.TH1F("bestJetEta","bestJetEta", 100, -2.5,2.5)

#plot: K1 isolation
h_K1_Iso = ROOT.TH1F("K1_Iso","K1_Iso", 100, 0.,1.)

#plot: K1 isolation ch
h_K1_Iso_ch = ROOT.TH1F("h_K1_Iso_ch","h_K1_Iso_ch", 100, 0.,1.)

#plot: K2 isolation
h_K2_Iso = ROOT.TH1F("K2_Iso","K2_Iso", 100, 0.,1.)

#plot: K2 isolation ch
h_K2_Iso_ch = ROOT.TH1F("h_K2_Iso_ch","h_K2_Iso_ch", 100, 0.,1.)

#plot: couple isolation
h_couple_Iso = ROOT.TH1F("couple_Iso","couple_Iso", 100, 0.,1.)

#plot: couple isolation ch
h_couple_Iso_ch = ROOT.TH1F("h_couple_Iso_ch","h_couple_Iso_ch", 100, 0.,1.)



#cuts
phi_min_invMass = 0.9
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
    deltaR = math.sqrt((mytree.firstCandEta - mytree.secondCandEta)**2 + (mytree.firstCandPhi - mytree.secondCandPhi)**2)

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
    h_bestCoupleDeltaR.Fill(deltaR, eventWeight)
    h_bestJetPt.Fill(mytree.bestJet_pT, eventWeight)
    h_bestJetEta.Fill(mytree.bestJet_eta, eventWeight)
    h_K1_Iso.Fill(mytree.iso_K1, eventWeight)
    h_K1_Iso_ch.Fill(mytree.iso_K1_ch, eventWeight)
    h_K2_Iso.Fill(mytree.iso_K2, eventWeight)
    h_K2_Iso_ch.Fill(mytree.iso_K2_ch, eventWeight)
    h_couple_Iso.Fill(mytree.iso_couple, eventWeight)
    h_couple_Iso_ch.Fill(mytree.iso_couple_ch, eventWeight)

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

h_bestCoupleDeltaR.GetXaxis().SetTitle("#DeltaR_{K^{+}K^{-}}")
h_bestCoupleDeltaR.SetTitle("Delta R of the couple")

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

c11 = ROOT.TCanvas()
c11.cd()
h_Events.Draw("E1")
c11.SaveAs("plots/h_Events.pdf")

c12 = ROOT.TCanvas()
c12.cd()
h_bestJetPt.Draw("E1")
c12.SaveAs("plots/h_bestJetPt.pdf")

c13 = ROOT.TCanvas()
c13.cd()
h_bestJetEta.Draw("E1")
c13.SaveAs("plots/h_bestJetEta.pdf")

c14 = ROOT.TCanvas()
c14.cd()
h_K1_Iso.Draw("E1")
c14.SaveAs("plots/h_K1_Iso.pdf")

c15 = ROOT.TCanvas()
c15.cd()
h_K1_Iso_ch.Draw("E1")
c15.SaveAs("plots/h_K1_Iso_ch.pdf")

c16 = ROOT.TCanvas()
c16.cd()
h_K2_Iso.Draw("E1")
c16.SaveAs("plots/h_K2_Iso.pdf")

c17 = ROOT.TCanvas()
c17.cd()
h_K2_Iso_ch.Draw("E1")
c17.SaveAs("plots/h_K2_Iso_ch.pdf")

c18 = ROOT.TCanvas()
c18.cd()
h_couple_Iso.Draw("E1")
c18.SaveAs("plots/h_couple_Iso.pdf")

c19 = ROOT.TCanvas()
c19.cd()
h_couple_Iso_ch.Draw("E1")
c19.SaveAs("plots/h_couple_Iso_ch.pdf")

c20 = ROOT.TCanvas()
c20.cd()
h_bestCoupleDeltaR.Draw("E1")
c20.SaveAs("plots/h_bestCoupleDeltaR.pdf")
