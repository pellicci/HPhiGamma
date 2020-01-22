import ROOT

fInput = ROOT.TFile("rootfiles/HPhiGamma_Signal.root")
mytree = fInput.Get("HPhiGammaAnalysis/mytree")

h_InvMass_TwoTrk_Photon = ROOT.TH1F("h_InvMass_TwoTrk_Photon","h_InvMass_TwoTrk_Photon",50,100.,150.)

nentries = mytree.GetEntriesFast()

print "This sample has ", mytree.GetEntriesFast(), " events"

for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        continue

    h_InvMass_TwoTrk_Photon.Fill(mytree.Hmass_From2K_Photon)



h_InvMass_TwoTrk_Photon.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^2]")

h_InvMass_TwoTrk_Photon.SetTitle("Tracks+Photon invariant mass")


c1 = ROOT.TCanvas()
c1.cd()
h_InvMass_TwoTrk_Photon.Draw("E1")
c1.SaveAs("plots/h_InvMass_TwoTrk_Photon.pdf")

