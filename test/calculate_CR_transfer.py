import ROOT

#Take data rootfiles for control regions
fileCR1 = ROOT.TFile("histos/latest_production/histos_CR1_Data.root")
fileCR2 = ROOT.TFile("histos/latest_production/histos_CR2_Data.root")
fileCR3 = ROOT.TFile("histos/latest_production/histos_CR3_Data.root")

#Take histos corresponding to variables 
h_CR1_photWP   = fileCR1.Get("h_photonWP90")
h_CR2_photWP   = fileCR2.Get("h_photonWP90")
h_CR3_photWP   = fileCR3.Get("h_photonWP90")

h_CR1_2trk_iso = fileCR1.Get("h_couple_AbsIsoCh")
h_CR2_2trk_iso = fileCR2.Get("h_couple_AbsIsoCh")
h_CR3_2trk_iso = fileCR3.Get("h_couple_AbsIsoCh")

#Histos rebinning
h_CR1_photWP.Rebin(1)
h_CR2_photWP.Rebin(1)
h_CR3_photWP.Rebin(1)
h_CR1_2trk_iso.Rebin(4)
h_CR2_2trk_iso.Rebin(4)
h_CR3_2trk_iso.Rebin(4)

#Histos dividing for transfer factor calculation
h_CR1_photWP.Divide(h_CR3_photWP)
h_CR2_photWP.Divide(h_CR3_photWP)
h_CR1_2trk_iso.Divide(h_CR3_2trk_iso)
h_CR2_2trk_iso.Divide(h_CR3_2trk_iso)

canva = ROOT.TCanvas("canva","",200,106,600,600)
canva.Divide(2,1)
canva.cd(1)
h_CR1_photWP.Draw("E1")
canva.cd(2)
h_CR2_photWP.Draw("E1")

canva.SaveAs("plots/CRfrac.gif")

#Change histos name
h_CR1_photWP.SetName("CR1_fraction_photWP")
h_CR2_photWP.SetName("CR2_fraction_photWP")
h_CR1_2trk_iso.SetName("CR1_fraction_2trkIso")
h_CR2_2trk_iso.SetName("CR2_fraction_2trkIso")

#Output file creation
fOut = ROOT.TFile("histos/latest_production/CRfraction.root","RECREATE")
fOut.cd()

#Histos writing
h_CR1_photWP.Write()
h_CR2_photWP.Write()
h_CR1_2trk_iso.Write()
h_CR2_2trk_iso.Write()

fOut.Close()
