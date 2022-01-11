import ROOT

#Take data rootfiles for control regions
fileCR1 = ROOT.TFile("histos/latest_production/histos_CR1_Data.root")
fileCR2 = ROOT.TFile("histos/latest_production/histos_CR2_Data.root")
fileCR3 = ROOT.TFile("histos/latest_production/histos_CR3_Data.root")

#Take histos corresponding to variables 
h_CR1_photEt   = fileCR1.Get("h_photon_energy")
h_CR2_photEt   = fileCR2.Get("h_photon_energy")
h_CR3_photEt   = fileCR3.Get("h_photon_energy")
h_CR1_2trk_iso = fileCR1.Get("h_couple_Iso_ch")
h_CR2_2trk_iso = fileCR2.Get("h_couple_Iso_ch")
h_CR3_2trk_iso = fileCR3.Get("h_couple_Iso_ch")

#Histos rebinning
h_CR1_photEt.Rebin(4)
h_CR2_photEt.Rebin(4)
h_CR3_photEt.Rebin(4)
h_CR1_2trk_iso.Rebin(100) #rebin 100 because it is the number of bins of iso histo and I want a monobin histo
h_CR2_2trk_iso.Rebin(100)
h_CR3_2trk_iso.Rebin(100)

#Histos dividing for transfer factor calculation
h_CR1_photEt.Divide(h_CR3_photEt)
h_CR2_photEt.Divide(h_CR3_photEt)
h_CR1_2trk_iso.Divide(h_CR3_2trk_iso)
h_CR2_2trk_iso.Divide(h_CR3_2trk_iso)

canva_eT = ROOT.TCanvas("canva_eT","",200,106,600,600)
canva_eT.Divide(2,1)
canva_eT.cd(1)
h_CR1_photEt.Draw("E1")
canva_eT.cd(2)
h_CR2_photEt.Draw("E1")

#canva_iso = ROOT.TCanvas("canva_iso","",200,106,600,600)
#canva_iso.Divide(2,1)
#canva_iso.cd(1)
#h_CR1_2trk_iso.Draw("E1")
#canva_iso.cd(2)
#h_CR2_2trk_iso.Draw("E1")

canva_eT.SaveAs("plots/CRfrac_photEt.gif")
#canva_iso.SaveAs("plots/CRfrac_2trkIso.gif")

#Change histos name
h_CR1_photEt.SetName("CR1_fraction_photET")
h_CR2_photEt.SetName("CR2_fraction_photET")
h_CR1_2trk_iso.SetName("CR1_fraction_2trkIso")
h_CR2_2trk_iso.SetName("CR2_fraction_2trkIso")

#Output file creation
fOut = ROOT.TFile("histos/latest_production/CRfraction.root","RECREATE")
fOut.cd()

#Histos writing
h_CR1_photEt.Write()
h_CR2_photEt.Write()
h_CR1_2trk_iso.Write()
h_CR2_2trk_iso.Write()

fOut.Close()
