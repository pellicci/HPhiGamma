import ROOT

#Take data rootfiles for control regions
DataSR = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
Sidebands = ROOT.TFile("histos/latest_production/histos_CR_Sidebands.root")

#Output file creation
fOut = ROOT.TFile("histos/latest_production/histos_CR_SidebandsNorm.root","RECREATE")
fOut.cd()

list_histos = ["h_InvMass_TwoTrk_Photon","h_InvMass_TwoTrk_Photon_NoPhiMassCut","h_meson_InvMass_TwoTrk","h_firstTrk_pT","h_secondTrk_pT","h_firstTrk_Eta","h_secondTrk_Eta","h_firstTrk_Phi","h_secondTrk_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestCoupleDeltaR","h_bestJetPt","h_bestJetEta","h_firstTrk_Iso","h_firstTrk_Iso_ch","h_secondTrk_Iso","h_secondTrk_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons","h_photonWP90","h_decayChannel","h_couple_Iso_neutral"]#,"h_BDT_out"]

for histo_name in list_histos:
	
	print "###############"
	print "histo_name = ",histo_name

	histoSR = DataSR.Get(histo_name)
	histoCR = Sidebands.Get(histo_name)

	if histo_name == "h_InvMass_TwoTrk_Photon" or histo_name == "h_InvMass_TwoTrk_Photon_NoPhiMassCut":

		CRintegral = histoCR.Integral() - histoCR.Integral(histoCR.GetXaxis().FindBin(115.),histoCR.GetXaxis().FindBin(135.)) #since in this plot there is the blind window for data in SR, this trick is to make the divide properly. Remember to bypass it for the unblinding

	else:
		CRintegral = histoCR.Integral()

	SRintegral = histoSR.Integral()	
	
	print "histo SR integral = ", SRintegral	
	print "histo CR integral = ", CRintegral		

	histoCR.Scale(SRintegral/CRintegral)
	histoCR.Write()

	print "histo CR integral after the normalization = ", histoCR.Integral()		
	print "###############"

fOut.Close()
