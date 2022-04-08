import ROOT

#Take data rootfiles for control regions
DataSR = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
Sidebands = ROOT.TFile("histos/latest_production/histos_CR_Sidebands.root")

#Output file creation
fOut = ROOT.TFile("histos/latest_production/histos_CR_SidebandsNorm.root","RECREATE")
fOut.cd()

list_histos = ["h_InvMass_TwoTrk_Photon","h_InvMass_TwoTrk_Photon_NoPhiMassCut","h_phi_InvMass_TwoTrk","h_firstKCand_pT","h_secondKCand_pT","h_firstKCand_Eta","h_secondKCand_Eta","h_firstKCand_Phi","h_secondKCand_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestCoupleDeltaR","h_bestJetPt","h_bestJetEta","h_K1_Iso","h_K1_Iso_ch","h_K2_Iso","h_K2_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons","h_photonWP90","h_couple_AbsIsoCh"]#,"h_BDT_out"]

for histo_name in list_histos:
	
	print "###############"
	print "histo_name = ",histo_name

	histoSR = DataSR.Get(histo_name)
	histoCR = Sidebands.Get(histo_name)

	if histo_name == "h_InvMass_TwoTrk_Photon" or histo_name == "h_InvMass_TwoTrk_Photon_NoPhiMassCut":

		CRintegral = histoCR.Integral() - histoCR.Integral(histoCR.GetXaxis().FindBin(120.),histoCR.GetXaxis().FindBin(130.)) #since in this plot there is the blind window for data in SR, this trick is to make the divide properly. Remember to bypass it for the unblinding

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
