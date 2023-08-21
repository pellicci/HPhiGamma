import ROOT
import tdrstyle, CMS_lumi

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)  

fileNominal     = ROOT.TFile.Open("histos/systematics/histos_SR_BDTcat0_SignalggH_PhotonSigma0.root")
fileVariationUP = ROOT.TFile.Open("histos/systematics/histos_SR_BDTcat0_SignalggH_PhotonSigmaUP.root")
fileVariationDW = ROOT.TFile.Open("histos/systematics/histos_SR_BDTcat0_SignalggH_PhotonSigmaDW.root")

syst_name = "Photon energy smearing variation"
pdf_name  = "photonEnSmear"

list_histos = []
keylist = fileNominal.GetListOfKeys()
key = ROOT.TKey()
for key in keylist :
    obj_class = ROOT.gROOT.GetClass(key.GetClassName())
    if not obj_class.InheritsFrom("TH1") :
        continue
    if not (key.ReadObj().GetName() == "h_efficiency" or key.ReadObj().GetName() == "h_cutOverflow"): #h_efficiency and h_cutOverflow is a plot plotted in other way   
        list_histos.append( key.ReadObj().GetName() )

list_histos = ["h_InvMass_TwoTrk_Photon"]#,"h_photon_energy"]

HmassMin = 115.
HmassMax = 135.

for histo_name in list_histos:

	nominal     = fileNominal.Get(histo_name)
	variationUP = fileVariationUP.Get(histo_name)
	variationDW = fileVariationDW.Get(histo_name)

	# Creazione del canvas e divisone in 2 parti
	canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
	canvas.Divide(1, 2)

	#CMS-style plotting 
	tdrstyle.setTDRStyle()

	# Primo pad
	canvas.cd(1)
	ROOT.gStyle.SetOptStat(0)
	ROOT.gPad.SetMargin(0.15, 0.05, 0.01, 0.1)  # margini L, R, B, T

	nominal.SetLineColor(ROOT.kBlack)
	variationUP.SetLineColor(ROOT.kBlue)
	variationDW.SetLineColor(ROOT.kRed)

	nominal.SetTitle(syst_name)
	variationUP.SetTitle(syst_name)
	variationDW.SetTitle(syst_name)

	if histo_name == "h_InvMass_TwoTrk_Photon":
		nominal.GetXaxis().SetRangeUser(HmassMin,HmassMax)
		variationUP.GetXaxis().SetRangeUser(HmassMin,HmassMax)
		variationDW.GetXaxis().SetRangeUser(HmassMin,HmassMax)

	nominal.Draw()
	variationUP.Draw("same")
	variationDW.Draw("same")

	ROOT.gPad.RedrawAxis("g")
	ROOT.gPad.Update()

	# Creazione della legenda
	legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
	legend.AddEntry(0,pdf_name,"")
	legend.AddEntry(nominal, "nominal", "l")
	legend.AddEntry(variationUP, "variationUP", "l")
	legend.AddEntry(variationDW, "variationDW", "l")
	legend.Draw()

	# Secondo pad
	canvas.cd(2)
	ROOT.gPad.SetMargin(0.15, 0.05, 0.3, 0.01)  # margini L, R, B, T

	# Calcolo dei ratio e disegno
	ratioUP = variationUP.Clone()
	ratioUP.Add(nominal, -1)  # sottrae nominal da variationUP
	ratioUP.Divide(nominal)   # divide per nominal
	ratioUP.SetLineColor(ROOT.kBlue)
	#ratioUP.GetXaxis().SetTitle("Invariant Mass [GeV]")
	ratioUP.GetYaxis().SetTitle("(var - nominal) / nominal")
	if histo_name == "h_InvMass_TwoTrk_Photon": ratioUP.GetXaxis().SetRangeUser(HmassMin,HmassMax)

	ratioUP.Draw()

	ratioDW = variationDW.Clone()
	ratioDW.Add(nominal, -1)  # sottrae nominal da variationDW
	ratioDW.Divide(nominal)   # divide per nominal
	ratioDW.SetLineColor(ROOT.kRed)
	if histo_name == "h_InvMass_TwoTrk_Photon": ratioDW.GetXaxis().SetRangeUser(HmassMin,HmassMax)
	ratioDW.Draw("same")

	# Aumenta la dimensione dei marker
	# Imposta lo stile e il colore dei marker
	nominal.SetMarkerStyle(20)
	nominal.SetMarkerColor(ROOT.kBlack)

	variationUP.SetMarkerStyle(20)
	variationUP.SetMarkerColor(ROOT.kBlue)

	variationDW.SetMarkerStyle(20)
	variationDW.SetMarkerColor(ROOT.kRed)

	ratioUP.SetMarkerStyle(20)
	ratioUP.SetMarkerColor(ROOT.kBlue)

	ratioDW.SetMarkerStyle(20)
	ratioDW.SetMarkerColor(ROOT.kRed)
	#nominal.SetMarkerSize(1.5)
	#variationUP.SetMarkerSize(1.5)
	#variationDW.SetMarkerSize(1.5)
	#ratioUP.SetMarkerSize(1.5)
	#ratioDW.SetMarkerSize(1.5)

	# Disegno di una linea orizzontale centrata su 0
	line = ROOT.TLine(ratioDW.GetXaxis().GetXmin(), 0, ratioDW.GetXaxis().GetXmax(), 0)  # i valori HmassMin e HmassMax sono l'inizio e la fine dell'intervallo x
	if histo_name == "h_InvMass_TwoTrk_Photon":line = ROOT.TLine(HmassMin, 0, HmassMax, 0)  # i valori HmassMin e HmassMax sono l'inizio e la fine dell'intervallo x
	line.SetLineColor(ROOT.kBlack)
	line.SetLineWidth(2)
	line.Draw()

	ROOT.gPad.RedrawAxis("g")
	ROOT.gPad.Update()

	canvas.SaveAs("/eos/user/g/gumoret/www/latest_production/systematic_uncertainties/comparison_signal_"+histo_name+"_"+pdf_name+".png")
	canvas.SaveAs("/eos/user/g/gumoret/www/latest_production/systematic_uncertainties/comparison_signal_"+histo_name+"_"+pdf_name+".pdf")
