import ROOT
import argparse
import tdrstyle, CMS_lumi

#BOOL
doBiasStudy = True

#This script does the fit of a combined PDF (signal pdf + bkg pdf) to a generated dataset

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('Decay_channel_option', help='Type <<Phi>> for Phi, <<Rho>> for Rho') #flag for bkg estimation
args = p.parse_args()

if args.Decay_channel_option == "Phi":
	isPhiGammaAnalysis = True
	CHANNEL = "Phi"
	print "H -> PhiGamma analysis"
else: 
    isRhoGammaAnalysis = True
    CHANNEL = "Rho"
    print "H -> RhoGamma analysis"


fileInput_ggH = ROOT.TFile("histos/latest_production/histos_SR_SignalggH.root")
fileInput_VBF = ROOT.TFile("histos/latest_production/histos_SR_SignalVBF.root")


#Some useful things before starting --------------------------------------------------------------------
ROOT.gROOT.SetBatch(True) #Supress the opening of many Canvas's
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+") #import the Doube Crystal Ball PDF

#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

#Define the observable --------------------------
mass = ROOT.RooRealVar("mesonGammaMass","mesonGammaMass",100.,170.,"GeV")

#Data input -------------------------------------------------------------
fileInputData = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
fileInputData.cd()
mytreeData = fileInputData.Get("tree_output")

#Signal input -------------------------------------------------------------
fileInputSignal_ggH = ROOT.TFile("histos/latest_production/histos_SR_SignalggH.root")
fileInputSignal_ggH.cd()
mytreeSignal_ggH = fileInputSignal_ggH.Get("tree_output")

fileInputSignal_VBF = ROOT.TFile("histos/latest_production/histos_SR_SignalVBF.root")
fileInputSignal_VBF.cd()
mytreeSignal_VBF = fileInputSignal_VBF.Get("tree_output")

#Import from workspaces --------------------------------------------------------------------
fInputSignalPDF = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_2018.root") #there's only one ws for both ggH and VBF 
fInputSignalPDF.cd()
workspaceSignal = fInputSignalPDF.Get("workspace_STAT_"+CHANNEL+"_GFcat_2018")

print "################################ WORKSPACE SIGNAL #############################################################"
workspaceSignal.Print() #print some info of the workspace content
print "###############################################################################################################"
print ""

#Signal PDF
workspaceSignal.var("dCB_nL_"+CHANNEL+"_GFcat_ggH").setConstant(1) 
workspaceSignal.var("dCB_nR_"+CHANNEL+"_GFcat_ggH").setConstant(1)
workspaceSignal.var("dCB_aL_"+CHANNEL+"_GFcat_ggH").setConstant(1)
workspaceSignal.var("dCB_aR_"+CHANNEL+"_GFcat_ggH").setConstant(1)
workspaceSignal.var("dCB_pole_"+CHANNEL+"_GFcat_ggH").setConstant(1)
workspaceSignal.var("dCB_width_"+CHANNEL+"_GFcat_ggH").setConstant(1)

workspaceSignal.var("dCB_nL_"+CHANNEL+"_GFcat_VBF").setConstant(1) 
workspaceSignal.var("dCB_nR_"+CHANNEL+"_GFcat_VBF").setConstant(1)
workspaceSignal.var("dCB_aL_"+CHANNEL+"_GFcat_VBF").setConstant(1)
workspaceSignal.var("dCB_aR_"+CHANNEL+"_GFcat_VBF").setConstant(1)
workspaceSignal.var("dCB_pole_"+CHANNEL+"_GFcat_VBF").setConstant(1)
workspaceSignal.var("dCB_width_"+CHANNEL+"_GFcat_VBF").setConstant(1)

signalPDF_ggH = workspaceSignal.pdf("crystal_ball_"+CHANNEL+"_GFcat_ggH") #SIGNAL PDF from the workspace
signalPDF_VBF = workspaceSignal.pdf("crystal_ball_"+CHANNEL+"_GFcat_VBF") #SIGNAL PDF from the workspace

#Background PDF
fInputSidebandsPDF = ROOT.TFile("workspaces/ws_sidebands.root")
fInputSidebandsPDF.cd()
workspaceSidebands   = fInputSidebandsPDF.Get("myworkspace")
bkgPDF_chebychev     = workspaceSidebands.pdf("chebychev_GFcat_bkg") #BKG PDF from the workspace
bkgPDF_exponential   = workspaceSidebands.pdf("exponential_GFcat_bkg") #BKG PDF from the workspace
bkgPDF_bernstein     = workspaceSidebands.pdf("bernstein_GFcat_bkg") #BKG PDF from the workspace
bkgPDF_polynomial    = workspaceSidebands.pdf("polynomial_GFcat_bkg") #BKG PDF from the workspace

print "############################### WORKSPACE BKG ##############################################################"
workspaceSidebands.Print()
print "############################################################################################################"
print ""

#Define the POI -------------------------------------------------------------------------------------
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",10**(-5),-0.1,10**(-2)) #parameter of interest
#B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",0.) #parameter of interest

#Number of events in m_KKg
fInput_ggH = ROOT.TFile("rootfiles/latest_production/MC/signals/HPhiGammaAnalysis_Signal_"+CHANNEL+"_ggH.root")
mytree_ggH         = fInput_ggH.Get("HPhiGammaAnalysis/mytree")
h_Events_ggH       = fInput_ggH.Get("HPhiGammaAnalysis/h_Events")
NsigPassed_ggH     = mytreeSignal_ggH.GetEntriesFast()
totalSigEvents_ggH = h_Events_ggH.GetBinContent(1)
sigEfficiency_ggH  = float(NsigPassed_ggH/totalSigEvents_ggH)

fInput_VBF = ROOT.TFile("rootfiles/latest_production/MC/signals/HPhiGammaAnalysis_Signal_"+CHANNEL+"_VBF.root")
mytree_VBF    = fInput_VBF.Get("HPhiGammaAnalysis/mytree")
h_Events_VBF  = fInput_VBF.Get("HPhiGammaAnalysis/h_Events")
NsigPassed_VBF = mytreeSignal_VBF.GetEntriesFast()
totalSigEvents_VBF = h_Events_VBF.GetBinContent(1)
sigEfficiency_VBF = float(NsigPassed_VBF/totalSigEvents_VBF)

Ndata = mytreeData.GetEntriesFast()

#N bkg estimation
Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata, 10., 10000.)

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
cross_sig_ggH = ROOT.RooRealVar("cross_sig","The signal cross section",46.87)#46.87pb (gluon fusion)
cross_sig_VBF = ROOT.RooRealVar("cross_sig","The signal cross section",3.78)#3.78pb (VBF)
lumi = ROOT.RooRealVar("lumi","The luminosity",39540)#pb^-1
eff_ggH  = ROOT.RooRealVar("eff_ggH","efficiency ggH",sigEfficiency_ggH) # = (n sig over tightselection)/49750
eff_VBF  = ROOT.RooRealVar("eff_VBF","efficiency VBF",sigEfficiency_VBF) # = (n sig over tightselection)/49750
if CHANNEL == "Phi": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",0.489) # = BR(phi -> K+K-)
if CHANNEL == "Rho": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",1.) # = BR(rho -> PiPi)
Nsig_ggH = ROOT.RooFormulaVar("Nsig_ggH","@0*@1*@2*@3*@4",ROOT.RooArgList(lumi,cross_sig_ggH,eff_ggH,B_R_meson,B_R_))
Nsig_VBF = ROOT.RooFormulaVar("Nsig_VBF","@0*@1*@2*@3*@4",ROOT.RooArgList(lumi,cross_sig_VBF,eff_VBF,B_R_meson,B_R_))

#ADD PDFs --------------------------------------------------------------------------------------------------
argListN = ROOT.RooArgList(Nsig_ggH,Nsig_VBF,Nbkg)
print "arg list N size = ",argListN.getSize()
argListPDFs_chebychev   = ROOT.RooArgList(signalPDF_ggH,signalPDF_VBF,bkgPDF_chebychev)
argListPDFs_exponential = ROOT.RooArgList(signalPDF_ggH,signalPDF_VBF,bkgPDF_exponential)
argListPDFs_bernstein   = ROOT.RooArgList(signalPDF_ggH,signalPDF_VBF,bkgPDF_bernstein)
argListPDFs_polynomial  = ROOT.RooArgList(signalPDF_ggH,signalPDF_VBF,bkgPDF_polynomial)

print "arg lists created!"
totPDF_chebychev   = ROOT.RooAddPdf("totPDF_chebychev","Cheby",argListPDFs_chebychev,argListN)
totPDF_exponential = ROOT.RooAddPdf("totPDF_exponential","Exp",argListPDFs_exponential,argListN)
totPDF_bernstein   = ROOT.RooAddPdf("totPDF_bernstein","Bern",argListPDFs_bernstein,argListN)
totPDF_polynomial  = ROOT.RooAddPdf("totPDF_polynomial","Pol",argListPDFs_polynomial,argListN)

print "PDFs added!"

#Generate dataset for blind analysis -----------------------------------------------------------------------------------
datasetGenerated_chebychev   = totPDF_chebychev.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_exponential = totPDF_exponential.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_bernstein   = totPDF_bernstein.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_polynomial  = totPDF_polynomial.generate(ROOT.RooArgSet(mass),Ndata)

print "Dataset generated!"

#FIT
totPDF_chebychev.fitTo(datasetGenerated_chebychev)
#totPDF_exponential.fitTo(datasetGenerated_exponential)
totPDF_bernstein.fitTo(datasetGenerated_bernstein)
#totPDF_polynomial.fitTo(datasetGenerated_polynomial)

#PLOT
xframe_chebychev = mass.frame(28)
xframe_chebychev.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_chebychev.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
datasetGenerated_chebychev.plotOn(xframe_chebychev)
totPDF_chebychev.plotOn(xframe_chebychev)

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe_chebychev.Draw()
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_chebychev.pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_chebychev.png")

'''
xframe_exponential = mass.frame(28)
xframe_exponential.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_exponential.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
datasetGenerated_exponential.plotOn(xframe_exponential)
totPDF_exponential.plotOn(xframe_exponential)

c2 = ROOT.TCanvas()
c2.cd()
c2.SetTitle("")
xframe_exponential.Draw()
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_exponential.pdf")
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_exponential.png")
'''
xframe_bernstein = mass.frame(28)
xframe_bernstein.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_bernstein.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
datasetGenerated_bernstein.plotOn(xframe_bernstein)
totPDF_bernstein.plotOn(xframe_bernstein)

c3 = ROOT.TCanvas()
c3.cd()
c3.SetTitle("")
xframe_bernstein.Draw()
c3.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_bernstein.pdf")
c3.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_bernstein.png")
'''
xframe_polynomial = mass.frame(28)
xframe_polynomial.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_polynomial.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
datasetGenerated_polynomial.plotOn(xframe_polynomial)
totPDF_polynomial.plotOn(xframe_polynomial)

c4 = ROOT.TCanvas()
c4.cd()
c4.SetTitle("")
xframe_polynomial.Draw()
c4.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_polynomial.pdf")
c4.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_polynomial.png")
'''

#Final useful information for debugging
print ""
print "######################################################"
print "Ndata          = ", Ndata
print "totalSigEvents _ggH = ",totalSigEvents_ggH
print "NsigPassed ggH    = ",NsigPassed_ggH
print "Signal efficiency after tightselection ggH = ", sigEfficiency_ggH
print "totalSigEvents _VBF = ",totalSigEvents_VBF
print "NsigPassed VBF    = ",NsigPassed_VBF
print "Signal efficiency after tightselection VBF = ", sigEfficiency_VBF
print "######################################################"
print ""

#create Workspace
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totPDF_chebychev)
getattr(workspace,'import')(datasetGenerated_chebychev)
#getattr(workspace,'import')(totPDF_exponential)
#getattr(workspace,'import')(datasetGenerated_exponential)
getattr(workspace,'import')(totPDF_bernstein)
getattr(workspace,'import')(datasetGenerated_bernstein)
#getattr(workspace,'import')(totPDF_polynomial)

fOut = ROOT.TFile("workspaces/ws_data_blind.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()


if doBiasStudy:
	
	pdfList = [totPDF_chebychev,totPDF_bernstein]#,totPDF_chebychev]#,totPDF_bernstein]

	multicanvas = ROOT.TCanvas()
	multicanvas.cd()
	multicanvas.Divide(len(pdfList),len(pdfList))
	canvas_index = 0

	for fitPDF in pdfList:
		for genPDF in pdfList:
			canvas_index += 1

			mcstudy = ROOT.RooMCStudy(genPDF, ROOT.RooArgSet(mass), ROOT.RooFit.Silence(),ROOT.RooFit.FitModel(fitPDF), ROOT.RooFit.Extended(1), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))
			mcstudy.generateAndFit(1000)
	
			#set the frame

			leg2 = ROOT.TLegend(0.12,0.75,0.44,0.97) #left positioning
			leg2.SetHeader(" ")
			leg2.SetNColumns(1)
			leg2.SetFillColorAlpha(0,0.)
			leg2.SetBorderSize(0)
			leg2.SetLineColor(1)
			leg2.SetLineStyle(1)
			leg2.SetLineWidth(1)
			leg2.SetFillStyle(1001)
			leg2.AddEntry(0,"gen pdf: "+genPDF.GetTitle(),"")
			leg2.AddEntry(0,"fit pdf: "+fitPDF.GetTitle(),"")

			BRpull_frame = mcstudy.plotPull(B_R_, ROOT.RooFit.Bins(28), ROOT.RooFit.FitGauss(1))
			BRpull_frame.SetTitle("")
			BRpull_frame.SetTitleOffset(1.5,"y")
			BRpull_frame.SetXTitle("Gen: "+genPDF.GetTitle()+" Vs Fit: "+fitPDF.GetTitle())
			BRpull_frame.SetMaximum(1.1*BRpull_frame.GetMaximum())
			
			multicanvas.cd(canvas_index)
			BRpull_frame.Draw()
			leg2.Draw()


	multicanvas.Draw()

	multicanvas.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull.png")
	multicanvas.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull.pdf")

	#NLLframe.SetTitle("")

	#Actually plot
	#c4 = ROOT.TCanvas()
#	c4.cd()
#	BRval_frame.Draw()
#	c4.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRval.pdf")
#	c4.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRval.png")#

#	c5 = ROOT.TCanvas()
#	c5.cd()
#	BRerr_frame.Draw()
#	c5.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRerr.pdf")
#	c5.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRerr.png")#

#	c6 = ROOT.TCanvas()
#	c6.cd()
#	BRpull_frame.Draw()
#	c6.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull.pdf")
#	c6.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull.png")#

#	c7 = ROOT.TCanvas()
#	c7.cd()
#	NLLframe.Draw()
#	c7.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/NLL.pdf")
#	c7.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/NLL.png")

