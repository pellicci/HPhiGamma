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
p.add_argument('Sample', help='Type <<ggH>> for ggH, <<VBF>> for VBF')
args = p.parse_args()

if args.Decay_channel_option == "Phi":
	isPhiGammaAnalysis = True
	CHANNEL = "Phi"
	print "H -> PhiGamma analysis"
else: 
    isRhoGammaAnalysis = True
    CHANNEL = "Rho"
    print "H -> RhoGamma analysis"

if args.Sample == "ggH":
	SAMPLE = "ggH"
	fileInput = ROOT.TFile("histos/latest_production/histos_SR_SignalggH.root")
	print "ggH production mechanism"
else: 
	SAMPLE = "VBF"
	fileInput = ROOT.TFile("histos/latest_production/histos_SR_SignalVBF.root")
	print "VBF production mechanism"

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
fileInputSignal = ROOT.TFile("histos/latest_production/histos_SR_Signal"+SAMPLE+".root")
fileInputSignal.cd()
mytreeSignal = fileInputSignal.Get("tree_output")

#fileInputSignal_VBF = ROOT.TFile("histos/latest_production/histos_SR_SignalVBF.root")
#fileInputSignal_VBF.cd()
#mytreeSignal_VBF = fileInputSignal_VBF.Get("tree_output")

#Import from workspaces --------------------------------------------------------------------
fInputSignalPDF = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_2018.root")
fInputSignalPDF.cd()
workspaceSignal = fInputSignalPDF.Get("workspace_STAT_"+CHANNEL+"_GFcat_"+SAMPLE+"_2018")
print "################################ WORKSPACE SIGNAL #############################################################"
workspaceSignal.Print() #print some info of the workspace content
print "###############################################################################################################"
print ""

#Signal PDF
workspaceSignal.var("dCB_nL_"+CHANNEL+"_GFcat_"+SAMPLE).setConstant(1) 
workspaceSignal.var("dCB_nR_"+CHANNEL+"_GFcat_"+SAMPLE).setConstant(1)
workspaceSignal.var("dCB_aL_"+CHANNEL+"_GFcat_"+SAMPLE).setConstant(1)
workspaceSignal.var("dCB_aR_"+CHANNEL+"_GFcat_"+SAMPLE).setConstant(1)
workspaceSignal.var("dCB_pole_"+CHANNEL+"_GFcat_"+SAMPLE).setConstant(1)
workspaceSignal.var("dCB_width_"+CHANNEL+"_GFcat_"+SAMPLE).setConstant(1)

signalPDF1 = workspaceSignal.pdf("crystal_ball_"+CHANNEL+"_GFcat_"+SAMPLE) #SIGNAL PDF from the workspace
signalPDF2 = workspaceSignal.pdf("crystal_ball_"+CHANNEL+"_GFcat_"+SAMPLE) #SIGNAL PDF from the workspace

#Background PDF
fInputSidebandsPDF = ROOT.TFile("workspaces/ws_sidebands.root")
fInputSidebandsPDF.cd()
workspaceSidebands   = fInputSidebandsPDF.Get("myworkspace")
bkgPDF_chebychev     = workspaceSidebands.pdf("chebychev_GFcat_bkg") #BKG PDF from the workspace
bkgPDF_exponential   = workspaceSidebands.pdf("exponential_GFcat_bkg") #BKG PDF from the workspace
print "############################### WORKSPACE BKG ##############################################################"
workspaceSidebands.Print()
print "############################################################################################################"
print ""

#Define the POI -------------------------------------------------------------------------------------
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",10**(-5),-0.1,10**(-2)) #parameter of interest

#Number of events in m_KKg
fInput = ROOT.TFile("rootfiles/latest_production/MC/signals/HPhiGammaAnalysis_Signal_"+CHANNEL+"_"+SAMPLE+".root")
mytree         = fInput.Get("HPhiGammaAnalysis/mytree")
h_Events       = fInput.Get("HPhiGammaAnalysis/h_Events")
NsigPassed     = mytreeSignal.GetEntriesFast()
totalSigEvents = h_Events.GetBinContent(1)
sigEfficiency  = float(NsigPassed/totalSigEvents)

#mytree_VBF    = fInput_VBF.Get("HPhiGammaAnalysis/mytree")
#h_Events_VBF  = fInput_VBF.Get("HPhiGammaAnalysis/h_Events")
#NsigPassed_VBF = mytreeSignal_VBF.GetEntriesFast()
#totalSigEvents_VBF = h_Events_VBF.GetBinContent(1)
#sigEfficiency_VBF = float(NsigPassed_VBF/totalSigEvents_VBF)

Ndata = mytreeData.GetEntriesFast()

#N bkg estimation
Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata, 10., 10000.)

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
if SAMPLE == "ggH": cross_sig = ROOT.RooRealVar("cross_sig","The signal cross section",46.87)#46.87pb (gluon fusion)
if SAMPLE == "VBF": cross_sig = ROOT.RooRealVar("cross_sig","The signal cross section",3.78)#3.78pb (VBF)
lumi = ROOT.RooRealVar("lumi","The luminosity",39540)#pb^-1
eff  = ROOT.RooRealVar("eff","efficiency",sigEfficiency) # = (n sig over tightselection)/49750
if CHANNEL == "Phi": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",0.489) # = BR(phi -> K+K-)
if CHANNEL == "Rho": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",1.) # = BR(rho -> PiPi)
#B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",10**(-5),-0.01,10**(-2))
Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3*@4",ROOT.RooArgList(lumi,cross_sig,eff,B_R_meson,B_R_))
#Nsig = ROOT.RooFormulaVar("Nsig","@0*(@1*@2+@3*@4)*@5*@6",ROOT.RooArgList(lumi,cross_sig_ggH,eff_ggH,cross_sig_VBF,eff_VBF,B_R_meson,B_R_)) #this is if you want to combine the two Higgs production mechanisms

#ADD PDFs --------------------------------------------------------------------------------------------------
argListN = ROOT.RooArgList(Nsig,Nbkg)
print "arg list N size = ",argListN.getSize()
argListPDFs_chebychev   = ROOT.RooArgList(signalPDF1,bkgPDF_chebychev)
argListPDFs_exponential = ROOT.RooArgList(signalPDF2,bkgPDF_exponential)

print "arg lists created!"
totPDF_chebychev   = ROOT.RooAddPdf("totPDF_chebychev","The total PDF",argListPDFs_chebychev,argListN)
totPDF_exponential = ROOT.RooAddPdf("totPDF_exponential","The total PDF",argListPDFs_exponential,argListN)

print "PDFs added!"

#Generate dataset for blind analysis -----------------------------------------------------------------------------------
datasetGenerated_chebychev   = totPDF_chebychev.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_exponential = totPDF_exponential.generate(ROOT.RooArgSet(mass),Ndata)

print "Dataset generated!"

#FIT
totPDF_chebychev.fitTo(datasetGenerated_chebychev)
totPDF_exponential.fitTo(datasetGenerated_exponential)

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
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_"+SAMPLE+"_chebychev.pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_"+SAMPLE+"_chebychev.png")

xframe_exponential = mass.frame(28)
xframe_exponential.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_exponential.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
datasetGenerated_exponential.plotOn(xframe_exponential)
totPDF_exponential.plotOn(xframe_exponential)

c2 = ROOT.TCanvas()
c2.cd()
c2.SetTitle("")
xframe_exponential.Draw()
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_"+SAMPLE+"_exponential.pdf")
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_to_generated_dataset_"+SAMPLE+"_exponential.png")

#Final useful information for debugging
print ""
print "######################################################"
print "Ndata          = ", Ndata
print "totalSigEvents = ",totalSigEvents
print "NsigPassed     = ",NsigPassed
print "Signal efficiency after tightselection = ", sigEfficiency
print "######################################################"
print ""

#create Workspace
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totPDF_chebychev)
#getattr(workspace,'import')(datasetGenerated_chebychev)
getattr(workspace,'import')(totPDF_exponential)
#getattr(workspace,'import')(datasetGenerated_exponential)

fOut = ROOT.TFile("workspaces/ws_data_blind.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()


if doBiasStudy:
	genPDF = totPDF_exponential
	fitPDF = totPDF_chebychev

	mcstudy = ROOT.RooMCStudy(genPDF, ROOT.RooArgSet(mass), ROOT.RooFit.Silence(),ROOT.RooFit.FitModel(genPDF), ROOT.RooFit.Extended(), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))

	#mcstudy = ROOT.RooMCStudy(totPDF, ROOT.RooArgSet(mass), ROOT.RooFit.Silence(), ROOT.RooFit.Extended(), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))

	#mcstudy = ROOT.RooMCStudy(totPDF, ROOT.RooArgSet(mass), ROOT.RooFit.Silence(), ROOT.RooFit.ProtoData(datasetGenerated), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))
	print "mcstudy set"
	mcstudy.generateAndFit(1000)
	print "toy generated"

	#Plot the distributions of the fitted parameter, the error and the pull
	BRval_frame = mcstudy.plotParam(B_R_, ROOT.RooFit.Bins(28))
	print "BRval_frame set!"
	BRerr_frame = mcstudy.plotError(B_R_, ROOT.RooFit.Bins(28))
	print "BRerr_frame set!"
	BRpull_frame = mcstudy.plotPull(B_R_, ROOT.RooFit.Bins(28), ROOT.RooFit.FitGauss(1))
	print "fit done"

	#print "RMS=", BRpull_frame.getHist().GetRMS()
	#print "mean=", BRpull_frame.getHist().GetMean()

	#Plot distribution of minimized likelihood
	NLLframe = mcstudy.plotNLL(ROOT.RooFit.Bins(28))

	BRval_frame.SetTitle("")
	BRval_frame.SetTitleOffset(1.5,"y")
	BRval_frame.SetXTitle("BR(H#rightarrow meson,#gamma)")
	BRerr_frame.SetTitle("")
	BRerr_frame.SetTitleOffset(1.5,"y")
	BRerr_frame.SetXTitle("#sigma_{BR(H#rightarrow meson,#gamma)}")
	BRpull_frame.SetTitle("")
	#BRpull_frame.SetMaximum(900)
	BRpull_frame.SetTitleOffset(1.5,"y")
	BRpull_frame.SetXTitle("PULL_{BR(H#rightarrow meson,#gamma)}")
	NLLframe.SetTitle("")

	#Actually plot
	c3 = ROOT.TCanvas()
	c3.cd()
	BRval_frame.Draw()
	c3.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRval.pdf")
	c3.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRval.png")

	c4 = ROOT.TCanvas()
	c4.cd()
	BRerr_frame.Draw()
	c4.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRerr.pdf")
	c4.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRerr.png")

	c5 = ROOT.TCanvas()
	c5.cd()
	BRpull_frame.Draw()
	c5.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull.pdf")
	c5.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull.png")

	c6 = ROOT.TCanvas()
	c6.cd()
	NLLframe.Draw()
	c6.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/NLL.pdf")
	c6.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/NLL.png")

