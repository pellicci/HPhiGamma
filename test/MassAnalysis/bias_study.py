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


fileInput_ggH = ROOT.TFile("histos/latest_production/histos_SR_BDTcat0_SignalggH.root")


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
fileInputData = ROOT.TFile("histos/latest_production/histos_SR_BDTcat0_Data.root")
fileInputData.cd()
mytreeData = fileInputData.Get("tree_output")

#Signal input -------------------------------------------------------------
fileInputSignal_ggH = ROOT.TFile("histos/latest_production/histos_SR_BDTcat0_SignalggH.root")
fileInputSignal_ggH.cd()
mytreeSignal_ggH = fileInputSignal_ggH.Get("tree_output")


#Import from workspaces --------------------------------------------------------------------
fInputSignalPDF = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_2018.root") #there's only one ws for both ggH and VBF 
fInputSignalPDF.cd()
workspaceSignal = fInputSignalPDF.Get("workspace_STAT_"+CHANNEL+"_GFcat_2018")

print "################################ WORKSPACE SIGNAL #############################################################"
workspaceSignal.Print() #print some info of the workspace content
print "###############################################################################################################"
print ""

#Signal PDF
signalPDF_ggH = workspaceSignal.pdf("crystal_ball_"+CHANNEL+"_GFcat_ggH") #SIGNAL PDF from the workspace
workspaceSignal.var("dCB_nL_"+CHANNEL+"_GFcat_ggH").setConstant(1) 
workspaceSignal.var("dCB_nR_"+CHANNEL+"_GFcat_ggH").setConstant(1)
workspaceSignal.var("dCB_aL_"+CHANNEL+"_GFcat_ggH").setConstant(1)
workspaceSignal.var("dCB_aR_"+CHANNEL+"_GFcat_ggH").setConstant(1)
workspaceSignal.var("dCB_pole_"+CHANNEL+"_GFcat_ggH").setConstant(1)
workspaceSignal.var("dCB_width_"+CHANNEL+"_GFcat_ggH").setConstant(1)

#Background PDF
fInputSidebandsPDF = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_2018.root")
fInputSidebandsPDF.cd()
workspaceSidebands   = fInputSidebandsPDF.Get("workspace_STAT_"+CHANNEL+"_GFcat_2018")
bkgPDF_chebychev     = workspaceSidebands.pdf("chebychev_GFcat_bkg") #BKG PDF from the workspace
bkgPDF_exponential   = workspaceSidebands.pdf("exponential_GFcat_bkg") #BKG PDF from the workspace
bkgPDF_bernstein     = workspaceSidebands.pdf("bernstein_GFcat_bkg") #BKG PDF from the workspace
bkgPDF_chebychev4    = workspaceSidebands.pdf("chebychev4_GFcat_bkg") #BKG PDF from the workspace

print "############################### WORKSPACE BKG ##############################################################"
workspaceSidebands.Print()
print "############################################################################################################"
print ""

#Define the POI -------------------------------------------------------------------------------------
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",0.001,-0.1,0.1) #parameter of interest
#B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",0.) #parameter of interest

#Number of events in m_KKg
fInput_ggH = ROOT.TFile("rootfiles/latest_production/MC/signals/HPhiGammaAnalysis_Signal_"+CHANNEL+"_ggH.root")
mytree_ggH         = fInput_ggH.Get("HPhiGammaAnalysis/mytree")
h_Events_ggH       = fInput_ggH.Get("HPhiGammaAnalysis/h_Events")
NsigPassed_ggH     = mytreeSignal_ggH.GetEntriesFast()
totalSigEvents_ggH = h_Events_ggH.GetBinContent(1)
sigEfficiency_ggH  = float(NsigPassed_ggH/totalSigEvents_ggH)

Ndata = mytreeData.GetEntriesFast()

#N bkg estimation
Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata, 0.5*Ndata, 2*Ndata)

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
cross_sig_ggH = ROOT.RooRealVar("cross_sig","The signal cross section",48.85)#46.87pb (gluon fusion)
lumi = ROOT.RooRealVar("lumi","The luminosity",39540)#pb^-1
eff_ggH  = ROOT.RooRealVar("eff_ggH","efficiency ggH",sigEfficiency_ggH) # = (n sig over tightselection)/49750
if CHANNEL == "Phi": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",0.489) # = BR(phi -> K+K-)
if CHANNEL == "Rho": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",1.) # = BR(rho -> PiPi)
Nsig_ggH = ROOT.RooFormulaVar("Nsig_ggH","@0*@1*@2*@3*@4",ROOT.RooArgList(lumi,cross_sig_ggH,eff_ggH,B_R_meson,B_R_))

#ADD PDFs --------------------------------------------------------------------------------------------------
argListN = ROOT.RooArgList(Nsig_ggH,Nbkg)
print "arg list N size = ",argListN.getSize()
argListPDFs_chebychev   = ROOT.RooArgList(signalPDF_ggH,bkgPDF_chebychev)
argListPDFs_exponential = ROOT.RooArgList(signalPDF_ggH,bkgPDF_exponential)
argListPDFs_bernstein   = ROOT.RooArgList(signalPDF_ggH,bkgPDF_bernstein)
argListPDFs_chebychev4  = ROOT.RooArgList(signalPDF_ggH,bkgPDF_chebychev4)

print "arg lists created!"
totPDF_chebychev   = ROOT.RooAddPdf("totPDF_chebychev","Cheby",argListPDFs_chebychev,argListN)
totPDF_exponential = ROOT.RooAddPdf("totPDF_exponential","Exp",argListPDFs_exponential,argListN)
totPDF_bernstein   = ROOT.RooAddPdf("totPDF_bernstein","Bern",argListPDFs_bernstein,argListN)
totPDF_chebychev4  = ROOT.RooAddPdf("totPDF_chebychev4","Cheby4",argListPDFs_chebychev,argListN)

print "PDFs added!"

#Generate dataset for blind analysis -----------------------------------------------------------------------------------
datasetGenerated_chebychev   = totPDF_chebychev.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_exponential = totPDF_exponential.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_bernstein   = totPDF_bernstein.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_chebychev4  = totPDF_chebychev4.generate(ROOT.RooArgSet(mass),Ndata)

print "Dataset generated!"

#Final useful information for debugging
print ""
print "######################################################"
print "Ndata              = ",Ndata
print "totalSigEvents ggH = ",totalSigEvents_ggH
print "NsigPassed ggH     = ",NsigPassed_ggH
print "######################################################"
print ""

multicanvas = ROOT.TCanvas()
multicanvas.cd()

genPDF = totPDF_chebychev
fitPDF = totPDF_exponential

mcstudy = ROOT.RooMCStudy(genPDF, ROOT.RooArgSet(mass), ROOT.RooFit.Silence(),ROOT.RooFit.FitModel(fitPDF), ROOT.RooFit.Extended(), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))
mcstudy.generateAndFit(5000)

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

BRpull_frame = mcstudy.plotPull(B_R_, ROOT.RooFit.Bins(56), ROOT.RooFit.FitGauss(1))
BRpull_frame.SetTitle("")
BRpull_frame.SetTitleOffset(1.5,"y")
BRpull_frame.SetXTitle("Gen: "+genPDF.GetTitle()+" Vs Fit: "+fitPDF.GetTitle())
BRpull_frame.SetMaximum(1.1*BRpull_frame.GetMaximum())

BRpull_frame.Draw()
leg2.Draw()


multicanvas.Draw()

multicanvas.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull_"+genPDF.GetTitle()+"Vs"+fitPDF.GetTitle()+".png")
multicanvas.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull_"+genPDF.GetTitle()+"Vs"+fitPDF.GetTitle()+".pdf")

