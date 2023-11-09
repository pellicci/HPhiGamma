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
p.add_argument('CAT_option', help='Type the category') #flag for bkg estimation
args = p.parse_args()

if args.Decay_channel_option == "Phi":
	isPhiGammaAnalysis = True
	CHANNEL = "Phi"
	print "H -> PhiGamma analysis"
if args.Decay_channel_option == "Rho":
    isRhoGammaAnalysis = True
    CHANNEL = "Rho"
    print "H -> RhoGamma analysis"
if args.Decay_channel_option == "K0s":
    isK0sGammaAnalysis = True
    CHANNEL = "K0s"
    print "H -> K0sGamma analysis"

CATEGORY = args.CAT_option


if CATEGORY == "bdt0": fileInput_ggH = ROOT.TFile("histos/latest_production/histos_SR_BDTcat0_SignalggH.root")
if CATEGORY == "bdt1": fileInput_ggH = ROOT.TFile("histos/latest_production/histos_SR_BDTcat1_SignalggH.root")


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
mass = ROOT.RooRealVar("mesonGammaMass","mesonGammaMass",110.,160.,"GeV")

#Data input -------------------------------------------------------------
if CATEGORY == "bdt0": fileInputData = ROOT.TFile("histos/latest_production/histos_SR_BDTcat0_Data.root")
if CATEGORY == "bdt1": fileInputData = ROOT.TFile("histos/latest_production/histos_SR_BDTcat1_Data.root")
fileInputData.cd()
mytreeData = fileInputData.Get("tree_output")

#Signal input -------------------------------------------------------------
if CATEGORY == "bdt0": fileInputSignal_ggH = ROOT.TFile("histos/latest_production/histos_SR_BDTcat0_SignalggH.root")
if CATEGORY == "bdt1": fileInputSignal_ggH = ROOT.TFile("histos/latest_production/histos_SR_BDTcat1_SignalggH.root")
fileInputSignal_ggH.cd()
mytreeSignal_ggH = fileInputSignal_ggH.Get("tree_output")


#Import from workspaces --------------------------------------------------------------------
fInputSignalPDF = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_"+CATEGORY+"_2018.root") #there's only one ws for both ggH and VBF 
fInputSignalPDF.cd()
workspaceSignal = fInputSignalPDF.Get("workspace_STAT_"+CHANNEL+"_GFcat_"+CATEGORY+"_2018")

print "################################ WORKSPACE SIGNAL #############################################################"
workspaceSignal.Print() #print some info of the workspace content
print "###############################################################################################################"
print ""

#Signal PDF
signalPDF_ggH = workspaceSignal.pdf("crystal_ball_"+CHANNEL+"_GFcat_"+CATEGORY+"_ggH") #SIGNAL PDF from the workspace
workspaceSignal.var("dCB_nL_"+CHANNEL+"_GFcat_"+CATEGORY+"_ggH").setConstant(1) 
workspaceSignal.var("dCB_nR_"+CHANNEL+"_GFcat_"+CATEGORY+"_ggH").setConstant(1)
workspaceSignal.var("dCB_aL_"+CHANNEL+"_GFcat_"+CATEGORY+"_ggH").setConstant(1)
workspaceSignal.var("dCB_aR_"+CHANNEL+"_GFcat_"+CATEGORY+"_ggH").setConstant(1)
workspaceSignal.var("dCB_pole_"+CHANNEL+"_GFcat_"+CATEGORY+"_ggH").setConstant(1)
workspaceSignal.var("dCB_width_"+CHANNEL+"_GFcat_"+CATEGORY+"_ggH").setConstant(1)

#Background PDF
fInputSidebandsPDF = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_"+CATEGORY+"_2018.root")
fInputSidebandsPDF.cd()
workspaceSidebands   = fInputSidebandsPDF.Get("workspace_STAT_"+CHANNEL+"_GFcat_"+CATEGORY+"_2018")
bkgPDF_chebychev     = workspaceSidebands.pdf("chebychev_GFcat_"+CATEGORY+"_bkg") #BKG PDF from the workspace
bkgPDF_bernstein     = workspaceSidebands.pdf("bernstein_GFcat_"+CATEGORY+"_bkg") #BKG PDF from the workspace
#bkgPDF_chebychev4    = workspaceSidebands.pdf("chebychev4_GFcat_"+CATEGORY+"_bkg") #BKG PDF from the workspace
bkgPDF_exponential   = workspaceSidebands.pdf("exponential_GFcat_"+CATEGORY+"_bkg") #BKG PDF from the workspace

print "############################### WORKSPACE BKG ##############################################################"
workspaceSidebands.Print()
print "############################################################################################################"
print ""

#Define the POI -------------------------------------------------------------------------------------
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",0.05,-.1,.1) #parameter of interest
#B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",0.) #parameter of interest

#Number of events in m_KKg
fInput_ggH = ROOT.TFile("rootfiles/latest_production/MC/signals/HPhiGammaAnalysis_Signal_"+CHANNEL+"_ggH.root")
mytree_ggH         = fInput_ggH.Get("HPhiGammaAnalysis/mytree")
h_Events_ggH       = fInput_ggH.Get("HPhiGammaAnalysis/h_Events")
NsigPassed_ggH     = mytreeSignal_ggH.GetEntriesFast()
totalSigEvents_ggH = h_Events_ggH.GetBinContent(1)
sigEfficiency_ggH  = float(NsigPassed_ggH/totalSigEvents_ggH)

Ndata = mytreeData.GetEntriesFast()
#Ndata = 1000
#N bkg estimation
Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata, 0.5*Ndata, 3*Ndata)

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
cross_sig_ggH = ROOT.RooRealVar("cross_sig","The signal cross section",48.85)#46.87pb (gluon fusion)
lumi = ROOT.RooRealVar("lumi","The luminosity",39540)#pb^-1
eff_ggH  = ROOT.RooRealVar("eff_ggH","efficiency ggH",sigEfficiency_ggH) # = (n sig over tightselection)/49750
if CHANNEL == "Phi": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",0.489) # = BR(phi -> K+K-)
if CHANNEL == "Rho": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",1.) # = BR(rho -> PiPi)
if CHANNEL == "K0s": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",1.) # = BR(rho -> KPi)
Nsig_ggH = ROOT.RooFormulaVar("Nsig_ggH","@0*@1*@2*@3*@4",ROOT.RooArgList(lumi,cross_sig_ggH,eff_ggH,B_R_meson,B_R_))

#ADD PDFs --------------------------------------------------------------------------------------------------
argListN = ROOT.RooArgList(Nsig_ggH,Nbkg)
print "arg list N size = ",argListN.getSize()
argListPDFs_chebychev   = ROOT.RooArgList(signalPDF_ggH,bkgPDF_chebychev)
argListPDFs_bernstein   = ROOT.RooArgList(signalPDF_ggH,bkgPDF_bernstein)
argListPDFs_exponential = ROOT.RooArgList(signalPDF_ggH,bkgPDF_exponential)
#argListPDFs_chebychev4  = ROOT.RooArgList(signalPDF_ggH,bkgPDF_chebychev4)

print "arg lists created!"
totPDF_chebychev   = ROOT.RooAddPdf("totPDF_chebychev","Cheby",argListPDFs_chebychev,argListN)
totPDF_bernstein   = ROOT.RooAddPdf("totPDF_bernstein","Bern",argListPDFs_bernstein,argListN)
totPDF_exponential = ROOT.RooAddPdf("totPDF_exponential","Exp",argListPDFs_exponential,argListN)
#totPDF_chebychev4  = ROOT.RooAddPdf("totPDF_chebychev4","Cheby5",argListPDFs_chebychev4,argListN)

print "PDFs added!"

#Generate dataset for blind analysis -----------------------------------------------------------------------------------
datasetGenerated_chebychev   = totPDF_chebychev.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_bernstein   = totPDF_bernstein.generate(ROOT.RooArgSet(mass),Ndata)
#datasetGenerated_chebychev4  = totPDF_chebychev4.generate(ROOT.RooArgSet(mass),Ndata)
datasetGenerated_exponential = totPDF_exponential.generate(ROOT.RooArgSet(mass),Ndata)

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

genPDF = totPDF_bernstein
fitPDF = totPDF_chebychev

for fitPDF in [totPDF_chebychev]:
    for genPDF in [totPDF_bernstein,totPDF_exponential]:  #ROOT.RooFit.Silence()
        mcstudy = ROOT.RooMCStudy(genPDF, ROOT.RooArgSet(mass),ROOT.RooFit.FitModel(fitPDF), ROOT.RooFit.Extended(1), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))
        mcstudy.generateAndFit(500)

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

        BRpull_frame = mcstudy.plotPull(B_R_, ROOT.RooFit.Bins(100), ROOT.RooFit.FitGauss(1))
        BRpull_frame.SetTitle("")
        BRpull_frame.SetTitleOffset(1.5,"y")
        BRpull_frame.SetXTitle("Gen: "+genPDF.GetTitle()+" Vs Fit: "+fitPDF.GetTitle())
        BRpull_frame.SetMaximum(1.3*BRpull_frame.GetMaximum())

        BRpull_frame.Draw()
        leg2.Draw()


        multicanvas.Draw()

        multicanvas.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull_"+genPDF.GetTitle()+"Vs"+fitPDF.GetTitle()+"_"+CATEGORY+".png")
        multicanvas.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull_"+genPDF.GetTitle()+"Vs"+fitPDF.GetTitle()+"_"+CATEGORY+".pdf")

