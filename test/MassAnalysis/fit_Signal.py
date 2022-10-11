import ROOT
import tdrstyle, CMS_lumi
import argparse

#bool initialization
isRhoGammaAnalysis = False
isPhiGammaAnalysis = False

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

#import the Doube Crystal Ball PDF
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+")

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

#Define the observable
mass = ROOT.RooRealVar("mesonGammaMass","mesonGammaMass",125.,100.,170.,"GeV/c^2")

#Double Crystal Ball definition
if SAMPLE == "ggH":
	dCB_pole  = ROOT.RooRealVar("dCB_pole_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB pole", 125.,120.,130.)
	dCB_width = ROOT.RooRealVar("dCB_width_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB width",1.,0.5,10.)
	dCB_aL    = ROOT.RooRealVar("dCB_aL_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB alpha left", 3., 0.1, 50.)
	dCB_aR    = ROOT.RooRealVar("dCB_aR_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB alpha right", 1., 0.1, 50.)
	dCB_nL    = ROOT.RooRealVar("dCB_nL_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB n left", 3., 0.1, 50.)
	dCB_nR    = ROOT.RooRealVar("dCB_nR_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB n right", 1., 0.1, 50.)

else:
	dCB_pole  = ROOT.RooRealVar("dCB_pole_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB pole", 125.,120.,130.)
	dCB_width = ROOT.RooRealVar("dCB_width_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB width",1.,0.5,10.)
	dCB_aL    = ROOT.RooRealVar("dCB_aL_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB alpha left", 1., 0.1, 50.)
	dCB_aR    = ROOT.RooRealVar("dCB_aR_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB alpha right", 1., 0.1, 50.)
	dCB_nL    = ROOT.RooRealVar("dCB_nL_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB n left", 3., 0.1, 50.)
	dCB_nR    = ROOT.RooRealVar("dCB_nR_"+CHANNEL+"_GFcat_"+SAMPLE, "Double CB n right", 1., 0.1, 50.)

signalPDF = ROOT.RooDoubleCBFast("crystal_ball_"+CHANNEL+"_GFcat_"+SAMPLE, "Double Crystal Ball", mass, dCB_pole, dCB_width, dCB_aL, dCB_nL, dCB_aR, dCB_nR)

#Input file and tree
fileInput.cd()
tree = fileInput.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also
dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))

#Do the fit
fitResult = signalPDF.fitTo(dataset)

#Plot
xframe = mass.frame(40)

#nParam   = fitResult.floatParsFinal().getSize()
#chi2     = xframe.chiSquare(nParam)#Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
#cut_chi2 = "{:.2f}".format(chi2) #Crop the chi2 to 2 decimal digits

dataset.plotOn(xframe)
signalPDF.plotOn(xframe)
xframe.SetTitle("")
ROOT.gStyle.SetOptFit(1111)
xframe.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe.GetXaxis().SetRangeUser(100.,170.)
xframe.GetYaxis().SetMaxDigits(2)

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe.Draw()

c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fitsignal_"+CHANNEL+"_GFcat_"+SAMPLE+".pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fitsignal_"+CHANNEL+"_GFcat_"+SAMPLE+".png")                                                                                                                                                          


#create Workspace
workspace = ROOT.RooWorkspace("workspace_STAT_"+CHANNEL+"_GFcat_"+SAMPLE+"_2018")
getattr(workspace,'import')(signalPDF)

fOut = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_2018.root","RECREATE")
fOut.cd()
workspace.Write()
xframe.Write("xframe")
fOut.Close()	
