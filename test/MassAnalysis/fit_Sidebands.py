import ROOT
import tdrstyle, CMS_lumi
import argparse

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#bool initialization
isRhoGammaAnalysis = False
isPhiGammaAnalysis = False

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


#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

#Parameters of the PDF
mass = ROOT.RooRealVar("mesonGammaMass","mesonGammaMass",100.,170.,"GeV")
mass.setRange("LowSideband",100.,115.)
mass.setRange("HighSideband",135.,170.)

#Initialize a Chebychev pdf
a_bkg = ROOT.RooRealVar("a_bkg_chebychev_"+CHANNEL+"_GFcat","a_bkg",-1.08,-1.1,-1.05)
b_bkg = ROOT.RooRealVar("b_bkg_chebychev_"+CHANNEL+"_GFcat","b_bkg",0.4,-1.,1.)
c_bkg = ROOT.RooRealVar("c_bkg_chebychev_"+CHANNEL+"_GFcat","c_bkg",0.01,-0.1,0.1)
d_bkg = ROOT.RooRealVar("d_bkg_chebychev_"+CHANNEL+"_GFcat","d_bkg",0.,-1.,1.)
bkgPDF_chebychev = ROOT.RooChebychev("chebychev_GFcat_bkg","bkgPDF",mass,ROOT.RooArgList(a_bkg,b_bkg))

#Initialize a exponential pdf
e_bkg = ROOT.RooRealVar("e_bkg_exponential_"+CHANNEL+"_GFcat","e_bkg",-10.,10.)
bkgPDF_exponential = ROOT.RooExponential("exponential_GFcat_bkg","bkgPDF",mass,e_bkg)

#Initialize a Bernstein pdf
#if PDFindex == 2:
#	a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",0.,100.)
#	b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",-2.,100.)
#	c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",-1.,-2.,2.)
#	d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",0.,-1.,1.)
#	bkgPDF = ROOT.RooBernstein("bkgPDF","bkgPDF",mass,ROOT.RooArgList(a_bkg,b_bkg,c_bkg))


#Input file and tree
fileInput = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also
dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))
nEntries = dataset.numEntries() 

#Do the fit
fitResult_chebychev   = bkgPDF_chebychev.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Save())
fitResult_exponential = bkgPDF_exponential.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Save())

print "minNll = ", fitResult_chebychev.minNll()
print "2Delta_minNll = ", 2*(1226.50106409-fitResult_chebychev.minNll()) # If 2*(NLL(N)-NLL(N+1)) > 3.85 -> N+1 is significant improvement

#Give the blind range
data_blinded = dataset.reduce("mesonGammaMass < 115. || mesonGammaMass > 135.")

#Plot
c1 = ROOT.TCanvas()
c1.cd()
#c1.SetTitle("")

xframe = mass.frame(28)
data_blinded.plotOn(xframe)
bkgPDF_chebychev.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_chebychev"),ROOT.RooFit.LineColor(ROOT.kBlue))
bkgPDF_exponential.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_exponential"), ROOT.RooFit.LineColor(ROOT.kRed))
xframe.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe.SetMaximum(1.3*xframe.GetMaximum())

xframe.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

leg1 = ROOT.TLegend(0.65,0.62,0.87,0.90) #right positioning
leg1.SetHeader(" ")
leg1.SetNColumns(1)
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.AddEntry("bkgPDF_chebychev","Chebychev pdf","l")
leg1.AddEntry("bkgPDF_exponential","Exponential pdf","l")

CMS_lumi.CMS_lumi(c1, iPeriod, iPos) #Print integrated lumi and energy information

c1.Update()
leg1.Draw()

c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_sidebands.pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_sidebands.png")
	


#create Workspace
print "************************************** n. events = ",nEntries
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(bkgPDF_chebychev)
getattr(workspace,'import')(bkgPDF_exponential)

fOut = ROOT.TFile("workspaces/ws_sidebands.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()
