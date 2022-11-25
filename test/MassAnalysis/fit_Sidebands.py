import ROOT
import tdrstyle, CMS_lumi
import argparse

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#bool initialization
isRhoGammaAnalysis = False
isPhiGammaAnalysis = False
doBiasStudy = False

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
mass.setRange("HighSideband",130.,170.)

#Initialize a Chebychev pdf
a_bkg = ROOT.RooRealVar("a_bkg_chebychev_"+CHANNEL+"_GFcat","a_bkg",-1.005,-1.01,-1.003)
b_bkg = ROOT.RooRealVar("b_bkg_chebychev_"+CHANNEL+"_GFcat","b_bkg",0.3,0.1,1.)
c_bkg = ROOT.RooRealVar("c_bkg_chebychev_"+CHANNEL+"_GFcat","c_bkg",0.,-0.1,0.)
d_bkg = ROOT.RooRealVar("d_bkg_chebychev_"+CHANNEL+"_GFcat","d_bkg",-0.05,-0.1,0.)
bkgPDF_chebychev = ROOT.RooChebychev("chebychev_GFcat_bkg","bkgPDF",mass,ROOT.RooArgList(a_bkg,b_bkg,c_bkg))

#Initialize a exponential pdf
e_bkg = ROOT.RooRealVar("e_bkg_exponential_"+CHANNEL+"_GFcat","e_bkg",-0.031,-0.4,0.2)
bkgPDF_exponential = ROOT.RooExponential("exponential_GFcat_bkg","bkgPDF",mass,e_bkg)

#Initialize a Bernstein pdf
bern_c0 = ROOT.RooRealVar('bern_c0', 'bern_c0', 9.5, 20.)
bern_c1 = ROOT.RooRealVar('bern_c1', 'bern_c1', 0.5, 1.)
bern_c2 = ROOT.RooRealVar('bern_c2', 'bern_c2', 0., 1.)
bern_c3 = ROOT.RooRealVar('bern_c3', 'bern_c3', 0.5, 0., 5.)
bern_c4 = ROOT.RooRealVar('bern_c4', 'bern_c4', 0.5, 0., 5.)
bern_c5 = ROOT.RooRealVar('bern_c5', 'bern_c5', 1e-2, 0., 0.1)
bkgPDF_bernstein = ROOT.RooBernstein("bernstein_GFcat_bkg", "bkgPDF", mass, ROOT.RooArgList(bern_c0,bern_c1,bern_c2))

#Initialize a polynomial pdf
a_pol = ROOT.RooRealVar("a_pol_"+CHANNEL+"_GFcat","a_pol",0.,-0.01,0.01)
b_pol = ROOT.RooRealVar("b_pol_"+CHANNEL+"_GFcat","b_pol",0.,-0.01,0.01)
c_pol = ROOT.RooRealVar("c_pol_"+CHANNEL+"_GFcat","c_pol",0.,-0.01,0.01)
d_pol = ROOT.RooRealVar("d_pol_"+CHANNEL+"_GFcat","d_pol",-0.05,-0.1,0.)
bkgPDF_polynomial = ROOT.RooPolynomial("polynomial_GFcat_bkg", "bkgPDF", mass, ROOT.RooArgList(a_pol,b_pol,c_pol))


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
#fitResult_exponential = bkgPDF_exponential.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Save())
fitResult_bernstein   = bkgPDF_bernstein.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Save())
fitResult_polynomial  = bkgPDF_polynomial.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Save())


print "minNll = ", fitResult_chebychev.minNll()
print "2Delta_minNll = ", 2*(13091.0575078-fitResult_chebychev.minNll()) # If 2*(NLL(N)-NLL(N+1)) > 3.85 -> N+1 is significant improvement

#Give the blind range
data_blinded = dataset.reduce("mesonGammaMass < 115. || mesonGammaMass > 130.")

#Plot
c1 = ROOT.TCanvas()
c1.cd()
#c1.SetTitle("")

xframe = mass.frame(28)
data_blinded.plotOn(xframe)
bkgPDF_chebychev.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_chebychev"),ROOT.RooFit.LineColor(ROOT.kBlue))
#bkgPDF_exponential.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_exponential"), ROOT.RooFit.LineColor(ROOT.kRed))
bkgPDF_bernstein.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_bernstein"), ROOT.RooFit.LineColor(ROOT.kGreen))
bkgPDF_polynomial.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_polynomial"), ROOT.RooFit.LineColor(ROOT.kOrange))
xframe.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe.SetMaximum(1.3*xframe.GetMaximum())

xframe.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

nParam = fitResult_chebychev.floatParsFinal().getSize()
chi2 = xframe.chiSquare(nParam)#Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2 = "{:.2f}".format(chi2) #Crop the chi2 to 2 decimal digits
print "Chi square = ",chi2
print "n param = ",nParam

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
leg1.AddEntry("bkgPDF_bernstein","Bernstein pdf","l")
leg1.AddEntry("bkgPDF_polynomial","Polynomial pdf","l")

CMS_lumi.CMS_lumi(c1, iPeriod, iPos) #Print integrated lumi and energy information

c1.Update()
leg1.Draw()

c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_sidebands.pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_sidebands.png")
	
if doBiasStudy:

	multicanvas = ROOT.TCanvas()
	multicanvas.cd()
	multicanvas.Divide(2,2)
	canvas_index = 0
	
	pdfList = [bkgPDF_chebychev,bkgPDF_bernstein]

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

			BRpull_frame = mcstudy.plotPull(mass, ROOT.RooFit.Bins(28), ROOT.RooFit.FitGauss(1))
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

#create Workspace
norm     = fileInput.Get("h_InvMass_TwoTrk_Photon").Integral() #get the normalization of ggH signal (area under ggH signal)
bkg_norm = ROOT.RooRealVar(bkgPDF_chebychev.GetName()+ "_norm", bkgPDF_chebychev.GetName()+ "_norm", norm)#,0.5*norm, 2*norm)

print "************************************** n. events = ",nEntries
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(bkgPDF_chebychev)
getattr(workspace,'import')(bkgPDF_exponential)
getattr(workspace,'import')(bkgPDF_bernstein)
getattr(workspace,'import')(bkgPDF_polynomial)
getattr(workspace,'import')(dataset)
getattr(workspace,'import')(bkg_norm)
print("integral BKG",bkg_norm.Print())
workspace.Print()

fOut = ROOT.TFile("workspaces/ws_sidebands.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()


