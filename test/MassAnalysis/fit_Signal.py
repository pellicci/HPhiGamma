import ROOT

#BOOLs for signal PDF
isCBplusGauss = False 
isDoubleCB = True

#import the Doube Crystal Ball PDF
if isDoubleCB:
#	ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+")
	ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/DCB.C+")


#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#Define the observable
mass = ROOT.RooRealVar("mass_KKg","mass_KKg",125.,100.,150.,"GeV/c^2")

#Define the signal lineshape
if isCBplusGauss:
	Gauss_pole = ROOT.RooRealVar("Gauss_pole","The gaussian pole", 125.,120.,130.)
	Gauss_sigma = ROOT.RooRealVar("Gauss_sigma","The gaussian sigma",5.,0.01,10.)
	Gauss_pdf = ROOT.RooGaussian("Gauss_pdf","The Gaussian",mass,Gauss_pole,Gauss_sigma)

	fracSig_prime = ROOT.RooRealVar("fracSig_prime","Partial fraction",0.5,0.,1.)
	fracSig = ROOT.RooRealVar("fracSig","Fraction",0.5,0.,1.)

	#Crystal Ball PDF parameters
	mean = ROOT.RooRealVar("mean","mean",125.,120.,130.)
	width = ROOT.RooRealVar("width","width",2.5,0.1,5.)
	sigma = ROOT.RooRealVar("sigma","sigma",2.,0.1,5.)
	alpha = ROOT.RooRealVar("alpha","alpha",1.,0.1,10.)
	enne = ROOT.RooRealVar("enne","enne",2.,0.1,20.)

	#Initialize a Crystal Ball pdf
	CBpdf = ROOT.RooCBShape("signalPDF","signalPDF",mass,mean,sigma,alpha,enne)

	#Add the Gaussian and the CB
	signalPDF = ROOT.RooAddPdf("totSignal","Total signal PDF",ROOT.RooArgList(CBpdf,Gauss_pdf),ROOT.RooArgList(fracSig))

#Double Crystal Ball definition
if isDoubleCB:
	dCB_pole  = ROOT.RooRealVar("dCB_pole", "Double CB pole", 125.,100.,150.)
	dCB_width = ROOT.RooRealVar("dCB_width", "Double CB width",1.,0.01,10.)
	dCB_aL    = ROOT.RooRealVar("dCB_aL", "Double CB alpha left", 3., 0.1, 50.)
	dCB_aR    = ROOT.RooRealVar("dCB_aR", "Double CB alpha right", 1., 0.1, 50.)
	dCB_nL    = ROOT.RooRealVar("dCB_nL", "Double CB n left", 3., 0.1, 50.)
	dCB_nR    = ROOT.RooRealVar("dCB_nR", "Double CB n right", 1., 0.1, 50.)
	signalPDF = ROOT.RooDoubleCBFast("dCB", "Double Crystal Ball", mass, dCB_pole, dCB_width, dCB_aL, dCB_nL, dCB_aR, dCB_nR)


#Input file and tree
fileInput = ROOT.TFile("histos/latest_production/histos_SR_Signal.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also
dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))

#Do the fit
signalPDF.fitTo(dataset)

#Plot
xframe = mass.frame(100)
dataset.plotOn(xframe)
signalPDF.plotOn(xframe)
xframe.SetTitle("")
xframe.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
xframe.GetXaxis().SetRangeUser(100.,150.)

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe.Draw()
c1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fitsignal.pdf")
c1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fitsignal.png")

#create Workspace
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(signalPDF)

fOut = ROOT.TFile("workspaces/ws_signal.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()	
