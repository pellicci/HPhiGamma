import ROOT

isPhiGammaAnalysis = False

#import the Doube Crystal Ball PDF
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+")

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#Define the observable
mass = ROOT.RooRealVar("mesonGammaMass","mesonGammaMass",125.,100.,150.,"GeV/c^2")

#Double Crystal Ball definition
dCB_pole  = ROOT.RooRealVar("dCB_pole", "Double CB pole", 125.,120.,130.)
dCB_width = ROOT.RooRealVar("dCB_width", "Double CB width",1.,0.01,10.)
dCB_aL    = ROOT.RooRealVar("dCB_aL", "Double CB alpha left", 3., 0.1, 50.)
dCB_aR    = ROOT.RooRealVar("dCB_aR", "Double CB alpha right", 1., 0.1, 50.)
dCB_nL    = ROOT.RooRealVar("dCB_nL", "Double CB n left", 3., 0.1, 50.)
dCB_nR    = ROOT.RooRealVar("dCB_nR", "Double CB n right", 1., 0.1, 50.)

signalPDF = ROOT.RooDoubleCBFast("signalPDF", "Double Crystal Ball", mass, dCB_pole, dCB_width, dCB_aL, dCB_nL, dCB_aR, dCB_nR)

#Input file and tree
fileInput = ROOT.TFile("histos/latest_production/histos_SR_Signal.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also
dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))

#Do the fit
signalPDF.fitTo(dataset)

#Plot
xframe = mass.frame(40)
dataset.plotOn(xframe)
signalPDF.plotOn(xframe)
xframe.SetTitle("")
xframe.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
xframe.GetXaxis().SetRangeUser(100.,150.)

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe.Draw()

if isPhiGammaAnalysis:
	c1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fitsignal.pdf")
	c1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fitsignal.png")
else:
	c1.SaveAs("~/cernbox/www/MyAnalysis/HRhoGamma/MassAnalysis/latest_production/fitsignal.pdf")
	c1.SaveAs("~/cernbox/www/MyAnalysis/HRhoGamma/MassAnalysis/latest_production/fitsignal.png")

#create Workspace
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(signalPDF)

fOut = ROOT.TFile("workspaces/ws_signal.root","RECREATE")
fOut.cd()
workspace.Write()
xframe.Write("xframe")
fOut.Close()	