import ROOT

mass = ROOT.RooRealVar("mass_KKg","mass_KKg",125.,100.,150.,"GeV/c^2")
mean = ROOT.RooRealVar("mean","mean",125.,100.,150.)
width = ROOT.RooRealVar("width","width",2.5,0.1,5.)
sigma = ROOT.RooRealVar("sigma","sigma",2.,0.1,5.)
alpha = ROOT.RooRealVar("alpha","alpha",1.,0.1,10.)
enne = ROOT.RooRealVar("enne","enne",2.,0.1,20.)

signalPDF = ROOT.RooCBShape("signalPDF","signalPDF",mass,mean,sigma,alpha,enne)

fileInput = ROOT.TFile("histos/latest_production/histos_Signal.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))

signalPDF.fitTo(dataset)

xframe = mass.frame(30)
dataset.plotOn(xframe)
signalPDF.plotOn(xframe)

c1 = ROOT.TCanvas()
c1.cd()
xframe.Draw()
c1.SaveAs("fitsignal.pdf")

#create Workspace
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(signalPDF)

fOut = ROOT.TFile("workspaces/ws_signal.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()
