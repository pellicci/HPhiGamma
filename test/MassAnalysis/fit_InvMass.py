import ROOT


mass = ROOT.RooRealVar("_HiggsMass","_HiggsMass",100.,150.,"GeV/c^2")
mass.setRange("LowSideband",100.,120.)
mass.setRange("HighSideband",130.,150.)

a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",0.1,-10.,10.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",0.,-10.,10.)
c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.,-10.,10.)
d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",0.,-10.,10.)
e_bkg = ROOT.RooRealVar("e_bkg","e_bkg",0.,-10.,10.)

bkgPDF = ROOT.RooChebychev("bkgPDF","bkgPDF",mass,ROOT.RooArgList(a_bkg,b_bkg,c_bkg,d_bkg,e_bkg))

fileInput = ROOT.TFile("histos/latest_production/histos_Data.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))

bkgPDF.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"))

data_blinded = dataset.reduce("_HiggsMass < 120. || _HiggsMass > 130.")
xframe = mass.frame(30)
data_blinded.plotOn(xframe)
bkgPDF.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"))

c1 = ROOT.TCanvas()
c1.cd()
xframe.Draw()
c1.SaveAs("fitdata.pdf")