import ROOT

mass = ROOT.RooRealVar("mass_KKg","mass_KKg",100.,170.,"GeV/c^2")
mass.setRange("LowSideband",100.,120.)
mass.setRange("HighSideband",130.,170.)

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
nEntries = dataset.numEntries() 

bkgPDF.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"))

data_blinded = dataset.reduce("mass_KKg < 120. || mass_KKg > 130.")
xframe = mass.frame(50)
data_blinded.plotOn(xframe)
bkgPDF.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"))

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe.Draw()
c1.SaveAs("fitdata.pdf")

#create Workspace
print "**************************************n. events = ",nEntries
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(bkgPDF)

fOut = ROOT.TFile("workspaces/ws_data.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()
