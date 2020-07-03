import ROOT

mass = ROOT.RooRealVar("mass_KKg","mass_KKg",100.,170.,"GeV")
mass.setRange("LowSideband",100.,120.)
mass.setRange("HighSideband",130.,170.)

a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",0.1,-10.,10.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",0.,-10.,10.)
c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.,-10.,10.)

bkgPDF = ROOT.RooChebychev("bkgPDF","bkgPDF",mass,ROOT.RooArgList(a_bkg,b_bkg,c_bkg))

fileInput = ROOT.TFile("histos/latest_production/histos_Data.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))
nEntries = dataset.numEntries() 

bkgPDF.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"))

data_blinded = dataset.reduce("mass_KKg < 120. || mass_KKg > 130.")
xframe = mass.frame(70)
data_blinded.plotOn(xframe)
bkgPDF.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"))
xframe.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe.Draw()
c1.SaveAs("fitdata_blinded.pdf")

#create Workspace
print "**************************************n. events = ",nEntries
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(bkgPDF)

fOut = ROOT.TFile("workspaces/ws_data.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()
