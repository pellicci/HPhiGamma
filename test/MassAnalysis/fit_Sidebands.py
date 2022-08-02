import ROOT

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

isPhiGammaAnalysis = False

#Parameters of the PDF
mass = ROOT.RooRealVar("mesonGammaMass","mesonGammaMass",100.,170.,"GeV")
mass.setRange("LowSideband",100.,115.)
mass.setRange("HighSideband",135.,170.)

a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",0.01,-1.,1.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",-0.2,-100.,100.)
#c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.,-10.,10.)
#d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",0.,-10.,10.)

#Initialize a Chebychev pdf
bkgPDF = ROOT.RooChebychev("bkgPDF","bkgPDF",mass,ROOT.RooArgList(a_bkg,b_bkg))

#Input file and tree
fileInput = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also
dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))
nEntries = dataset.numEntries() 

#Do the fit
bkgPDF.fitTo(dataset,ROOT.RooFit.Range("LowSideband,HighSideband"))

#Give the blind range
data_blinded = dataset.reduce("mesonGammaMass < 115. || mesonGammaMass > 135.")

#Plot
xframe = mass.frame(40)
data_blinded.plotOn(xframe)
bkgPDF.plotOn(xframe,ROOT.RooFit.Range("LowSideband,HighSideband"))
xframe.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe.Draw()

c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_sidebands.pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fit_sidebands.png")
	


#create Workspace
print "**************************************n. events = ",nEntries
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(bkgPDF)

fOut = ROOT.TFile("workspaces/ws_sidebands.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()
