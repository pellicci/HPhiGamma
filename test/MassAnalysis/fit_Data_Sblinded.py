import ROOT

#INPUT
fInput1 = ROOT.TFile("workspaces/ws_signal.root")
fInput1.cd()
workspaceSignal = fInput1.Get("myworkspace")
fInput2 = ROOT.TFile("histos/latest_production/histos_Data.root")
fInput2.cd()
tree = fInput2.Get("tree_output")
fInput3 = ROOT.TFile("nEvents_Container.root")
fInput3.cd()

#VARS and CONSTANTS
workspaceSignal.var("mass_KKg").setConstant(1)
workspaceSignal.var("mean").setConstant(1)
workspaceSignal.var("alpha").setConstant(1)
workspaceSignal.var("sigma").setConstant(1)
workspaceSignal.var("enne").setConstant(1)
mytree = fInput3.Get("tree_output")
mytree.Scan()
Ndata = mytree.nEvents_BKG
signalWeight = mytree.nEvents_Signal
print "Ndata = ",Ndata
print "Nsig = ",signalWeight

Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata,Ndata-0.5*Ndata,Ndata+0.5*Ndata)
mass = ROOT.RooRealVar("mass_KKg","mass_KKg",100.,170.,"GeV")

a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",0.1,-10.,10.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",0.,-10.,10.)
c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.,-10.,10.)

#BKG PDF
bkgPDF = ROOT.RooChebychev("bkgPDF","bkgPDF",mass,ROOT.RooArgList(a_bkg,b_bkg,c_bkg))

#RIPARAMETRIZATION
cross_sig = ROOT.RooRealVar("cross_sig","The signal cross section",52.36)#52.36pb (gluon fusion + vector boson fusion xsec)
lumi = ROOT.RooRealVar("lumi","The luminosity",39540)#pb^-1
eff = ROOT.RooRealVar("eff","efficiency",0.061727) # = 1707/27654=0.061727
B_R_phiKK = ROOT.RooRealVar("B_R_phiKK","b_r_phi_KK",0.489) # = BR(phi -> K+K-)
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",10**(-5),-10**(-2),10**(-2))
Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3*@4",ROOT.RooArgList(cross_sig,lumi,eff,B_R_phiKK,B_R_))

#SIGNAL PDF
sigPDF = workspaceSignal.pdf("signalPDF")

dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))
nEntries = dataset.numEntries() 

#ADD PDFs
argListPDFs = ROOT.RooArgList(sigPDF,bkgPDF)
argListN = ROOT.RooArgList(Nsig,Nbkg)
totPDF = ROOT.RooAddPdf("totPDF","The total PDF",argListPDFs,argListN)
print "pdf added"

#FIT
totPDF.fitTo(dataset)

#PLOT
xframe = mass.frame(50)
xframe.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
dataset.plotOn(xframe)
totPDF.plotOn(xframe)

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe.Draw()
c1.SaveAs("fitdataSblinded.pdf")

#create Workspace
print "**************************************n. events = ",nEntries
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totPDF)

fOut = ROOT.TFile("workspaces/ws_data.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()
