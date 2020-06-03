import ROOT

#INPUT SETTINGS
fInput1 = ROOT.TFile("workspaces/ws_signal.root")
fInput1.cd()
workspaceSignal = fInput1.Get("myworkspace")
workspaceSignal.Print() 
fInput2 = ROOT.TFile("workspaces/ws_dataCombined.root")
fInput2.cd()
workspaceDataset = fInput2.Get("myworkspace")
fInput3 = ROOT.TFile("nEvents_Container.root")
fInput3.cd()

#VAR and CONSTANTS SETTINGS
mass_KKg = workspaceDataset.var("mass_KKg")
workspaceSignal.var("mass_KKg").setConstant(1)
workspaceSignal.var("mean").setConstant(1)
workspaceSignal.var("alpha").setConstant(1)
workspaceSignal.var("sigma").setConstant(1)
workspaceSignal.var("enne").setConstant(1)
a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",0.1,-10.,10.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",0.,-10.,10.)
c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.,-10.,10.)
d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",0.,-10.,10.)
e_bkg = ROOT.RooRealVar("e_bkg","e_bkg",0.,-10.,10.)
mytree = fInput3.Get("tree_output")
mytree.Scan()
Ndata = mytree.nEvents_BKG
signalWeight = mytree.nEvents_Signal
print "Ndata = ",Ndata
print "Nsig = ",signalWeight

Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata,Ndata-0.5*Ndata,Ndata+0.5*Ndata)

#RIPARAMETRIZATION
cross_sig = ROOT.RooRealVar("cross_sig","The signal cross section",48.58)#pb
lumi = ROOT.RooRealVar("lumi","The luminosity",35340)#pb^-1
eff = ROOT.RooRealVar("eff","efficiency",0.0470818) # = 1302/27654
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",10**(-5),10**(-7),10**(-3))
Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3",ROOT.RooArgList(cross_sig,lumi,eff,B_R_))

#PDF
sigPDF = workspaceSignal.pdf("signalPDF")
bkgPDF = ROOT.RooChebychev("bkgPDF","bkgPDF",mass_KKg,ROOT.RooArgList(a_bkg,b_bkg,c_bkg,d_bkg,e_bkg))

#ADD PDFs
argListPDFs = ROOT.RooArgList(sigPDF,bkgPDF)
argListN = ROOT.RooArgList(Nsig,Nbkg)
totPDF = ROOT.RooAddPdf("totPDF","The total PDF",argListPDFs,argListN)

#DATASET
dataCombined = workspaceDataset.data("totPDFData")
print "dataCombined retrieved!"

#FIT
totPDF.fitTo(dataCombined)
print "dataCombined fitted!"

xframe = mass_KKg.frame(40)
dataCombined.plotOn(xframe)
totPDF.plotOn(xframe)

c1 = ROOT.TCanvas()
c1.cd()
xframe.Draw()
c1.SaveAs("fit_dataCombined.pdf")

#CREATE the OUTPUT
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totPDF)

fOut = ROOT.TFile("workspaces/ws_totPDF.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()

del workspace
