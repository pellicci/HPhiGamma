import ROOT

#This script do the fit of a combined PDF (signal pdf + bkg pdf) to a generated dataset

#BOOLs for signal PDF
isCBplusGauss = True 
isDoubleCB = False

#Some useful things before starting --------------------------------------------------------------------
ROOT.gROOT.SetBatch(True) #Supress the opening of many Canvas's
if isDoubleCB:
	ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+") #import the Doube Crystal Ball PDF
	ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/GBRMath.cc+") #import the Doube Crystal Ball PDF

#Define the observable --------------------------
mass = ROOT.RooRealVar("mass_KKg","mass_KKg",100.,150.,"GeV")

#Data input -------------------------------------------------------------
fileInputData = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
fileInputData.cd()
mytree = fileInputData.Get("tree_output")

#Import from workspaces --------------------------------------------------------------------
fInputSignalPDF = ROOT.TFile("workspaces/ws_signal.root")
fInputSignalPDF.cd()
workspaceSignal = fInputSignalPDF.Get("myworkspace")
workspaceSignal.Print() #print some info of the workspace content

#Signal PDF
if isDoubleCB:
	workspaceSignal.var("Gauss_pole").setConstant(1) #Set constant these parameters of the signal PDF, assuming the MC well describes the signal shape
	workspaceSignal.var("Gauss_sigma").setConstant(1)
	workspaceSignal.var("dCB_aR").setConstant(1)
	workspaceSignal.var("dCB_nR").setConstant(1)
	workspaceSignal.var("dCB_pole").setConstant(1)
	workspaceSignal.var("dCB_width").setConstant(1)

if isCBplusGauss:
	workspaceSignal.var("mean").setConstant(1) #Set constant these parameters of the signal PDF, assuming the MC well describes the signal shape
	workspaceSignal.var("mean").setConstant(1) 
	workspaceSignal.var("mean").setConstant(1) 
	workspaceSignal.var("sigma").setConstant(1)
	workspaceSignal.var("alpha").setConstant(1)
	workspaceSignal.var("enne").setConstant(1)

signalPDF = workspaceSignal.pdf("signalPDF") #SIGNAL PDF from the workspace

#Background PDF
fInputSidebandsPDF = ROOT.TFile("workspaces/ws_sidebands.root")
fInputSidebandsPDF.cd()
workspaceSidebands = fInputSidebandsPDF.Get("myworkspace")
bkgPDF = workspaceSidebands.pdf("bkgPDF") #BKG PDF from the workspace
workspaceSidebands.Print()

#Define the POI -------------------------------------------------------------------------------------
B_R = ROOT.RooRealVar("B_R","branching_ratio",10**(-5),-1,10**(-2)) #parameter of interest

#Number of events
Ndata = 4149 #number of events in m_KKg
signalWeight = 0.227113240117
print "Ndata = ", Ndata
#print "n of data events = ",mytree.mass_KKg.GetEntries()
print "Nsig = ", signalWeight

Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata, 0.1, 5000)

#BKG PDF
bkgPDF = workspaceSidebands.pdf("bkgPDF")

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
cross_sig = ROOT.RooRealVar("cross_sig","Th e signal cross section",52.36)#52.36pb (gluon fusion + vector boson fusion xsec)
lumi = ROOT.RooRealVar("lumi","The luminosity",39540)#pb^-1
eff = ROOT.RooRealVar("eff","efficiency",0.058251) # = 2898/49750 = 0.058251
B_R_phiKK = ROOT.RooRealVar("B_R_phiKK","b_r_phi_KK",0.489) # = BR(phi -> K+K-)
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",10**(-5),-1,10**(-2))
Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3*@4",ROOT.RooArgList(cross_sig,lumi,eff,B_R_phiKK,B_R_))

#ADD PDFs --------------------------------------------------------------------------------------------------
argListN = ROOT.RooArgList(Nsig,Nbkg)
argListPDFs = ROOT.RooArgList(signalPDF,bkgPDF)
totPDF = ROOT.RooAddPdf("totPDF","The total PDF",argListPDFs,argListN)
print "PDFs added!"

#Generate dataset for blind analysis --------------------------------------------------------------------------------------------------
datasetGenerated = totPDF.generate(ROOT.RooArgSet(mass),Ndata)

#FIT
totPDF.fitTo(datasetGenerated)

#PLOT
xframe = mass.frame(40)
xframe.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
datasetGenerated.plotOn(xframe)
totPDF.plotOn(xframe)

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe.Draw()
c1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fit_to_generated_dataset.pdf")
c1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fit_to_generated_dataset.png")

#create Workspace
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totPDF)

fOut = ROOT.TFile("workspaces/ws_data_blind.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()
