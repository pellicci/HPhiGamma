import ROOT

#This script do the fit of a combined PDF (signal pdf + bkg pdf) to a generated dataset

#Some useful things before starting --------------------------------------------------------------------
ROOT.gROOT.SetBatch(True) #Supress the opening of many Canvas's
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+") #import the Doube Crystal Ball PDF

#Define the observable --------------------------
mass = ROOT.RooRealVar("mass_KKg","mass_KKg",100.,170.,"GeV")

#Data input -------------------------------------------------------------
fileInputData = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
fileInputData.cd()
mytreeData = fileInputData.Get("tree_output")

#Signal input -------------------------------------------------------------
fileInputSignal = ROOT.TFile("histos/latest_production/histos_SR_Signal.root")
fileInputSignal.cd()
mytreeSignal = fileInputSignal.Get("tree_output")

#Import from workspaces --------------------------------------------------------------------
fInputSignalPDF = ROOT.TFile("workspaces/ws_signal.root")
fInputSignalPDF.cd()
workspaceSignal = fInputSignalPDF.Get("myworkspace")
workspaceSignal.Print() #print some info of the workspace content

#Signal PDF
workspaceSignal.var("dCB_nL").setConstant(1) 
workspaceSignal.var("dCB_nR").setConstant(1)
workspaceSignal.var("dCB_aL").setConstant(1)
workspaceSignal.var("dCB_aR").setConstant(1)
workspaceSignal.var("dCB_pole").setConstant(1)
workspaceSignal.var("dCB_width").setConstant(1)

signalPDF = workspaceSignal.pdf("signalPDF") #SIGNAL PDF from the workspace

#Background PDF
fInputSidebandsPDF = ROOT.TFile("workspaces/ws_sidebands.root")
fInputSidebandsPDF.cd()
workspaceSidebands = fInputSidebandsPDF.Get("myworkspace")
bkgPDF = workspaceSidebands.pdf("bkgPDF") #BKG PDF from the workspace
workspaceSidebands.Print()

#Define the POI -------------------------------------------------------------------------------------
B_R = ROOT.RooRealVar("B_R","branching_ratio",10**(-5),-0.1,10**(-2)) #parameter of interest

#Number of events in m_KKg
fInput = ROOT.TFile("rootfiles/latest_production/MC/signals/PhiGamma_Signal/HPhiGammaAnalysis_Signal.root")
mytree   = fInput.Get("HPhiGammaAnalysis/mytree")
h_Events = fInput.Get("HPhiGammaAnalysis/h_Events")
NsigPassed = mytreeSignal.GetEntriesFast()
totalSigEvents = h_Events.GetBinContent(1)
sigEfficiency = float(NsigPassed/totalSigEvents)
Ndata = mytreeData.GetEntriesFast()

#N bkg estimation
Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata, 10., 5000.)

#BKG PDF
bkgPDF = workspaceSidebands.pdf("bkgPDF")

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
cross_sig = ROOT.RooRealVar("cross_sig","Th e signal cross section",46.87)#46.87pb (gluon fusion)
lumi = ROOT.RooRealVar("lumi","The luminosity",39540)#pb^-1
eff = ROOT.RooRealVar("eff","efficiency",sigEfficiency) # = (n sig over tightselection)/49750
B_R_phiKK = ROOT.RooRealVar("B_R_phiKK","b_r_phi_KK",0.489) # = BR(phi -> K+K-)
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",10**(-5),-0.01,10**(-2))
Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3*@4",ROOT.RooArgList(cross_sig,lumi,eff,B_R_phiKK,B_R_))

#ADD PDFs --------------------------------------------------------------------------------------------------
argListN = ROOT.RooArgList(Nsig,Nbkg)
print "arg list N size = ",argListN.getSize()
argListPDFs = ROOT.RooArgList(signalPDF,bkgPDF)
print "arg list pdf size = ",argListPDFs.getSize()

print "arg lists created!"
totPDF = ROOT.RooAddPdf("totPDF","The total PDF",argListPDFs,argListN)

print "PDFs added!"

#Generate dataset for blind analysis -----------------------------------------------------------------------------------
datasetGenerated = totPDF.generate(ROOT.RooArgSet(mass),Ndata)
print "Dataset generated!"

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

#Final useful information for debugging
print ""
print "######################################################"
print "Ndata          = ", Ndata
print "totalSigEvents = ",totalSigEvents
print "NsigPassed     = ",NsigPassed
print "Signal efficiency after tightselection = ", sigEfficiency
print "######################################################"
print ""

#create Workspace
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totPDF)
getattr(workspace,'import')(datasetGenerated)

fOut = ROOT.TFile("workspaces/ws_data_blind.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()
