import ROOT

#Some useful things before starting --------------------------------------------------------------------
ROOT.gROOT.SetBatch(True) #Supress the opening of many Canvas's
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+") #import the Doube Crystal Ball PDF


#Define the observable and sidebands --------------------------
mass = ROOT.RooRealVar("mass_KKg","mass_KKg",100.,150.,"GeV")
mass.setRange("LowSideband",100.,120.)
mass.setRange("HighSideband",130.,150.)


#Data input -------------------------------------------------------------
fileInputData = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
fileInputData.cd()
mytree = fileInputData.Get("tree_output")


#Import dataset ---------------------------------------------------------------------------
data = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass),ROOT.RooFit.Import(mytree))
print "Using ", data.numEntries(), " events to fit"


#Import from workspaces --------------------------------------------------------------------
fInputSignalPDF = ROOT.TFile("workspaces/ws_signal.root")
fInputSignalPDF.cd()
workspaceSignal = fInputSignalPDF.Get("myworkspace")

#Signal PDF
workspaceSignal.var("dCB_aL").setConstant(1)
workspaceSignal.var("dCB_nL").setConstant(1)
workspaceSignal.var("dCB_aR").setConstant(1)
workspaceSignal.var("dCB_nR").setConstant(1)
workspaceSignal.var("dCB_pole").setConstant(1)
workspaceSignal.var("dCB_width").setConstant(1)		
sigPDF = workspaceSignal.pdf("signalPDF") #SIGNAL PDF

#Background PDF
fInputSidebandsPDF = ROOT.TFile("workspaces/ws_sidebands.root")
fInputSidebandsPDF.cd()
workspaceSidebands = fInputSidebandsPDF.Get("myworkspace")
bkgPDF = workspaceSidebands.pdf("bkgPDF") #BKG PDF

#Number of events
#Ndata = 3617 #number of events in the sidebands of m_KKg
#signalWeight = 0.227113240117 #sum of signal weights event by event
#print "Ndata = ", Ndata
#print "Nsig = ", signalWeight


#Define the POI and an offset for blind analysis -------------------------------------------------------------------------------------
B_R = ROOT.RooRealVar("B_R","branching_ratio",10**(-5),-1,10**(-2)) #parameter of interest
B_R_blind = ROOT.RooUnblindOffset("B_R_blind","B_R_blind","aSeedString",0.000001,B_R) #generate an unkown offset to be applied to the blind window

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
cross_sig = ROOT.RooRealVar("cross_sig","Th e signal cross section",52.36)#52.36pb (gluon fusion + vector boson fusion xsec)
lumi = ROOT.RooRealVar("lumi","The luminosity",39540)#pb^-1
eff = ROOT.RooRealVar("eff","efficiency",0.058251) # = 2898/49750 = 0.058251
B_R_phiKK = ROOT.RooRealVar("B_R_phiKK","b_r_phi_KK",0.489) # = BR(phi -> K+K-)
Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3*@4",ROOT.RooArgList(cross_sig,lumi,eff,B_R_phiKK,B_R)) #for Signal
Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",700, 300, 5000) #for BKG

#ADD PDFs
argListPDFs = ROOT.RooArgList(sigPDF,bkgPDF)
argListN = ROOT.RooArgList(Nsig,Nbkg)
totPDF = ROOT.RooAddPdf("totPDF","The total PDF",argListPDFs,argListN)
print "PDFs added!"

#Do the fit --------------------------------------------------------------------------------------------------------------
result_dataFit = totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Save() ) #For the signal region, I want the fit to be extended (Poisson fluctuation of number of events) 
																										#to take into account that the total number of events is the sum of signal and background events. 
																										#Either I do this, or I use a fraction frac*Nbkg+(1-frac)*Nsig, which will become a parameter of 
																										#the fit and will have a Gaussian behavior (whilst the extended fit preserves the natural Poisson behavior)

#Plots -------------------------------------------------------------------------------------------------------------------
xframe = mass.frame(55.,95.,15)
xframe.SetTitle(" ")
xframe.SetTitleOffset(1.4,"y")
xframe.SetMaximum(90)

#Exclude the control regions
data_reduced = dataset.reduce("mass_KKg < 120. || mass_KKg > 130.")
data_reduced.plotOn(xframe)
totPDF.plotOn(xframe)

#Draw on canvas
canvas = ROOT.TCanvas()
canvas.cd()
xframe.Draw()

#Save the plot
canvas.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fitdataBlind.pdf")
canvas.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fitdataBlind.png")

#Write to file -------------------------------------------------------------------------------------------------------------------
fOutput = ROOT.TFile("workspaces/ws_data_blind.root","RECREATE")
fOutput.cd()

workspace_out = ROOT.RooWorkspace("workspace")
getattr(workspace_out,'import')(data)
getattr(workspace_out,'import')(totPDF)

workspace_out.Write()

fOutput.Close()
del workspace_out

raw_input()

#PLOT
#xframe = mass.frame(50)
#xframe.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
#xframe.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
#dataset.plotOn(xframe)
#totPDF.plotOn(xframe)

#c1 = ROOT.TCanvas()
#c1.cd()
#c1.SetTitle("")
#xframe.Draw()
#c1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fitdataBlind.pdf")
#c1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/fitdataBlind.png")

#create Workspace
#print "**************************************n. events = ",nEntries
#workspace = ROOT.RooWorkspace("myworkspace")
#getattr(workspace,'import')(totPDF)

#fOut = ROOT.TFile("workspaces/ws_data_blind.root","RECREATE")
#fOut.cd()
#workspace.Write()
#fOut.Close()
