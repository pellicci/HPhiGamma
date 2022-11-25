import ROOT
import tdrstyle, CMS_lumi
import argparse

#bool initialization
isRhoGammaAnalysis = False
isPhiGammaAnalysis = False

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('Decay_channel_option', help='Type <<Phi>> for Phi, <<Rho>> for Rho') #flag for bkg estimation
args = p.parse_args()

if args.Decay_channel_option == "Phi":
	isPhiGammaAnalysis = True
	CHANNEL = "Phi"
	print "H -> PhiGamma analysis"
else: 
    isRhoGammaAnalysis = True
    CHANNEL = "Rho"
    print "H -> RhoGamma analysis"

fileInput_ggH = ROOT.TFile("histos/latest_production/histos_SR_SignalggH.root")
fileInput_VBF = ROOT.TFile("histos/latest_production/histos_SR_SignalVBF.root")

#import the Doube Crystal Ball PDF
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+")

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

#Define the observable
mass = ROOT.RooRealVar("mesonGammaMass","mesonGammaMass",125.,100.,170.,"GeV/c^2")

#Double Crystal Ball definition
dCB_pole_ggH  = ROOT.RooRealVar("dCB_pole_"+CHANNEL+"_GFcat_ggH", "Double CB pole", 125.,120.,130.)
dCB_width_ggH = ROOT.RooRealVar("dCB_width_"+CHANNEL+"_GFcat_ggH", "Double CB width",1.,0.5,10.)
dCB_aL_ggH    = ROOT.RooRealVar("dCB_aL_"+CHANNEL+"_GFcat_ggH", "Double CB alpha left", 3., 0.1, 50.)
dCB_aR_ggH    = ROOT.RooRealVar("dCB_aR_"+CHANNEL+"_GFcat_ggH", "Double CB alpha right", 1., 0.1, 50.)
dCB_nL_ggH    = ROOT.RooRealVar("dCB_nL_"+CHANNEL+"_GFcat_ggH", "Double CB n left", 3., 0.1, 50.)
dCB_nR_ggH    = ROOT.RooRealVar("dCB_nR_"+CHANNEL+"_GFcat_ggH", "Double CB n right", 1., 0.1, 50.)

dCB_pole_VBF  = ROOT.RooRealVar("dCB_pole_"+CHANNEL+"_GFcat_VBF", "Double CB pole", 125.,120.,130.)
dCB_width_VBF = ROOT.RooRealVar("dCB_width_"+CHANNEL+"_GFcat_VBF", "Double CB width",1.,0.5,10.)
dCB_aL_VBF    = ROOT.RooRealVar("dCB_aL_"+CHANNEL+"_GFcat_VBF", "Double CB alpha left", 1., 0.1, 50.)
dCB_aR_VBF    = ROOT.RooRealVar("dCB_aR_"+CHANNEL+"_GFcat_VBF", "Double CB alpha right", 1., 0.1, 50.)
dCB_nL_VBF    = ROOT.RooRealVar("dCB_nL_"+CHANNEL+"_GFcat_VBF", "Double CB n left", 3., 0.1, 50.)
dCB_nR_VBF    = ROOT.RooRealVar("dCB_nR_"+CHANNEL+"_GFcat_VBF", "Double CB n right", 1., 0.1, 50.)

signalPDF_ggH = ROOT.RooDoubleCBFast("crystal_ball_"+CHANNEL+"_GFcat_ggH", "Double Crystal Ball", mass, dCB_pole_ggH, dCB_width_ggH, dCB_aL_ggH, dCB_nL_ggH, dCB_aR_ggH, dCB_nR_ggH)
signalPDF_VBF = ROOT.RooDoubleCBFast("crystal_ball_"+CHANNEL+"_GFcat_VBF", "Double Crystal Ball", mass, dCB_pole_VBF, dCB_width_VBF, dCB_aL_VBF, dCB_nL_VBF, dCB_aR_VBF, dCB_nR_VBF)

#Input file and tree
fileInput_ggH.cd()
tree_ggH = fileInput_ggH.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also
dataset_ggH = ROOT.RooDataSet("dataset_ggH","dataset_ggH", ROOT.RooArgSet(mass), ROOT.RooFit.Import(tree_ggH))
#dataset_ggH = ROOT.RooDataHist("dataset_ggH","dataset_ggH", ROOT.RooArgList(mass), ROOT.RooFit.Import(fileInput_ggH.Get("h_InvMass_TwoTrk_Photon")))

#Do the fit
fitResult_ggH = signalPDF_ggH.fitTo(dataset_ggH)

#Plot
xframe_ggH = mass.frame(40)

#ChiSquare test
#nParam_ggH   = fitResult_ggH.floatParsFinal().getSize()
#chi2_ggH     = xframe_ggH.chiSquare()#Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
#cut_chi2_ggH = "{:.2f}".format(chi2_ggH) #Crop the chi2 to 2 decimal digits


dataset_ggH.plotOn(xframe_ggH)
signalPDF_ggH.plotOn(xframe_ggH)
xframe_ggH.SetTitle("")
ROOT.gStyle.SetOptFit(1111)
xframe_ggH.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_ggH.GetXaxis().SetRangeUser(100.,170.)
xframe_ggH.GetYaxis().SetMaxDigits(2)

c1 = ROOT.TCanvas()
c1.cd()
c1.SetTitle("")
xframe_ggH.Draw()
                                                                                                                                                        
'''
#######################################################################################################
#Pads --------------------------------------
pad1 = ROOT.TPad("pad1","pad1",0,0.28,1,1.)
pad2 = ROOT.TPad("pad2","pad2",0,0.01,1,0.27)
pad1.Draw()
pad2.Draw()
xframe_ggH.GetYaxis().SetTitle("Events")
xframe_ggH.GetYaxis().SetTitleOffset(0.9)
xframe_ggH.GetXaxis().SetLabelOffset(999)
xframe_ggH.GetXaxis().SetRangeUser(100.,170.)
xframe_ggH.GetYaxis().SetMaxDigits(2)
#xframe_ggH.SetMaximum(1.9 * histo_map["h_MesonMass"].GetMaximum())
pad1.SetBottomMargin(0.02)
pad1.SetTopMargin(0.07)
pad1.SetBorderMode(0)
pad1.SetBorderSize(0)
pad1.SetFrameBorderSize(0)
pad2.SetBorderSize(0)
pad2.SetFrameBorderSize(0)
pad2.SetBottomMargin(0.3)
pad2.SetBorderMode(0)
pad2.SetTopMargin(0.03)
pad1.SetRightMargin(0.05)
pad2.SetRightMargin(0.05)
pad1.SetLeftMargin(0.15)
pad2.SetLeftMargin(0.15)
pad2.SetFillColor(0)
pad2.SetFillStyle(0)
#ROOT.gStyle.SetOptStat(0)

CMS_lumi.CMS_lumi(c1, iPeriod, iPos) #Print integrated lumi and energy information

#Pull sub plot ------------------------------------
pullHist  = xframe_ggH.pullHist()
framePull = mass.frame()
framePull.addPlotable(pullHist,"P")
framePull.SetXTitle("")
framePull.SetTitle("")
framePull.SetMarkerStyle(8)
framePull.SetMarkerColor(1)
framePull.SetLineColor(1)
framePull.GetYaxis().SetLabelSize(0.1)
framePull.GetYaxis().SetTitle("Pull")
framePull.GetYaxis().SetTitleSize(0.16)
framePull.GetYaxis().SetTitleOffset(0.2)
framePull.GetYaxis().CenterTitle(True)
framePull.GetYaxis().SetRangeUser(-50.,50.)
framePull.GetYaxis().SetNdivisions(502,ROOT.kFALSE)
framePull.GetXaxis().SetLabelSize(0.10)
framePull.GetXaxis().SetTitleSize(0.12)
framePull.GetXaxis().SetTitleOffset(1.0)
framePull.GetXaxis().SetTitle("m_{KK#gamma} [Gev]")

#Dashed blue line of the sub plot ---------------------------------------
line_on_one = ROOT.TLine(framePull.GetXaxis().GetXmin(),0.,framePull.GetXaxis().GetXmax(),0.)
line_on_one.SetLineColor(4)
line_on_one.SetLineStyle(2)

pad1.cd()
xframe_ggH.Draw()
pad2.cd()
framePull.Draw("E1,X0")
line_on_one.Draw()
c1.Update()
'''
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fitsignal_"+CHANNEL+"_GFcat_ggH.pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fitsignal_"+CHANNEL+"_GFcat_ggH.png")  
#######################################################################################################

#the same for VBF
fileInput_VBF.cd()
tree_VBF = fileInput_VBF.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also
dataset_VBF = ROOT.RooDataSet("dataset_VBF","dataset_VBF",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree_VBF))

#Do the fit
fitResult_VBF = signalPDF_VBF.fitTo(dataset_VBF)

#Plot
xframe_VBF = mass.frame(40)
dataset_VBF.plotOn(xframe_VBF)
signalPDF_VBF.plotOn(xframe_VBF)
xframe_VBF.SetTitle("")
ROOT.gStyle.SetOptFit(1111)
xframe_VBF.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_VBF.GetXaxis().SetRangeUser(100.,170.)
xframe_VBF.GetYaxis().SetMaxDigits(2)

c2 = ROOT.TCanvas()
c2.cd()
c2.SetTitle("")
xframe_VBF.Draw()

c2.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fitsignal_"+CHANNEL+"_GFcat_VBF.pdf")
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/fitsignal_"+CHANNEL+"_GFcat_VBF.png")                                                                                                                                                          


#create Workspace
norm_ggH     = fileInput_ggH.Get("h_InvMass_TwoTrk_Photon").Integral() #get the normalization of ggH signal (area under ggH signal)
sig_norm_ggH = ROOT.RooRealVar(signalPDF_ggH.GetName()+ "_norm", signalPDF_ggH.GetName()+ "_norm", norm_ggH)

norm_VBF     = fileInput_VBF.Get("h_InvMass_TwoTrk_Photon").Integral() #get the normalization of VBF signal (area under VBF signal)
sig_norm_VBF = ROOT.RooRealVar(signalPDF_VBF.GetName()+ "_norm", signalPDF_VBF.GetName()+ "_norm", norm_VBF)

dCB_pole_ggH.setConstant()  
dCB_width_ggH.setConstant()  
dCB_aL_ggH.setConstant()     
dCB_aR_ggH.setConstant()    
dCB_nL_ggH.setConstant()    
dCB_nR_ggH.setConstant()     
sig_norm_ggH.setConstant()

dCB_pole_VBF.setConstant()  
dCB_width_VBF.setConstant()  
dCB_aL_VBF.setConstant()     
dCB_aR_VBF.setConstant()    
dCB_nL_VBF.setConstant()    
dCB_nR_VBF.setConstant()     
sig_norm_VBF.setConstant()

workspace = ROOT.RooWorkspace("workspace_STAT_"+CHANNEL+"_GFcat_2018")
getattr(workspace,'import')(signalPDF_ggH)
getattr(workspace,'import')(signalPDF_VBF)
getattr(workspace,'import')(sig_norm_ggH)
getattr(workspace,'import')(sig_norm_VBF)
getattr(workspace,'import')(dataset_ggH)
getattr(workspace,'import')(dataset_VBF)
print "ggH Signal integral = ",sig_norm_ggH.Print()
print "VBF Signal integral = ",sig_norm_VBF.Print()

workspace.Print()

fOut = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_2018.root","RECREATE")
fOut.cd()
workspace.Write()
xframe_ggH.Write("xframe_ggH")
xframe_VBF.Write("xframe_VBF")
fOut.Close()	

'''
#Number of events in m_KKg
fileInputSignal_ggH = ROOT.TFile("histos/latest_production/histos_SR_SignalggH.root")
fileInputSignal_ggH.cd()
mytreeSignal_ggH = fileInputSignal_ggH.Get("tree_output")
fInput_ggH = ROOT.TFile("rootfiles/latest_production/MC/signals/HPhiGammaAnalysis_Signal_"+CHANNEL+"_ggH.root")
mytree_ggH         = fInput_ggH.Get("HPhiGammaAnalysis/mytree")
h_Events_ggH       = fInput_ggH.Get("HPhiGammaAnalysis/h_Events")
NsigPassed_ggH     = mytreeSignal_ggH.GetEntriesFast()
totalSigEvents_ggH = h_Events_ggH.GetBinContent(1)
sigEfficiency_ggH  = float(NsigPassed_ggH/totalSigEvents_ggH)

fileInputSignal_VBF = ROOT.TFile("histos/latest_production/histos_SR_SignalVBF.root")
fileInputSignal_VBF.cd()
mytreeSignal_VBF = fileInputSignal_VBF.Get("tree_output")
fInput_VBF = ROOT.TFile("rootfiles/latest_production/MC/signals/HPhiGammaAnalysis_Signal_"+CHANNEL+"_VBF.root")
mytree_VBF    = fInput_VBF.Get("HPhiGammaAnalysis/mytree")
h_Events_VBF  = fInput_VBF.Get("HPhiGammaAnalysis/h_Events")
NsigPassed_VBF = mytreeSignal_VBF.GetEntriesFast()
totalSigEvents_VBF = h_Events_VBF.GetBinContent(1)
sigEfficiency_VBF = float(NsigPassed_VBF/totalSigEvents_VBF)

lumi = 39.54
x_sec_ggH = 46870.
x_sec_VBF =  3780.
BR = 10**(-5)

print "----------------------------------------------------------------------------------------------------------------"
print "For: Lumi = 39.54 fb-1, eff = ", sigEfficiency_ggH,", ggH cross section = 46870 fb and a BR injected of = 10^-5"
print "N expected events for ggH = ",lumi*sigEfficiency_ggH*x_sec_ggH*BR
print "For: Lumi = 39.54 fb-1, eff = ", sigEfficiency_ggH,", ggH cross section = 46870 fb and a BR injected of = 10^-5"
print "N expected events for VBF = ",lumi*sigEfficiency_VBF*x_sec_VBF*BR
print "Integral ggh = ",fileInput_ggH.Get("h_InvMass_TwoTrk_Photon").Integral()
'''

