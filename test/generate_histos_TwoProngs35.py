import ROOT
import argparse
import math
import numpy as np
import sys
import copy
import tdrstyle, CMS_lumi
from ROOT import gROOT
import os, sys


#bools
debug = False #Bool for verbose
isPhiAnalysis = True   

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)  

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('is_Data', help='Type which kind of sample it is')
p.add_argument('binNumber', help='Type the pT bin you want to calculate the SF to')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
mytree = fInput.Get("HPhiGammaTwoProngsTriggerAnalysis/mytree")

if args.is_Data == "Data":
    SAMPLE_NAME = "Data"
else:
    SAMPLE_NAME = "MC"

BIN_NUMBER = args.binNumber
binIndex = int(BIN_NUMBER)

samplename =(args.rootfile_name.split("HPhiGammaTwoProngsTriggerAnalysis_output_")[1])[:-5] 
print "samplename =", samplename

#Output
output_filename = args.outputfile_option
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

#Normalization for MC dataset ################################################################################
if SAMPLE_NAME == "MC":
    normalization_InputFile = open("rootfiles/latest_production/MC/normalizations/Normalizations_table_TwoProngs.txt","r")
    norm_map = dict()
    for line in normalization_InputFile:
        data_norm = line.split()
        norm_map[data_norm[0]] = float(data_norm[1])

    normalization_weight = norm_map[samplename]
    print "normalization_weight = ",normalization_weight

#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_mass_mumu","h_TwoProngs_pT","h_nTwoProngsTriggered_pT","h_nTwoProngsOffline_pT","h_triggerEff_TwoProngsPt","h_MesonMass","h_MesonMass_triggered"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{#mu#mu} in "+SAMPLE_NAME+" (TwoProngs leg)",100,65.,115.)
histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 100, 30.,160.)
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"n. events TwoProngs35 triggered as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"n. events two prongs over offline selection as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"Trigger efficiency as function of p_{T}^{prongs} in "+SAMPLE_NAME+" (TwoProngs leg)", 10, 35.,90.)
if isPhiAnalysis:
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+SAMPLE_NAME+" for events at the denominator (TwoProngs leg)", 30, 1.,1.05)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"Inv mass of the pair in "+SAMPLE_NAME+" for events at the numerator (TwoProngs leg)", 30, 1.,1.05)
else:
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+SAMPLE_NAME+" for events at the denominator (TwoProngs leg)", 30, 0.5,1.)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"Inv mass of the pair in "+SAMPLE_NAME+" for events at the numerator (TwoProngs leg)", 30, 0.5,1.)

#    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+SAMPLE_NAME+" (TwoProngs leg)", 100, 0.5,1.)

#VARIABLES INITIALIZATION ##################################################################################################
nEventsIsoMuTrigger       = 0
nEventsTwoProngsTriggered = 0
nEventsOfflineTwoProngs   = 0
triggerFraction           = 0

#EVENTS LOOP ##################################################################################################################### 

#Create lists of bin min e max values
binList = [38.,43.,48.,53.,100.]
#define pT bin
pTmin = binList[binIndex]
pTmax = binList[binIndex + 1]

print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        continue

    #Retrieve variables from the tree and fill the histos
    MuMuMass              = mytree.MuMuMass
    TwoProngsPt           = mytree.bestCouplePt
    isIsoMuTrigger        = mytree.isIsoMuTrigger
    isTwoProngsTrigger    = mytree.isTwoProngsTrigger      
    isBestPairFound       = mytree.isBestCoupleOfTheEvent_Found
    isPhi                 = mytree.isPhi
    isRho                 = mytree.isRho
    MesonMass             = mytree.MesonMass
    bestPairIsoCh         = mytree.iso_couple_ch
    if SAMPLE_NAME == "MC": 
        PUweight = mytree.PU_Weight


    #calculate PUweights
    if args.is_Data == "Data":
        weight = 1
    elif args.is_Data == "MC":
        weight = PUweight 

    #weight = 1. #FIXME
    
    #print "isBestPairFound = ",isBestPairFound
    if not isBestPairFound: continue
    
    if isPhiAnalysis:
        if isRho : 
            continue
    else:        
        if isPhi : 
            continue

    if bestPairIsoCh < 0.9: continue

    if debug:
        print ""
        print "Processing EVENT n.",jentry+1,"================"
        print "isIsoMuTrigger = ",isIsoMuTrigger
        print "isTwoProngsTrigger = ",isTwoProngsTrigger
        print "isBestPairFound = ",isBestPairFound
        print "Event weight = ",weight
    
    if isIsoMuTrigger:
        nEventsIsoMuTrigger       = nEventsIsoMuTrigger  + 1
    if isTwoProngsTrigger:
        nEventsTwoProngsTriggered = nEventsTwoProngsTriggered + 1
    if isBestPairFound:
        nEventsOfflineTwoProngs   = nEventsOfflineTwoProngs + 1
    
    if debug:
        print "total nEventsTwoProngsTriggered = ",nEventsTwoProngsTriggered
        print "total nEventsOfflineTwoProngs = ",nEventsOfflineTwoProngs

    histo_map["h_mass_mumu"].Fill(MuMuMass, weight)        
    histo_map["h_TwoProngs_pT"].Fill(TwoProngsPt, weight)
    if (TwoProngsPt >= pTmin and TwoProngsPt < pTmax):
        histo_map["h_MesonMass"].Fill(MesonMass, weight)
    
    if isTwoProngsTrigger: 
            histo_map["h_nTwoProngsTriggered_pT"].Fill(TwoProngsPt, weight)
            histo_map["h_triggerEff_TwoProngsPt"].Fill(TwoProngsPt, weight) #these two histos contains the same entries. h_triggerEff_TwoProngsPt is the numerator of the trigger efficiency ratio so far.
            if (TwoProngsPt >= pTmin and TwoProngsPt < pTmax):
                histo_map["h_MesonMass_triggered"].Fill(MesonMass, weight)

    if isBestPairFound: 
        histo_map["h_nTwoProngsOffline_pT"].Fill(TwoProngsPt,weight)


#FIT DITRACK INV MASS ###################################################################################################
#observable
mass = ROOT.RooRealVar("mass","K^{+}K^{-} invariant mass",1.,1.05)

#signal parametrization
mean  = ROOT.RooRealVar("mean","The mean of the Crystal Ball",1.02,1.0,1.05)
sigma = ROOT.RooRealVar("sigma","The width of the Crystal Ball",0.002,0.0001,0.1)
alpha = ROOT.RooRealVar("alpha","The alpha of the Crystal Ball",1.2,-2.,2.)
n     = ROOT.RooRealVar("n","The n of the Crystal Ball",4.,0.,10.)
signalPDF = ROOT.RooCBShape("signalPDF","The Crystall Ball",mass,mean,sigma,alpha,n)

#Background parametrization: 
a1 = ROOT.RooRealVar("a1","The a1 of background",0.3,-2.,2.)
a2 = ROOT.RooRealVar("a2","The a2 of background",0.3,-2.,2.)
a3 = ROOT.RooRealVar("a3","The a3 of background",-0.03,-2.,2.)
#a4 = ROOT.RooRealVar("a4","The a4 of background",-0.03,-2.,2.)
backgroundPDF = ROOT.RooChebychev("backgroundPDF","The background PDF",mass,ROOT.RooArgList(a1,a2,a3))

#FIT TO THE OFFLINE SELECTION ---------------------------------------------------------------------------------------------------------
#Define the yields
NsigOffline   = ROOT.RooRealVar("Nsig","The Phi signal events",300.,0.,1000.)
NbkgOffline   = ROOT.RooRealVar("Nbkg","The bkg events",1500.,0.,10000.)

#Compose the total PDF
totPDFOffline = ROOT.RooAddPdf("totPDFOffline","The total PDF offline",ROOT.RooArgList(signalPDF,backgroundPDF),ROOT.RooArgList(NsigOffline,NbkgOffline))

#dataset
datasetOffline = ROOT.RooDataHist("datasetOffline","datasetOffline",ROOT.RooArgList(mass),histo_map["h_MesonMass"])

result_dataFitOffline = totPDFOffline.fitTo(datasetOffline,ROOT.RooFit.Extended(1),ROOT.RooFit.Save())

xframeOffline = mass.frame()

datasetOffline.plotOn(xframeOffline)

#One can also plot the single components of the total PDF, like the background component
totPDFOffline.plotOn(xframeOffline, ROOT.RooFit.Components("backgroundPDF"), ROOT.RooFit.Name("backgroundPDF"),ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

totPDFOffline.plotOn(xframeOffline, ROOT.RooFit.Name("totPDFOffline"),ROOT.RooFit.LineColor(ROOT.kBlue))

nParam = result_dataFitOffline.floatParsFinal().getSize()
chi2 = xframeOffline.chiSquare(nParam)#Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2 = "{:.2f}".format(chi2) #Crop the chi2 to 2 decimal digits

#DRAW THE RESULTS -----------------------------------------------------------------------------------------
canvasOffline = ROOT.TCanvas()
canvasOffline.cd()

#Pads --------------------------------------
pad1 = ROOT.TPad("pad1","pad1",0,0.28,1,1.)
pad2 = ROOT.TPad("pad2","pad2",0,0.01,1,0.27)
pad1.Draw()
pad2.Draw()
xframeOffline.GetYaxis().SetTitle("Events")
xframeOffline.GetYaxis().SetTitleOffset(0.9)
xframeOffline.GetXaxis().SetLabelOffset(999)
xframeOffline.SetMaximum(1.9 * histo_map["h_MesonMass"].GetMaximum())
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

CMS_lumi.CMS_lumi(canvasOffline, iPeriod, iPos) #Print integrated lumi and energy information

#Pull sub plot ------------------------------------
pullHist  = xframeOffline.pullHist()
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
framePull.GetYaxis().SetRangeUser(-5.,5.)
framePull.GetYaxis().SetNdivisions(502,ROOT.kFALSE)
framePull.GetXaxis().SetLabelSize(0.10)
framePull.GetXaxis().SetTitleSize(0.12)
framePull.GetXaxis().SetTitleOffset(1.0)
framePull.GetXaxis().SetTitle("m_{KK} [Gev]")

#Dashed blue line of the sub plot ---------------------------------------
line_on_one = ROOT.TLine(framePull.GetXaxis().GetXmin(),0.,framePull.GetXaxis().GetXmax(),0.)
line_on_one.SetLineColor(4)
line_on_one.SetLineStyle(2)

pad1.cd()
xframeOffline.Draw()

#Legend ----------------------------------------
leg1 = ROOT.TLegend(0.65,0.62,0.87,0.90) #right positioning
leg1.SetHeader(" ")
leg1.SetNColumns(1)
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.AddEntry(histo_map["h_MesonMass"],"Data","elp")
leg1.AddEntry("totPDFOffline","Fit","l")
leg1.AddEntry("backgroundPDF","Bkg only fit","l")
leg1.AddEntry(cut_chi2,"#chi^{2}/ndof = " + cut_chi2,"brNDC")
leg1.AddEntry(NsigOffline,"N_{sig} = "+str(round(NsigOffline.getValV(),1))+" #pm "+str(round(NsigOffline.getError(),1)),"brNDC")
leg1.Draw()

leg2 = ROOT.TLegend(0.12,0.65,0.44,0.84) #left positioning
leg2.SetHeader(" ")
leg2.SetNColumns(1)
leg2.SetFillColorAlpha(0,0.)
leg2.SetBorderSize(0)
leg2.SetLineColor(1)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillStyle(1001)
leg2.AddEntry(0,"bin "+BIN_NUMBER,"")
leg2.AddEntry(0,str(pTmin)+" GeV < p_{T}^{KK} < "+str(pTmax)+" GeV","")
leg2.Draw()

pad2.cd()
framePull.Draw("E1,X0")
line_on_one.Draw()
canvasOffline.Update()
#labelOffline.Draw("SAME")
canvasOffline.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/fit_mKK_bin"+BIN_NUMBER+"_"+SAMPLE_NAME+"_Offline.png")
canvasOffline.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/fit_mKK_bin"+BIN_NUMBER+"_"+SAMPLE_NAME+"_Offline.pdf")

#Print fit results ---------------------
NsigOffline.Print()
denominator = NsigOffline.getValV()
denError    = NsigOffline.getError()
print "NsigOffline = ",denominator
print "error = ",denError


#FIT TO TRIGGERING EVENTS ---------------------------------------------------------------------------------------------------------
#Define the yields
#observable
mass = ROOT.RooRealVar("mass","K^{+}K^{-} invariant mass",1.,1.05)

#signal parametrization
mean  = ROOT.RooRealVar("mean","The mean of the Crystal Ball",1.02,1.0,1.05)
sigma = ROOT.RooRealVar("sigma","The width of the Crystal Ball",0.002,0.0001,0.1)
alpha = ROOT.RooRealVar("alpha","The alpha of the Crystal Ball",1.2,-2.,2.)
n     = ROOT.RooRealVar("n","The n of the Crystal Ball",4.,0.,10.)
signalPDF = ROOT.RooCBShape("signalPDF","The Crystall Ball",mass,mean,sigma,alpha,n)

#Background parametrization: 
a1 = ROOT.RooRealVar("a1","The a1 of background",0.3,-2.,2.)
a2 = ROOT.RooRealVar("a2","The a2 of background",0.3,-2.,2.)
a3 = ROOT.RooRealVar("a3","The a3 of background",-0.03,-2.,2.)
#a4 = ROOT.RooRealVar("a4","The a4 of background",-0.03,-2.,2.)
backgroundPDF = ROOT.RooChebychev("backgroundPDF","The background PDF",mass,ROOT.RooArgList(a1,a2,a3))

NsigTriggered = ROOT.RooRealVar("Nsig","The Phi signal events",300.,0.,1000.)
NbkgTriggered = ROOT.RooRealVar("Nbkg","The bkg events",1500.,0.,10000.)

#Compose the total PDF
totPDFTriggered = ROOT.RooAddPdf("totPDFTriggered","The total PDF Triggered",ROOT.RooArgList(signalPDF,backgroundPDF),ROOT.RooArgList(NsigTriggered,NbkgTriggered))

#dataset
datasetTriggered = ROOT.RooDataHist("datasetTriggered","datasetTriggered",ROOT.RooArgList(mass),histo_map["h_MesonMass_triggered"])

result_dataFitTriggered = totPDFTriggered.fitTo(datasetTriggered,ROOT.RooFit.Extended(1))

xframeTriggered = mass.frame()
datasetTriggered.plotOn(xframeTriggered)

#One can also plot the single components of the total PDF, like the background component
totPDFTriggered.plotOn(xframeTriggered, ROOT.RooFit.Components("backgroundPDF"),ROOT.RooFit.Name("backgroundPDF"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

totPDFTriggered.plotOn(xframeTriggered,ROOT.RooFit.Name("totPDFTriggered"),ROOT.RooFit.LineColor(ROOT.kBlue))

#Draw the results
canvasTriggered = ROOT.TCanvas()
canvasTriggered.cd()

#Pads --------------------------------------
pad1 = ROOT.TPad("pad1","pad1",0,0.28,1,1.)
pad2 = ROOT.TPad("pad2","pad2",0,0.01,1,0.27)
pad1.Draw()
pad2.Draw()
xframeTriggered.GetYaxis().SetTitle("Events")
xframeTriggered.GetYaxis().SetTitleOffset(0.9)
xframeTriggered.GetXaxis().SetLabelOffset(999)
xframeTriggered.SetMaximum(1.9 * histo_map["h_MesonMass"].GetMaximum())
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

CMS_lumi.CMS_lumi(canvasTriggered, iPeriod, iPos) #Print integrated lumi and energy information

#Pull sub plot ------------------------------------
pullHist  = xframeTriggered.pullHist()
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
framePull.GetYaxis().SetRangeUser(-5.,5.)
framePull.GetYaxis().SetNdivisions(502,ROOT.kFALSE)
framePull.GetXaxis().SetLabelSize(0.10)
framePull.GetXaxis().SetTitleSize(0.12)
framePull.GetXaxis().SetTitleOffset(1.0)
framePull.GetXaxis().SetTitle("m_{KK} [Gev]")

#Dashed blue line of the sub plot ---------------------------------------
line_on_one = ROOT.TLine(framePull.GetXaxis().GetXmin(),0.,framePull.GetXaxis().GetXmax(),0.)
line_on_one.SetLineColor(4)
line_on_one.SetLineStyle(2)

pad1.cd()
xframeTriggered.Draw()

#Legend ----------------------------------------
leg1 = ROOT.TLegend(0.65,0.62,0.87,0.90) #right positioning
leg1.SetHeader(" ")
leg1.SetNColumns(1)
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.AddEntry(histo_map["h_MesonMass"],"Data","elp")
leg1.AddEntry("totPDFTriggered","Fit","l")
leg1.AddEntry("backgroundPDF","Bkg only fit","l")
leg1.AddEntry(cut_chi2,"#chi^{2}/ndof = " + cut_chi2,"brNDC")
leg1.AddEntry(NsigTriggered,"N_{sig} = "+str(round(NsigTriggered.getValV(),1))+" #pm "+str(round(NsigTriggered.getError(),1)),"brNDC")
leg1.Draw()

leg2 = ROOT.TLegend(0.12,0.65,0.44,0.84) #left positioning
leg2.SetHeader(" ")
leg2.SetNColumns(1)
leg2.SetFillColorAlpha(0,0.)
leg2.SetBorderSize(0)
leg2.SetLineColor(1)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillStyle(1001)
leg2.AddEntry(0,"bin "+BIN_NUMBER,"")
leg2.AddEntry(0,str(pTmin)+" GeV < p_{T}^{KK} < "+str(pTmax)+" GeV","")
leg2.Draw()

pad2.cd()
framePull.Draw("E1,X0")
line_on_one.Draw()
canvasTriggered.Update()
#labelTriggered.Draw("SAME")
canvasTriggered.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/fit_mKK_bin"+BIN_NUMBER+"_"+SAMPLE_NAME+"_Triggered.png")
canvasTriggered.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/fit_mKK_bin"+BIN_NUMBER+"_"+SAMPLE_NAME+"_Triggered.pdf")
NsigTriggered.Print()

numerator = NsigTriggered.getValV()
numError  = NsigTriggered.getError()
print "NsigTriggered = ",numerator
print "error = ",numError

#Trigger eff calculation -----------------------------------------------------
triggerEff = numerator/denominator
#effError   = math.sqrt(triggerEff * (1 - triggerEff) * (1/numerator))
effErrorAlternative = math.sqrt((numError/numerator)**2 + (denError/denominator)**2)

print "triggerEff = ",triggerEff
#print "effError = ",effError
print "effErrorAlternative = ",effErrorAlternative

#TRIGGER EFFICIENCY CALCULATION ###################################################################################################
triggerFraction = float(nEventsTwoProngsTriggered)/float(nEventsOfflineTwoProngs)

#trigger efficiency as function of two prongs transverse momentum -------------------------------------------------------------
histo_map["h_triggerEff_TwoProngsPt"].Divide(histo_map["h_nTwoProngsOffline_pT"]) #divide it by the denominator
#histo_map["h_triggerEff_TwoProngsPt"].Scale(100.) #multiply for the percentage


#CALCULATE EFFICIENCY ERRORS ########################################################
print "##########################################"
for bin_index in range(1,10):
    
    N = histo_map["h_nTwoProngsOffline_pT"].GetBinContent(bin_index)
    epsilon = histo_map["h_triggerEff_TwoProngsPt"].GetBinContent(bin_index)
    if N == 0 : 
        error = 0
    else: 
        error = math.sqrt(epsilon * (1 - epsilon) * (1/N))
    
    print "For bin ",bin_index,": N = ",N,", eff = ",epsilon, " and efficiency error = ",error

    histo_map["h_triggerEff_TwoProngsPt"].SetBinError(bin_index,error)
print "##########################################"

if debug:
    print "########################################"
    print "nEventsTwoProngsTriggered = ",nEventsTwoProngsTriggered
    print "nEventsOfflineTwoProngs = ",nEventsOfflineTwoProngs
    print "Trigger fraction = ",triggerFraction
    print "########################################"

#HISTO LABELS #####################################################################################################################
histo_map["h_mass_mumu"].GetXaxis().SetTitle("m_{#mu^{+}#mu^{-}} [GeV]")
histo_map["h_TwoProngs_pT"].GetXaxis().SetTitle("p_{T}^{meson} [GeV]")

histo_map["h_nTwoProngsTriggered_pT"].GetXaxis().SetTitle("p_{T}^{meson} [GeV]")
histo_map["h_nTwoProngsTriggered_pT"].GetYaxis().SetTitle("Events")

histo_map["h_nTwoProngsOffline_pT"].GetXaxis().SetTitle("p_{T} ^{meson}[GeV]")
histo_map["h_nTwoProngsOffline_pT"].GetYaxis().SetTitle("Events")

histo_map["h_triggerEff_TwoProngsPt"].GetXaxis().SetTitle("p_{T} ^{meson} [GeV]")
histo_map["h_triggerEff_TwoProngsPt"].GetYaxis().SetTitle("Events")

histo_map["h_MesonMass"].GetXaxis().SetTitle("m_{meson} [GeV]")
histo_map["h_MesonMass"].GetYaxis().SetTitle("Events")

histo_map["h_MesonMass_triggered"].GetXaxis().SetTitle("m_{meson} [GeV]")
histo_map["h_MesonMass_triggered"].GetYaxis().SetTitle("Events")


# HISTOS WRITING #########################################################################
#fOut.cd()
for hist_name in list_histos:
    #histo_map[hist_name].Write()
    histo_map[hist_name].GetYaxis().SetTitleSize(0.07)
    histo_map[hist_name].GetYaxis().SetTitleOffset(0.7)
    histo_map[hist_name].GetYaxis().SetTitle("Events")
    histo_map[hist_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
    histo_map[hist_name].GetYaxis().SetMaxDigits(3)

c1 = ROOT.TCanvas()
c1.cd()
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
ROOT.gStyle.SetOptStat(0)
histo_map["h_mass_mumu"].SetMarkerSize(1.4)
histo_map["h_mass_mumu"].Draw("")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_TwoProngsLeg_mass_mumu_"+SAMPLE_NAME+".pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_TwoProngsLeg_mass_mumu_"+SAMPLE_NAME+".png")  


c2 = ROOT.TCanvas()
c2.cd()
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
ROOT.gStyle.SetOptStat(0)
histo_map["h_TwoProngs_pT"].SetMarkerSize(1.4)
histo_map["h_TwoProngs_pT"].Draw("")
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_TwoProngs_pT_"+SAMPLE_NAME+".png") 
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_TwoProngs_pT_"+SAMPLE_NAME+".pdf") 


c3 = ROOT.TCanvas()
c3.cd()
ROOT.gStyle.SetOptStat(0)
histo_map["h_nTwoProngsTriggered_pT"].Draw("")
c3.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_nTwoProngsTriggered_pT_"+SAMPLE_NAME+".pdf")
c3.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_nTwoProngsTriggered_pT_"+SAMPLE_NAME+".png") 


c4 = ROOT.TCanvas()
c4.cd()
ROOT.gStyle.SetOptStat(0)
histo_map["h_nTwoProngsOffline_pT"].Draw("")
c4.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_nTwoProngsOffline_pT_"+SAMPLE_NAME+".pdf")
c4.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_nTwoProngsOffline_pT_"+SAMPLE_NAME+".png") 


#c5 = ROOT.TCanvas()
#c5.cd()
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
#histo_map["h_triggerEff_TwoProngsPt"].GetYaxis().SetRangeUser(0.,1.)
#histo_map["h_triggerEff_TwoProngsPt"].Draw("")
#c5.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_triggerEff_TwoProngsPt_"+SAMPLE_NAME+".pdf")
#c5.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_triggerEff_TwoProngsPt_"+SAMPLE_NAME+".png") 

#c6 = ROOT.TCanvas()
#c6.cd()
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
#histo_map["h_MesonMass"].SetMinimum(0.)
#histo_map["h_MesonMass"].Draw("")
#c6.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_MesonMass_"+SAMPLE_NAME+".pdf")
#c6.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_MesonMass_"+SAMPLE_NAME+".png")

#c7 = ROOT.TCanvas()
#c7.cd()
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
#histo_map["h_MesonMass_triggered"].SetMinimum(0.)
#histo_map["h_MesonMass_triggered"].Draw("")
#c7.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_MesonMass_triggered_"+SAMPLE_NAME+".pdf")
#c7.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_MesonMass_triggered_"+SAMPLE_NAME+".png")

#OUTPUT FILE ################################################
#It contains the trigger efficiency bin per bin
fOut   = ROOT.TFile("histos/histos_trigger_efficiency/histos_TwoProngs/histos_trigger_efficency_"+SAMPLE_NAME+".root","UPDATE")
#h_trigger_efficiency = ROOT.TH1F("h_trigger_efficiency","Trigger efficiency as function of p_{T}^{prongs} in (TwoProngs leg)", len(binList)-1, np.array(binList))

if ROOT.gSystem.AccessPathName("histos/histos_trigger_efficiency/histos_TwoProngs/histos_trigger_efficency_"+SAMPLE_NAME+".root"): #ATTENTION: AccessPathName is true if the file does not exists
    print "File does not exists: I create it."
    h_trigger_efficiency = ROOT.TH1F("h_trigger_efficiency","Trigger efficiency as function of p_{T}^{prongs} in (TwoProngs leg)", len(binList)-1, np.array(binList))

else: 
    print "File exists: I update it."
    #fOut   = ROOT.TFile("histos/histos_trigger_efficiency/histos_TwoProngs/histos_trigger_efficency_"+SAMPLE_NAME+".root","UPDATE")
    fOut.Print()
    h_trigger_efficiency = fOut.Get("h_trigger_efficiency")

triggerEff = float(numerator/denominator)

print "number of bins = ",len(binList)-1
print "binIndex = ",binIndex
print "triggerEff = ",triggerEff

if triggerEff > 1.: 
    triggerEff = 1.
    print "Trigger efficinecy set to 100 %"

h_trigger_efficiency.SetBinContent(binIndex+1,triggerEff)

h_trigger_efficiency.Write(h_trigger_efficiency.GetName(),ROOT.TObject.kOverwrite)

triggerEffCanvas = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
h_trigger_efficiency.SetTitle("TwoProngs trigger efficiency in "+SAMPLE_NAME)
h_trigger_efficiency.GetXaxis().SetTitle("p_{T} ^{meson} [GeV]")
h_trigger_efficiency.GetYaxis().SetTitle("trigger eff")
h_trigger_efficiency.GetYaxis().SetRangeUser(0.3,1.2)
h_trigger_efficiency.GetXaxis().SetLabelSize(0.035)
h_trigger_efficiency.GetYaxis().SetLabelSize(0.035)
triggerEffCanvas.SetRightMargin(0.05)
ROOT.gStyle.SetPaintTextFormat("4.2f")
h_trigger_efficiency.Draw("HIST TEXT0")
CMS_lumi.CMS_lumi(triggerEffCanvas, iPeriod, iPos) #Print integrated lumi and energy information
triggerEffCanvas.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_triggerEff_TwoProngsPt_"+SAMPLE_NAME+".pdf")
triggerEffCanvas.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_triggerEff_TwoProngsPt_"+SAMPLE_NAME+".png") 

fOut.cd()
fOut.Close()

