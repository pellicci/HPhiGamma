import ROOT
import argparse
import math
import numpy as np
import sys
import copy
import tdrstyle, CMS_lumi
import array as array


#bools
debug = False #Bool for verbose
isPhiAnalysis = True   

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('is_Data', help='Type which kind of sample it is')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
mytree = fInput.Get("HPhiGammaTwoProngsTriggerAnalysis/mytree")

if args.is_Data == "Data":
    CONST_NAME = "Data"
else:
    CONST_NAME = "MC"

samplename =(args.rootfile_name.split("HPhiGammaTwoProngsTriggerAnalysis_output_")[1])[:-5] 
print "samplename =", samplename

#Output
output_filename = args.outputfile_option
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Normalization for MC dataset ################################################################################
if CONST_NAME == "MC":
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

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{#mu#mu} in "+CONST_NAME+" (TwoProngs leg)",100,65.,115.)
histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"p_{T}^{prongs} in "+CONST_NAME+" (TwoProngs leg)", 100, 30.,160.)
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"n. events TwoProngs35 triggered as function of p_{T}^{prongs} in "+CONST_NAME+" (TwoProngs leg)", 10, 35.,90.)
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"n. events two prongs over offline selection as function of p_{T}^{prongs} in "+CONST_NAME+" (TwoProngs leg)", 10, 35.,90.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"Trigger efficiency as function of p_{T}^{prongs} in "+CONST_NAME+" (TwoProngs leg)", 10, 35.,90.)
if isPhiAnalysis:
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+CONST_NAME+" for events at the denominator (TwoProngs leg)", 30, 1.,1.05)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"Inv mass of the pair in "+CONST_NAME+" for events at the numerator (TwoProngs leg)", 30, 1.,1.05)
else:
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+CONST_NAME+" for events at the denominator (TwoProngs leg)", 30, 0.5,1.)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"Inv mass of the pair in "+CONST_NAME+" for events at the numerator (TwoProngs leg)", 30, 0.5,1.)

#    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"Inv mass of the pair in "+CONST_NAME+" (TwoProngs leg)", 100, 0.5,1.)

#VARIABLES INITIALIZATION ##################################################################################################
nEventsIsoMuTrigger       = 0
nEventsTwoProngsTriggered = 0
nEventsOfflineTwoProngs   = 0
triggerFraction           = 0

#EVENTS LOOP ##################################################################################################################### 
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
    if CONST_NAME == "MC": 
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
    histo_map["h_MesonMass"].Fill(MesonMass, weight)
    
    if isTwoProngsTrigger: 
            histo_map["h_nTwoProngsTriggered_pT"].Fill(TwoProngsPt, weight)
            histo_map["h_triggerEff_TwoProngsPt"].Fill(TwoProngsPt, weight) #these two histos contains the same entries. h_triggerEff_TwoProngsPt is the numerator of the trigger efficiency ratio so far.
            histo_map["h_MesonMass_triggered"].Fill(MesonMass, weight)

    if isBestPairFound: 
        histo_map["h_nTwoProngsOffline_pT"].Fill(TwoProngsPt,weight)


#FIT DITRACK INV MASS ###################################################################################################

#observable
mass = ROOT.RooRealVar("mass","K^{+}K^{-} invariant mass",1.,1.05)

#signal parametrization
mean  = ROOT.RooRealVar("mean","The mean of the Crystal Ball",1.02,1.0,1.05)
sigma = ROOT.RooRealVar("sigma","The width of the Crystal Ball",0.3,0.0001,1.)
alpha = ROOT.RooRealVar("alpha","The alpha of the Crystal Ball",1.5,-5.,5.)
n     = ROOT.RooRealVar("n","The alpha of the Crystal Ball",1.5,0.5,5.)
signalPDF = ROOT.RooCBShape("signalPDF","The Crystall Ball",mass,mean,sigma,alpha,n)

#Background parametrization: 
a1 = ROOT.RooRealVar("a1","The a1 of background",0.7,-2.,2.)
a2 = ROOT.RooRealVar("a2","The a2 of background",0.3,-2.,2.)
a3 = ROOT.RooRealVar("a3","The a3 of background",-0.03,-2.,2.)
backgroundPDF = ROOT.RooChebychev("backgroundPDF","The background PDF",mass,ROOT.RooArgList(a1,a2,a3))

#FIT TO THE OFFLINE SELECTION ---------------------------------------------------------------------------------------------------------
#Define the yields
NsigOffline   = ROOT.RooRealVar("Nsig","The Phi signal events",700.,0.,5000.)
NbkgOffline   = ROOT.RooRealVar("Nbkg","The bkg events",2400.,0.,10000.)
NsigTriggered = ROOT.RooRealVar("Nsig","The Phi signal events",700.,0.,5000.)
NbkgTriggered = ROOT.RooRealVar("Nbkg","The bkg events",2400.,0.,10000.)

#Compose the total PDF
totPDFOffline = ROOT.RooAddPdf("totPDFOffline","The total PDF offline",ROOT.RooArgList(signalPDF,backgroundPDF),ROOT.RooArgList(NsigOffline,NbkgOffline))
totPDFTriggered = ROOT.RooAddPdf("totPDFTriggered","The total PDF Triggered",ROOT.RooArgList(signalPDF,backgroundPDF),ROOT.RooArgList(NsigTriggered,NbkgTriggered))

#dataset
datasetOffline   = ROOT.RooDataHist("datasetOffline","datasetOffline",ROOT.RooArgList(mass),histo_map["h_MesonMass"])
datasetTriggered = ROOT.RooDataHist("datasetTriggered","datasetTriggered",ROOT.RooArgList(mass),histo_map["h_MesonMass_triggered"])

totPDFOffline.fitTo(datasetOffline,ROOT.RooFit.Extended(1))
totPDFTriggered.fitTo(datasetTriggered,ROOT.RooFit.Extended(1))

xframeOffline = mass.frame()
datasetOffline.plotOn(xframeOffline)
totPDFOffline.plotOn(xframeOffline)

xframeTriggered = mass.frame()
datasetTriggered.plotOn(xframeTriggered)
totPDFTriggered.plotOn(xframeTriggered)

#One can also plot the single components of the total PDF, like the background component
totPDFOffline.plotOn(xframeOffline, ROOT.RooFit.Components("backgroundPDF"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
totPDFTriggered.plotOn(xframeTriggered, ROOT.RooFit.Components("backgroundPDF"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

#Draw the results
canvasOffline = ROOT.TCanvas()
xframeOffline.Draw()
canvasOffline.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/fit_mKK_Offline_"+CONST_NAME+".png")
canvasOffline.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/fit_mKK_Offline_"+CONST_NAME+".pdf")
NsigOffline.Print()

denominator = NsigOffline.getValV()
denError    = NsigOffline.getError()
print "NsigOffline = ",denominator
print "error = ",denError

#FIT TO TRIGGERING EVENTS ---------------------------------------------------------------------------------------------------------
#Define the yields
NsigTriggered = ROOT.RooRealVar("Nsig","The Phi signal events",700.,0.,5000.)
NbkgTriggered = ROOT.RooRealVar("Nbkg","The bkg events",2400.,0.,10000.)

#Compose the total PDF
totPDFTriggered = ROOT.RooAddPdf("totPDFTriggered","The total PDF Triggered",ROOT.RooArgList(signalPDF,backgroundPDF),ROOT.RooArgList(NsigTriggered,NbkgTriggered))

#dataset
datasetTriggered = ROOT.RooDataHist("datasetTriggered","datasetTriggered",ROOT.RooArgList(mass),histo_map["h_MesonMass_triggered"])

totPDFTriggered.fitTo(datasetTriggered,ROOT.RooFit.Extended(1))

xframeTriggered = mass.frame()
datasetTriggered.plotOn(xframeTriggered)
totPDFTriggered.plotOn(xframeTriggered)

#One can also plot the single components of the total PDF, like the background component
totPDFTriggered.plotOn(xframeTriggered, ROOT.RooFit.Components("backgroundPDF"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

#Draw the results
canvasTriggered = ROOT.TCanvas()
xframeTriggered.Draw()
canvasTriggered.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/fit_mKK_Triggered_"+CONST_NAME+".png")
canvasTriggered.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/fit_mKK_Triggered_"+CONST_NAME+".pdf")
NsigTriggered.Print()

numerator = NsigTriggered.getValV()
numError  = NsigTriggered.getError()
print "NsigTriggered = ",numerator
print "error = ",numError

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
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()


c1 = ROOT.TCanvas()
c1.cd()
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
ROOT.gStyle.SetOptStat(0)
histo_map["h_mass_mumu"].SetMarkerSize(1.4)
histo_map["h_mass_mumu"].Draw("")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_TwoProngsLeg_mass_mumu_"+CONST_NAME+".pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_TwoProngsLeg_mass_mumu_"+CONST_NAME+".png")  


c2 = ROOT.TCanvas()
c2.cd()
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
ROOT.gStyle.SetOptStat(0)
histo_map["h_TwoProngs_pT"].SetMarkerSize(1.4)
histo_map["h_TwoProngs_pT"].Draw("")
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_TwoProngs_pT_"+CONST_NAME+".png") 
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_TwoProngs_pT_"+CONST_NAME+".pdf") 


c3 = ROOT.TCanvas()
c3.cd()
ROOT.gStyle.SetOptStat(0)
histo_map["h_nTwoProngsTriggered_pT"].Draw("")
c3.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_nTwoProngsTriggered_pT_"+CONST_NAME+".pdf")
c3.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_nTwoProngsTriggered_pT_"+CONST_NAME+".png") 


c4 = ROOT.TCanvas()
c4.cd()
ROOT.gStyle.SetOptStat(0)
histo_map["h_nTwoProngsOffline_pT"].Draw("")
c4.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_nTwoProngsOffline_pT_"+CONST_NAME+".pdf")
c4.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_nTwoProngsOffline_pT_"+CONST_NAME+".png") 


c5 = ROOT.TCanvas()
c5.cd()
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
histo_map["h_triggerEff_TwoProngsPt"].GetYaxis().SetRangeUser(0.,1.)
histo_map["h_triggerEff_TwoProngsPt"].Draw("")
c5.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_triggerEff_TwoProngsPt_"+CONST_NAME+".pdf")
c5.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_triggerEff_TwoProngsPt_"+CONST_NAME+".png") 

c6 = ROOT.TCanvas()
c6.cd()
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
histo_map["h_MesonMass"].SetMinimum(0.)
histo_map["h_MesonMass"].Draw("")
c6.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_MesonMass_"+CONST_NAME+".pdf")
c6.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_MesonMass_"+CONST_NAME+".png")

c7 = ROOT.TCanvas()
c7.cd()
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
histo_map["h_MesonMass_triggered"].SetMinimum(0.)
histo_map["h_MesonMass_triggered"].Draw("")
c7.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_MesonMass_triggered_"+CONST_NAME+".pdf")
c7.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_latest_production/h_MesonMass_triggered_"+CONST_NAME+".png")

fOut.Close()
