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
isNorm = False

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

# INPUT ###########################################################################################################
fileData = ROOT.TFile("histos/latest_production/histos_SR_preselection_SignalVBF.root")
fileMC   = ROOT.TFile("histos/latest_production/histos_SR_preselection_SignalVBF.root")

h_Data = fileData.Get("h_bestCouplePt")
h_MC   = fileMC.Get("h_genMesonPt")

if isNorm: h_MC.Scale(h_Data.Integral()/h_MC.Integral())

#PLOT ##################################################################

#Canvas
canvas = ROOT.TCanvas()
canvas.cd()

#Legend
leg1 = ROOT.TLegend(0.68,0.65,0.82,0.82) #right positioning
leg1.SetHeader(" ")
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.SetNColumns(1)
leg1.AddEntry(h_Data,"Reco meson","ep")
leg1.AddEntry(h_MC,"Gen meson","f")


#Pads
pad1 = ROOT.TPad("pad1","",0.,0.29,1,1)
pad2 = ROOT.TPad("pad2","",0,0,1,0.27)
pad1.SetBottomMargin(0.02)
pad1.SetBorderMode(0)
pad1.SetTicks(2,1) #ticks on the right and the upper axis are drawn inside
pad2.SetTopMargin(0.01)
pad2.SetBottomMargin(0.3)
pad2.SetBorderMode(0)
h_MC.GetYaxis().SetTitle("Events")
h_MC.SetFillColor(ROOT.kOrange)
h_MC.SetLineColor(ROOT.kOrange)
h_MC.SetBinErrorOption(ROOT.TH1.kPoisson)
#h_MC.SetFillColor(ROOT.kRed)

pad1.Draw()
pad2.Draw()

pad1.cd()
#h_MC.Rebin(10)
#h_Data.Rebin(10)

h_MC.GetYaxis().SetTitleSize(0.07)
h_MC.GetYaxis().SetTitleOffset(0.7)
h_MC.GetYaxis().SetTitle("Events")
h_MC.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
h_MC.GetXaxis().SetLabelOffset(999)
h_MC.GetYaxis().SetMaxDigits(2)

h_MC.SetMaximum(1.7 * h_MC.GetMaximum())

h_MC.Draw("hist")
h_Data.Draw("SAME,E1")

#MC errors
hMCErr = copy.deepcopy(h_MC)
hMCErr_size = hMCErr.GetSize() - 2
hMCErr.SetFillStyle(3005)
hMCErr.SetMarkerStyle(0)
hMCErr.SetFillColor(ROOT.kBlack)
hMCErr.SetLineColor(0)
hMCErr.Draw("sameE2")

leg1.Draw()

pad2.cd()
pad2.SetTopMargin(0)
pad2.SetFillColor(0)
pad2.SetFillStyle(0)
ROOT.gStyle.SetOptStat(0)
totalMC   = copy.deepcopy(hMCErr)
totalData = copy.deepcopy(h_Data)
totalData_forErrors = copy.deepcopy(h_Data)
totalData.Divide(totalMC)

for bin in range(1,hMCErr_size+1):

    #Set MC error band to MC relative uncertainty
    if not totalMC.GetBinContent(bin) == 0:
        new_MC_BinError = totalMC.GetBinError(bin)/totalMC.GetBinContent(bin)
    else:
        new_MC_BinError = 0.

    #Set data/MC ratio points error bar to data relative uncertainty
    if not totalData_forErrors.GetBinContent(bin) == 0:
        new_Data_BinError = totalData_forErrors.GetBinError(bin)/totalData_forErrors.GetBinContent(bin)
    else:
        new_Data_BinError = 0.

    totalMC.SetBinError(bin,new_MC_BinError)
    totalMC.SetBinContent(bin,1.)
    totalData.SetBinError(bin,new_Data_BinError)

totalData.SetTitle("")
totalData.SetMarkerStyle(8)
totalData.SetMarkerColor(1)
totalData.SetLineColor(1)
totalData.GetYaxis().SetLabelSize(0.10)
totalData.GetYaxis().SetTitle("Reco/Gen")
totalData.GetYaxis().SetTitleSize(0.08)
totalData.GetYaxis().SetTitleOffset(0.5)
#totalData.GetYaxis().SetRangeUser(0.4,1.6)

totalMC.SetTitle("")
totalMC.SetFillStyle(3002)

totalData.GetXaxis().SetTitle(h_Data.GetXaxis().GetTitle())
totalData.GetXaxis().SetLabelSize(0.10)
totalData.GetXaxis().SetTitleSize(0.12)
totalData.GetXaxis().SetTitleOffset(1.0)
totalData.SetMaximum(1.5)
totalData.SetMinimum(0.5)


line_on_one = ROOT.TLine(h_MC.GetXaxis().GetXmin(),1.,h_MC.GetXaxis().GetXmax(),1.)
line_on_one.SetLineColor(38)

totalData.Draw("E1,X0")
totalMC.Draw("sameE2")
line_on_one.Draw("SAME")


canvas.SaveAs("/eos/user/g/gumoret/www/testDir/ratio_plot_GenVsRecoMesonPt.pdf")
canvas.SaveAs("/eos/user/g/gumoret/www/testDir/ratio_plot_GenVsRecoMesonPt.png")
