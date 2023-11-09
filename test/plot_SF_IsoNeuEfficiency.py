import ROOT
import argparse
import math
import numpy as np
import sys
import copy
import tdrstyle
import CMS_lumi
from ROOT import gROOT
import os, sys

# Supress the opening of many Canvases
ROOT.gROOT.SetBatch(True)

p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('eta', help='Type the eta range')
args = p.parse_args()
ETA = args.eta

# CMS-style plotting
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
# CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}"

# Input
# --------------------------------------------------------------------------------------
inputFile    = ROOT.TFile("histos/IsoEfficiency/histos_IsoNeuEfficiencySimFit_"+ETA+".root")
inputFilePhi = ROOT.TFile("histos/IsoEfficiency/histos_IsoNeuEfficiencySimFit_Phi_"+ETA+".root")
inputFileRho = ROOT.TFile("histos/IsoEfficiency/histos_IsoNeuEfficiencySimFit_Rho_"+ETA+".root")
inputFileK0s = ROOT.TFile("histos/IsoEfficiency/histos_IsoNeuEfficiencySimFit_K0s_"+ETA+".root")

h_data   = inputFile.Get("h_efficiency_Data")
h_mc     = inputFile.Get("h_efficiency_MC")
h_phi    = inputFilePhi.Get("h_efficiency_Signal")
h_rho    = inputFileRho.Get("h_efficiency_Signal")
h_k0s    = inputFileK0s.Get("h_efficiency_Signal")

print("Number of bins in h_data    : ", h_data.GetNbinsX())
print("Number of bins in h_mc      : ", h_mc.GetNbinsX())
print("Number of bins in h_phi     : ", h_phi.GetNbinsX())
print("Number of bins in h_rho     : ", h_rho.GetNbinsX())

pT_list = [38, 40, 42, 44, 46, 50, 55]


def create_ratio_plot(h_data, h_mc, h_phi, h_rho, h_k0s):
    canvas = ROOT.TCanvas("canvas", "Efficiency Comparison", 800, 800)

    # Upper pad
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
    ROOT.gStyle.SetOptStat(0)
    pad1.SetBottomMargin(0.03)
    pad1.Draw()
    pad1.cd()

    # Remove numbers from x-axis
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetYaxis().SetTitle("IsoNeu Efficiency")

    # Draw data histogram
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetMarkerSize(1.)
    h_data.Draw("E")

    # Draw MC histogram
    h_mc.SetLineColor(ROOT.kBlue)
    h_mc.SetMarkerStyle(20)
    h_mc.SetMarkerColor(ROOT.kBlue)
    h_mc.SetMarkerSize(1.)
    h_mc.Draw("E SAME")

    # Draw Signal histogram
    h_phi.SetLineColor(ROOT.kCyan)
    h_phi.SetMarkerStyle(20)
    h_phi.SetMarkerColor(ROOT.kCyan)
    h_phi.SetMarkerSize(1.)
    #h_phi.Draw("E SAME")

    h_rho.SetLineColor(ROOT.kRed)
    h_rho.SetMarkerStyle(20)
    h_rho.SetMarkerColor(ROOT.kRed)
    h_rho.SetMarkerSize(1.)
    #h_rho.Draw("E SAME")

    h_k0s.SetLineColor(6)
    h_k0s.SetMarkerStyle(20)
    h_k0s.SetMarkerColor(6)
    h_k0s.SetMarkerSize(1.)
    #h_k0s.Draw("E SAME")

    # Create legend
    legend = ROOT.TLegend(0.2, 0.1, 0.4, 0.5)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextSize(0.04)
    legend.AddEntry(0,"IsoNeu efficiency","")
    legend.AddEntry(h_data, "Data", "lep")
    legend.AddEntry(h_mc, "MC", "lep")
    #legend.AddEntry(h_phi, "#phi#gamma", "lep")
    #legend.AddEntry(h_rho, "#rho#gamma", "lep")
    #legend.AddEntry(h_k0s, "K_{0}^{*}#gamma", "lep")
    legend.Draw()

    # Grid lines on the upper pad
    ROOT.gStyle.SetGridStyle(3)
    pad1.SetGrid()

    # Set y-axis range for the upper pad
    h_data.GetYaxis().SetRangeUser(0.4, 1.2)

    # Lower pad (ratio plot)
    canvas.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()

    # Calculate ratio
    h_ratio = h_data.Clone()
    h_ratio.Divide(h_mc)

    # Set ratio plot options
    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetMarkerColor(ROOT.kBlack)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(1.) 
    h_ratio.SetMinimum(0.5)
    h_ratio.SetMaximum(1.5)
    h_ratio.SetTitle("")
    h_ratio.GetYaxis().SetTitle("Data/MC SF")
    h_ratio.GetYaxis().SetTitleSize(0.1)
    h_ratio.GetYaxis().SetTitleOffset(0.3)
    h_ratio.GetYaxis().SetLabelSize(0.08)
    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetTitle("p_{T} [GeV]")    
    #h_ratio.GetXaxis().SetTitleOffset(1.0)
    h_ratio.GetXaxis().SetLabelSize(0.1)
    h_ratio.GetXaxis().SetTickLength(0.07)
    h_ratio.Draw("E")

    # Set the x-axis range for the ratio plot
    h_ratio.GetXaxis().SetRangeUser(38, 80)

    # Set the number of divisions on the x-axis
    h_ratio.GetXaxis().SetNdivisions(420)

    # Set the x-axis labels
    h_ratio.GetXaxis().SetLabelSize(0.1)
    h_ratio.GetXaxis().SetLabelOffset(0.02)

    # Draw the x-axis labels
    h_ratio.GetXaxis().Draw()

    # Grid lines on the upper pad
    ROOT.gStyle.SetGridStyle(3)
    pad1.SetGrid()

    # Draw horizontal line at ratio 1.0
    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.Draw()

    canvas.Update()
    canvas.SaveAs("~/cernbox/www/latest_production/IsoNeuEfficiency_Fit_"+ETA+"_latest_production/IsoNeuEfficiencySF_"+ETA+".png")
    canvas.SaveAs("~/cernbox/www/latest_production/IsoNeuEfficiency_Fit_"+ETA+"_latest_production/IsoNeuEfficiencySF_"+ETA+".pdf")


create_ratio_plot(h_data, h_mc, h_phi, h_rho, h_k0s)

fOut = ROOT.TFile("scale_factors/IsoNeuEfficiencySF_"+ETA+".root","RECREATE")
fOut.cd()

h_data.Divide(h_mc)
h_data.Write()

fOut.Close()