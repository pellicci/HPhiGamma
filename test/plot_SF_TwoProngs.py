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
inputFile    = ROOT.TFile("histos/histos_trigger_efficiency/histos_TwoProngs/histos_TwoProngsSimFit.root")

h_data   = inputFile.Get("h_efficiency_Data")
h_mc     = inputFile.Get("h_efficiency_MC")

print("Number of bins in h_data    : ", h_data.GetNbinsX())
print("Number of bins in h_mc      : ", h_mc.GetNbinsX())

def create_ratio_plot(h_data, h_mc):
    canvas = ROOT.TCanvas("canvas", "Efficiency Comparison", 800, 800)

    # Upper pad
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
    ROOT.gStyle.SetOptStat(0)
    pad1.SetBottomMargin(0.03)
    pad1.Draw()
    pad1.cd()

    # Remove numbers from x-axis
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetYaxis().SetTitle("Efficiency")

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

    # Create legend
    legend = ROOT.TLegend(0.2, 0.2, 0.4, 0.32)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextSize(0.04)
    legend.AddEntry(h_data, "Data", "lep")
    legend.AddEntry(h_mc, "MC", "lep")
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
    canvas.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/TwoProngsSF.png")
    canvas.SaveAs("/eos/user/g/gumoret/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/TwoProngsSF.pdf")

create_ratio_plot(h_data, h_mc)
