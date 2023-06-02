import ROOT
import argparse
import math
import numpy as np
import sys
import copy
import tdrstyle, CMS_lumi
from ROOT import gROOT, RooFit
import os, sys

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True) 

#Import the Doube Crystal Ball PDF -----------------------------------
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+")

# H i s t o s
# -------------------------------------------------------------
h_width_offline_Data = ROOT.TH1F("h_width_offline_Data","width",3,0.,3.) 
h_mean_offline_Data  = ROOT.TH1F("h_mean_offline_Data","mean",3,0.,3.) 
h_width_offline_MC   = ROOT.TH1F("h_width_offline_MC","width",3,0.,3.) 
h_mean_offline_MC    = ROOT.TH1F("h_mean_offline_MC","mean",3,0.,3.) 

h_width_triggered_Data = ROOT.TH1F("h_width_triggered_Data","width",3,0.,3.) 
h_mean_triggered_Data  = ROOT.TH1F("h_mean_triggered_Data","mean",3,0.,3.) 
h_width_triggered_MC   = ROOT.TH1F("h_width_triggered_MC","width",3,0.,3.) 
h_mean_triggered_MC    = ROOT.TH1F("h_mean_triggered_MC","mean",3,0.,3.) 

h_efficiency_Data = ROOT.TH1F("h_efficiency_Data","eff data",3,0.,3.)
h_efficiency_MC   = ROOT.TH1F("h_efficiency_MC","eff MC",3,0.,3.)


# B i n   L o o p
# -------------------------------------------------------------
binList = [0,1,2]

for SAMPLE in ["Data", "MC"]:

	for binIndex in binList:

		# I n p u t   f i l e s
		# -------------------------------------------------------------
		fileInput = ROOT.TFile("histos/histos_trigger_efficiency/histos_TwoProngs/histos_trigger_"+SAMPLE+".root")
		fileInput.cd()
		tree_den = fileInput.Get("tree_output_bin"+str(binIndex)+"_den")
		tree_num = fileInput.Get("tree_output_bin"+str(binIndex)+"_num")


		# C r e a t e   m o d e l  
		# -------------------------------------------------------------

		# Create observables
		mass = ROOT.RooRealVar("MesonMass","ditrack invariant mass",1.,1.05)

		# Construct signal pdf
		mean_central ,  mean_min,  mean_max = 1.02,1.,1.04
		sigma_central, sigma_min, sigma_max = 0.0012,0.0001,0.01
		alpha_central, alpha_min, alpha_max = 1.2,-5.,10.
		n_central    ,     n_min,     n_max = 4.,0.,20.

		#if double sided CB
		alphaL_central, alphaL_min, alphaL_max = 1.2,0.,15.
		alphaR_central, alphaR_min, alphaR_max = 1.2,0.,15.
		nL_central    ,     nL_min,     nL_max = 1.,0.,5.
		nR_central    ,     nR_min,     nR_max = 1.,0.,5.

		mean  = ROOT.RooRealVar("mean","The mean of the Crystal Ball", mean_central, mean_min, mean_max)
		sigma = ROOT.RooRealVar("sigma","The width of the Crystal Ball", sigma_central, sigma_min, sigma_max)
		alpha = ROOT.RooRealVar("alpha","The alpha of the Crystal Ball", alpha_central, alpha_min, alpha_max)
		n     = ROOT.RooRealVar("n","The n of the Crystal Ball", n_central, n_min, n_max)

		#if double sided CB
		alphaL = ROOT.RooRealVar("alphaL","The alpha left of the double Crystal Ball", alphaL_central, alphaL_min, alphaL_max)
		alphaR = ROOT.RooRealVar("alphaR","The alpha right of the double Crystal Ball", alphaR_central, alphaR_min, alphaR_max)
		nL     = ROOT.RooRealVar("nL","The n left of the double Crystal Ball", nL_central, nL_min, nL_max)
		nR     = ROOT.RooRealVar("nR","The n right of the double Crystal Ball", nR_central, nR_min, nR_max)
		
		signalPDF = ROOT.RooCBShape("signalPDF","CB pdf",mass,mean,sigma,alpha,n)
		#signalPDF = ROOT.RooGaussian("signalPDF","CB pdf",mass,mean,sigma)
		#signalPDF = ROOT.RooDoubleCBFast("signalPDF","CB pdf",mass,mean,sigma,alphaL,alphaR,nL,nR)


		# Construct background pdf
		a1_central, a1_min, a1_max = 0.1,-0.5,0.5
		a2_central, a2_min, a2_max = 0.3,-2.,2.
		a3_central, a3_min, a3_max = 0.3,-2.,2.

		a1_off = ROOT.RooRealVar("a1_off","The a1 of background", a1_central, a1_min, a1_max)
		a2_off = ROOT.RooRealVar("a2_off","The a2 of background", a2_central, a2_min, a2_max)
		a3_off = ROOT.RooRealVar("a3_off","The a3 of background", a3_central, a3_min, a3_max)

		a1_trg = ROOT.RooRealVar("a1_trg","The a1 of background", a1_central, a1_min, a1_max)
		a2_trg = ROOT.RooRealVar("a2_trg","The a2 of background", a2_central, a2_min, a2_max)
		a3_trg = ROOT.RooRealVar("a3_trg","The a3 of background", a3_central, a3_min, a3_max)

		backgroundPDFOffline   = ROOT.RooChebychev("backgroundPDFOffline","The background PDF for offline dataset",mass,ROOT.RooArgList(a1_off,a2_off))#,a3_off))
		backgroundPDFTriggered = ROOT.RooChebychev("backgroundPDFTriggered","The background PDF for triggered dataset",mass,ROOT.RooArgList(a1_off,a2_off))#,a3_off))

		# N   e v e n t s
		# ---------------------------------------------------------------

		if SAMPLE == "Data": nsig = ROOT.RooRealVar("nsig", "signal yield 1", 700., 300., 1200.)
		if SAMPLE == "MC":   nsig = ROOT.RooRealVar("nsig", "signal yield 1", 150., 50., 300.)

		nbkgOffline   = ROOT.RooRealVar("nbkgOffline", "background yield", 3000., 500., 10000.)
		nbkgTriggered = ROOT.RooRealVar("nbkgTriggered", "background yield", 3000., 500., 10000.)
		efficiency    = ROOT.RooRealVar("efficiency", "efficiency", 0.95, 0.9, 1.)

		# R e t r i e v e   d a t a s e t s
		# ---------------------------------------------------------------

		data1 = ROOT.RooDataSet("datasetOffline","datasetOffline",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree_den))
		data2 = ROOT.RooDataSet("datasetTriggered","datasetTriggered",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree_num))

		# C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s
		# ---------------------------------------------------------------------------

		# definisci le formule per i signal yield nei due dataset
		signal_yield1_formula = ROOT.RooFormulaVar("signal_yield1_formula", "@0", ROOT.RooArgList(nsig))
		signal_yield2_formula = ROOT.RooFormulaVar("signal_yield2_formula", "@0 * @1", ROOT.RooArgList(nsig, efficiency))

		# Create index category and join samples
		sample = ROOT.RooCategory("sample", "sample")
		sample.defineType("Offline")
		sample.defineType("Triggered")

		#allData = ROOT.RooFit.RooDataset("allData","all data", ROOT.RooArgSet(mass), ROOT.RooFit.Import("Offline", data1))

		jointData = ROOT.RooDataSet("jointData", "joint data", ROOT.RooArgSet(mass), ROOT.RooFit.Index(sample), ROOT.RooFit.Import("Offline", data1), ROOT.RooFit.Import("Triggered", data2))

		# C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
		# -----------------------------------------------------------------------------------
		simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)

		# Create model for each dataset
		modelOffline   = ROOT.RooAddPdf("modelOffline", "Signal and Background PDF for Offline", ROOT.RooArgList(signalPDF, backgroundPDFOffline), ROOT.RooArgList(signal_yield1_formula, nbkgOffline))
		modelTriggered = ROOT.RooAddPdf("modelTriggered", "Signal and Background PDF for Triggered", ROOT.RooArgList(signalPDF, backgroundPDFTriggered), ROOT.RooArgList(signal_yield2_formula, nbkgTriggered))

		simPdf.addPdf(modelOffline, "Offline")
		simPdf.addPdf(modelTriggered, "Triggered")
		
		'''
		allData = ROOT.RooDataSet("allData","all data",ROOT.RooArgSet(mass),ROOT.RooFit.Index(sample), ROOT.RooFit.Import("Offline", data1))
		modelOffline.fitTo(allData)

		mean.setConstant(1)
		sigma.setConstant(1)
		alpha.setConstant(1)
		n.setConstant(1)
		a1_off.setConstant(1)
		a2_off.setConstant(1)
		'''
		# P e r f o r m   a   s i m u l t a n e o u s   f i t
		# ---------------------------------------------------
		fitResult = simPdf.fitTo(jointData, ROOT.RooFit.Extended(True), ROOT.RooFit.Save())

		# Print results
		fitResult.Print()


		# P l o t s
		# ---------------------------------------------------

		#CMS-style plotting 
		tdrstyle.setTDRStyle()
		iPeriod = 4
		iPos = 11
		CMS_lumi.lumiTextSize = 0.7
		CMS_lumi.lumiTextOffset = 0.25
		CMS_lumi.cmsTextSize = 0.8
		#CMS_lumi.cmsTextOffset = 0.4
		CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

		# Plot results
		canvas = ROOT.TCanvas("canvas", "canvas", 800, 400)
		canvas.Divide(2, 1)

		#FRAME1
		frame1 = mass.frame(30)
		frame1.SetTitle("Offline selection")
		jointData.plotOn(frame1, ROOT.RooFit.Cut("sample==sample::Offline"), ROOT.RooFit.Name("data1"))
		simPdf.plotOn(frame1, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Offline"), ROOT.RooFit.Components("backgroundPDFOffline"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("curve1_bkg"))
		simPdf.plotOn(frame1, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Offline"), ROOT.RooFit.Name("curve1_total"))
		frame1.GetXaxis().SetLabelSize(0.033)
		frame1.GetXaxis().SetTitle("m_{KK} [Gev]")
		frame1.SetMaximum(1.8 * frame1.GetMaximum())

		#FRAME2
		frame2 = mass.frame(30)
		frame2.SetTitle("Offline selection and trigger")
		jointData.plotOn(frame2, ROOT.RooFit.Cut("sample==sample::Triggered"), ROOT.RooFit.Name("data2"))
		simPdf.plotOn(frame2, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Triggered"), ROOT.RooFit.Components("backgroundPDFTriggered"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("curve2_bkg"))
		simPdf.plotOn(frame2, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Triggered"), ROOT.RooFit.Name("curve2_total"))
		frame2.GetXaxis().SetLabelSize(0.033)
		frame2.GetXaxis().SetTitle("m_{KK} [Gev]")
		frame2.SetMaximum(1.8 * frame2.GetMaximum())

		# LEGEND
		# Get the number of floating parameters for each fit
		nFitParams1 = simPdf.getParameters(data1).selectByAttrib("Constant", ROOT.kFALSE).getSize()
		nFitParams2 = simPdf.getParameters(data2).selectByAttrib("Constant", ROOT.kFALSE).getSize()

		# Calculate reduced chi-square values
		chi2_ndf1 = frame1.chiSquare(nFitParams1)
		chi2_ndf2 = frame2.chiSquare(nFitParams2)

		nsig1     = signal_yield1_formula.getParameter(0).getValV()
		nsig1_err = signal_yield1_formula.getParameter(0).errorVar().getValV()

		eff       = signal_yield2_formula.getParameter(1).getValV()
		eff_err   = signal_yield2_formula.getParameter(1).errorVar().getValV()
		nsig2     = signal_yield2_formula.getParameter(0).getValV()
		nsig2_err = signal_yield2_formula.getParameter(0).errorVar().getValV()

		# Get the fitted values of the Crystal Ball width (sigma) for each fit
		sigma_val = sigma.getVal()
		sigma_err = sigma.getError()
		mean_val  = mean.getVal()
		mean_err  = mean.getError()

		#Build legends
		x0, y0, x1, y1 = 0.4, 0.6, 0.91, 0.91 #legend size
		legend1 = ROOT.TLegend(x0, y0, x1, y1)
		legend1.SetBorderSize(0)
		legend1.SetFillStyle(0)
		legend1.SetTextSize(0.04)
		legend1.AddEntry(0, "Offline", "")
		legend1.AddEntry(frame1.findObject("curve1_total"), "Fit", "l")
		legend1.AddEntry(frame1.findObject("curve1_bkg"), "Bkg-only", "l")
		legend1.AddEntry(0, "#chi^{{2}}/ndf = {:.2f}".format(chi2_ndf1), "")
		legend1.AddEntry(0,"N_{sig} = "+str(round(nsig1,1))+" #pm "+str(round(nsig1_err,1)),"brNDC")
		legend1.AddEntry(0,"N_{bkg} = "+str(round(nbkgOffline.getValV(),1))+" #pm "+str(round(nbkgOffline.errorVar().getValV(),1)),"brNDC")
		#legend1.AddEntry(0,"#epsilon = "+str(round(eff,3))+" #pm "+str(round(eff_err,3)),"brNDC")
		legend1.AddEntry(0,"width = "+str(round(sigma_val,5))+" #pm "+str(round(sigma_err,5)),"brNDC")

		legend2 = ROOT.TLegend(x0, y0, x1, y1)
		legend2.SetBorderSize(0)
		legend2.SetFillStyle(0)
		legend2.SetTextSize(0.04)
		legend2.AddEntry(0, "Triggered", "")
		legend2.AddEntry(frame2.findObject("curve2_total"), "Fit", "l")
		legend2.AddEntry(frame2.findObject("curve2_bkg"), "Bkg-only", "l")
		legend2.AddEntry(0, "#chi^{{2}}/ndf = {:.2f}".format(chi2_ndf2), "")
		legend2.AddEntry(0,"N_{sig} = "+str(round(nsig2,1))+" #pm "+str(round(nsig2_err,1)),"brNDC")
		legend2.AddEntry(0,"N_{bkg} = "+str(round(nbkgTriggered.getValV(),1))+" #pm "+str(round(nbkgTriggered.errorVar().getValV(),1)),"brNDC")
		legend2.AddEntry(0,"#epsilon = "+str(round(eff,3))+" #pm "+str(round(eff_err,3)),"brNDC")
		legend2.AddEntry(0,"width = "+str(round(sigma_val,5))+" #pm "+str(round(sigma_err,5)),"brNDC")

		canvas.cd(1)
		frame1.Draw()
		legend1.Draw()
		canvas.cd(2)
		frame2.Draw()
		legend2.Draw()

		print "---------------------"
		print "nsig1     = ",nsig1," +/- ",
		print "nsig1_err = ",nsig1_err
		print "nsig2     = ",nsig2
		print "nsig2_err = ",nsig2_err
		print "eff       = ",eff
		print "eff_err   = ",eff_err
		print "---------------------"

		canvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/simultaneous_fit_bin_"+str(binIndex)+"_"+SAMPLE+".png")
		canvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/simultaneous_fit_bin_"+str(binIndex)+"_"+SAMPLE+".pdf")

		#Fill width and mean histos
		if SAMPLE == "Data":
			h_width_offline_Data.Fill(binIndex + 0.5,sigma_val)
			h_width_offline_Data.SetBinError(binIndex + 1,sigma_err)
			h_mean_offline_Data.Fill(binIndex + 0.5,mean_val)
			h_mean_offline_Data.SetBinError(binIndex + 1,mean_err)

			h_width_triggered_Data.Fill(binIndex + 0.5,sigma_val)
			h_width_triggered_Data.SetBinError(binIndex + 1,sigma_err)
			h_mean_triggered_Data.Fill(binIndex + 0.5,mean_val)
			h_mean_triggered_Data.SetBinError(binIndex + 1,mean_err)

			h_efficiency_Data.Fill(binIndex + 0.5,eff)
			h_efficiency_Data.SetBinError(binIndex + 1,eff_err)

		else:
			h_width_offline_MC.Fill(binIndex + 0.5,sigma_val)
			h_width_offline_MC.SetBinError(binIndex + 1,sigma_err)
			h_mean_offline_MC.Fill(binIndex + 0.5,mean_val)
			h_mean_offline_MC.SetBinError(binIndex + 1,mean_err)

			h_width_triggered_MC.Fill(binIndex + 0.5,sigma_val)
			h_width_triggered_MC.SetBinError(binIndex + 1,sigma_err)
			h_mean_triggered_MC.Fill(binIndex + 0.5,mean_val)
			h_mean_triggered_MC.SetBinError(binIndex + 1,mean_err)

			h_efficiency_MC.Fill(binIndex + 0.5,eff)
			h_efficiency_MC.SetBinError(binIndex + 1,eff_err)

# P l o t s
# ---------------------------------------------------

#width offline 
widthOfflineCanvas = ROOT.TCanvas("c1", "c1", 1000, 1000)
widthOfflineCanvas.cd()

leg1 = ROOT.TLegend(0.2,0.75,0.58,0.92) #left positioning
leg1.SetHeader(" ")
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.SetNColumns(1)
leg1.AddEntry(0,"offline dataset","")
leg1.AddEntry(h_width_offline_Data,"Data","ep")
leg1.AddEntry(h_width_offline_MC,"MC","ep")

h_width_offline_MC.SetMarkerStyle(20)
h_width_offline_MC.SetMarkerSize(2.)
h_width_offline_MC.SetMarkerColor(ROOT.kRed)
h_width_offline_MC.SetLineColor(ROOT.kRed)
h_width_offline_MC.SetAxisRange(0.0005,0.007,"Y")
h_width_offline_MC.GetXaxis().SetBinLabel(1,"bin0")
h_width_offline_MC.GetXaxis().SetBinLabel(2,"bin1")
h_width_offline_MC.GetXaxis().SetBinLabel(3,"bin2")
h_width_offline_MC.GetYaxis().SetTitle("width")
h_width_offline_MC.Draw()
h_width_offline_Data.SetMarkerStyle(20)
h_width_offline_Data.SetMarkerSize(2.)
h_width_offline_Data.SetLineColor(ROOT.kBlack)
h_width_offline_Data.Draw("SAME")

leg1.Draw()

#mean offline 
meanOfflineCanvas = ROOT.TCanvas("c2", "c2", 1000, 1000)
meanOfflineCanvas.cd()

leg2 = ROOT.TLegend(0.2,0.75,0.58,0.92) #left positioning
leg2.SetHeader(" ")
leg2.SetFillColorAlpha(0,0.)
leg2.SetBorderSize(0)
leg2.SetLineColor(1)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillStyle(1001)
leg2.SetNColumns(1)
leg2.AddEntry(0,"offline dataset","")
leg2.AddEntry(h_mean_offline_Data,"Data","ep")
leg2.AddEntry(h_mean_offline_MC,"MC","ep")

h_mean_offline_MC.SetMarkerStyle(20)
h_mean_offline_MC.SetMarkerSize(2.)
h_mean_offline_MC.SetMarkerColor(ROOT.kRed)
h_mean_offline_MC.SetLineColor(ROOT.kRed)
h_mean_offline_MC.SetAxisRange(1.016,1.025,"Y")
h_mean_offline_MC.GetXaxis().SetBinLabel(1,"bin0")
h_mean_offline_MC.GetXaxis().SetBinLabel(2,"bin1")
h_mean_offline_MC.GetXaxis().SetBinLabel(3,"bin2")
h_mean_offline_MC.GetYaxis().SetTitle("mean")
h_mean_offline_MC.Draw()
h_mean_offline_Data.SetMarkerStyle(20)
h_mean_offline_Data.SetMarkerSize(2.)
h_mean_offline_Data.SetLineColor(ROOT.kBlack)
h_mean_offline_Data.Draw("SAME")

leg2.Draw()

#width triggered 
widthTriggeredCanvas = ROOT.TCanvas("c3", "c3", 1000, 1000)
widthTriggeredCanvas.cd()

leg3 = ROOT.TLegend(0.2,0.75,0.58,0.92) #left positioning
leg3.SetHeader(" ")
leg3.SetFillColorAlpha(0,0.)
leg3.SetBorderSize(0)
leg3.SetLineColor(1)
leg3.SetLineStyle(1)
leg3.SetLineWidth(1)
leg3.SetFillStyle(1001)
leg3.SetNColumns(1)
leg3.AddEntry(0,"triggered dataset","")
leg3.AddEntry(h_width_triggered_Data,"Data","ep")
leg3.AddEntry(h_width_triggered_MC,"MC","ep")

h_width_triggered_MC.SetMarkerStyle(20)
h_width_triggered_MC.SetMarkerSize(2.)
h_width_triggered_MC.SetMarkerColor(ROOT.kRed)
h_width_triggered_MC.SetLineColor(ROOT.kRed)
h_width_triggered_MC.SetAxisRange(0.000,0.007,"Y")
h_width_triggered_MC.GetXaxis().SetBinLabel(1,"bin0")
h_width_triggered_MC.GetXaxis().SetBinLabel(2,"bin1")
h_width_triggered_MC.GetXaxis().SetBinLabel(3,"bin2")
h_width_triggered_MC.GetYaxis().SetTitle("width")
h_width_triggered_MC.Draw()
h_width_triggered_Data.SetMarkerStyle(20)
h_width_triggered_Data.SetMarkerSize(2.)
h_width_triggered_Data.SetLineColor(ROOT.kBlack)
h_width_triggered_Data.Draw("SAME")

leg3.Draw()

#mean triggered 
meanTriggeredCanvas = ROOT.TCanvas("c4", "c4", 1000, 1000)
meanTriggeredCanvas.cd()

leg4 = ROOT.TLegend(0.2,0.75,0.58,0.92) #left positioning
leg4.SetHeader(" ")
leg4.SetFillColorAlpha(0,0.)
leg4.SetBorderSize(0)
leg4.SetLineColor(1)
leg4.SetLineStyle(1)
leg4.SetLineWidth(1)
leg4.SetFillStyle(1001)
leg4.SetNColumns(1)
leg4.AddEntry(0,"triggered dataset","")
leg4.AddEntry(h_mean_triggered_Data,"Data","ep")
leg4.AddEntry(h_mean_triggered_MC,"MC","ep")

h_mean_triggered_MC.SetMarkerStyle(20)
h_mean_triggered_MC.SetMarkerSize(2.)
h_mean_triggered_MC.SetMarkerColor(ROOT.kRed)
h_mean_triggered_MC.SetLineColor(ROOT.kRed)
h_mean_triggered_MC.SetAxisRange(1.016,1.025,"Y")
h_mean_triggered_MC.GetXaxis().SetBinLabel(1,"bin0")
h_mean_triggered_MC.GetXaxis().SetBinLabel(2,"bin1")
h_mean_triggered_MC.GetXaxis().SetBinLabel(3,"bin2")
h_mean_triggered_MC.GetYaxis().SetTitle("mean")
h_mean_triggered_MC.Draw()
h_mean_triggered_Data.SetMarkerStyle(20)
h_mean_triggered_Data.SetMarkerSize(2.)
h_mean_triggered_Data.SetLineColor(ROOT.kBlack)
h_mean_triggered_Data.Draw("SAME")

leg4.Draw()


widthOfflineCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/width.png")
meanOfflineCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/mean.png")
#widthTriggeredCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/widthTriggered.png")
#meanTriggeredCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/meanTriggered.png")
widthOfflineCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/width.pdf")
meanOfflineCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/mean.pdf")
#widthTriggeredCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/widthTriggered.pdf")
#meanTriggeredCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/meanTriggered.pdf")

fOut = ROOT.TFile("histos/histos_trigger_efficiency/histos_TwoProngs/histos_TwoProngsSimFit.root","RECREATE")
fOut.cd()

h_efficiency_MC.Write()
h_efficiency_Data.Write()

fOut.Close()



