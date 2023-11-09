import ROOT
import argparse
import math
import numpy as np
import sys
import copy
import tdrstyle, CMS_lumi
from array import array
from ROOT import gROOT, RooFit
import os, sys

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True) 

#Import the Doube Crystal Ball PDF -----------------------------------
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+")

# B i n   L o o p
# -------------------------------------------------------------
binList = [0,1,2]
pT_list = [38.,42.,50.,80.] #this list is just for the plot representation

# H i s t o s
# -------------------------------------------------------------
h_width_Fail_Data = ROOT.TH1F("h_width_Fail_Data","width",3,0.,3.) 
h_mean_Fail_Data  = ROOT.TH1F("h_mean_Fail_Data","mean",3,0.,3.) 
h_width_Fail_MC   = ROOT.TH1F("h_width_Fail_MC","width",3,0.,3.) 
h_mean_Fail_MC    = ROOT.TH1F("h_mean_Fail_MC","mean",3,0.,3.) 

h_width_Pass_Data = ROOT.TH1F("h_width_Pass_Data","width",3,0.,3.) 
h_mean_Pass_Data  = ROOT.TH1F("h_mean_Pass_Data","mean",3,0.,3.) 
h_width_Pass_MC   = ROOT.TH1F("h_width_Pass_MC","width",3,0.,3.) 
h_mean_Pass_MC    = ROOT.TH1F("h_mean_Pass_MC","mean",3,0.,3.) 

h_efficiency_Data = ROOT.TH1F("h_efficiency_Data","eff data",len(pT_list) - 1, array('d', pT_list))
h_efficiency_MC   = ROOT.TH1F("h_efficiency_MC","eff MC",len(pT_list) - 1, array('d', pT_list))
h_TwoProngsTrigger_SF = ROOT.TH1F("h_TwoProngsTriggerSF","TwoProngs trigger leg SFs",len(pT_list) - 1, array('d', pT_list))

for SAMPLE in ["Data", "MC"]:

	for binIndex in binList:

		# I n p u t   f i l e s
		# -------------------------------------------------------------
		fileInput = ROOT.TFile("histos/histos_trigger_efficiency/histos_TwoProngs/histos_trigger_"+SAMPLE+".root")
		fileInput.cd()
		tree_fail = fileInput.Get("tree_output_bin"+str(binIndex)+"_fail")
		tree_pass = fileInput.Get("tree_output_bin"+str(binIndex)+"_pass")

		nEntriesFail = tree_fail.GetEntriesFast()
		nEntriesPass = tree_pass.GetEntriesFast()

		# C r e a t e   m o d e l  
		# -------------------------------------------------------------

		# Create observables
		mass = ROOT.RooRealVar("MesonMass","ditrack invariant mass",1.,1.05)

		# Construct signal pdf
		mean_central ,  mean_min,  mean_max = 1.02,1.,1.04
		sigma_central, sigma_min, sigma_max = 0.0022,0.00001,0.003
		alpha_central, alpha_min, alpha_max = 1.2,-5.,10.
		n_central    ,     n_min,     n_max = 4.,0.,20.

		#if double sided CB
		alphaL_central, alphaL_min, alphaL_max = 1.2,-5.,10.
		alphaR_central, alphaR_min, alphaR_max = 1.2,-5.,10.
		nL_central    ,     nL_min,     nL_max = 4.,0.,20.
		nR_central    ,     nR_min,     nR_max = 4.,0.,20.

		mean  = ROOT.RooRealVar("mean","The mean of the Crystal Ball", mean_central, mean_min, mean_max)
		sigma = ROOT.RooRealVar("sigma","The width of the Crystal Ball", sigma_central, sigma_min, sigma_max)
		alpha = ROOT.RooRealVar("alpha","The alpha of the Crystal Ball", alpha_central, alpha_min, alpha_max)
		n     = ROOT.RooRealVar("n","The n of the Crystal Ball", n_central, n_min, n_max)

		#if double sided CB
		alphaL = ROOT.RooRealVar("alphaL","The alpha left of the double Crystal Ball", alphaL_central, alphaL_min, alphaL_max)
		alphaR = ROOT.RooRealVar("alphaR","The alpha right of the double Crystal Ball", alphaR_central, alphaR_min, alphaR_max)
		nL     = ROOT.RooRealVar("nL","The n left of the double Crystal Ball", nL_central, nL_min, nL_max)
		nR     = ROOT.RooRealVar("nR","The n right of the double Crystal Ball", nR_central, nR_min, nR_max)
		
		#signalPDF = ROOT.RooCBShape("signalPDF","CB pdf",mass,mean,sigma,alpha,n)
		signalPDF = ROOT.RooGaussian("signalPDF","CB pdf",mass,mean,sigma)
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

		backgroundPDFFail = ROOT.RooChebychev("backgroundPDFFail","The background PDF for Fail dataset",mass,ROOT.RooArgList(a1_off,a2_off,a3_off))
		backgroundPDFPass = ROOT.RooChebychev("backgroundPDFPass","The background PDF for Pass dataset",mass,ROOT.RooArgList(a1_off,a2_off,a3_off))

		#backgroundPDFFail = ROOT.RooChebychev("backgroundPDFFail","The background PDF for Fail dataset",mass,ROOT.RooArgList(a1_off,a2_off))
		#backgroundPDFPass = ROOT.RooChebychev("backgroundPDFPass","The background PDF for Pass dataset",mass,ROOT.RooArgList(a1_off,a2_off))

		# N   e v e n t s
		# ---------------------------------------------------------------

		nsig = ROOT.RooRealVar("nsig", "signal yield 1", nEntriesFail, nEntriesFail/10, nEntriesPass*1.5)
		#if SAMPLE == "MC":   nsig = ROOT.RooRealVar("nsig", "signal yield 1", 150., 0., 300.)

		nbkgFail   = ROOT.RooRealVar("nbkgFail", "background yield", nEntriesFail, nEntriesFail/10, nEntriesFail*2)
		nbkgPass   = ROOT.RooRealVar("nbkgPass", "background yield", nEntriesPass, nEntriesPass/10, nEntriesPass*2)
		efficiency = ROOT.RooRealVar("efficiency", "efficiency", 0.95, 0., 1.)

		# R e t r i e v e   d a t a s e t s
		# ---------------------------------------------------------------

		data1 = ROOT.RooDataSet("datasetFail","datasetFail",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree_fail))
		data2 = ROOT.RooDataSet("datasetPass","datasetPass",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree_pass))

		# C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s
		# ---------------------------------------------------------------------------

		# definisci le formule per i signal yield nei due dataset
		signal_yield1_formula = ROOT.RooFormulaVar("signal_yield1_formula", "@0 * (1 - @1)", ROOT.RooArgList(nsig,efficiency))
		signal_yield2_formula = ROOT.RooFormulaVar("signal_yield2_formula", "@0 * @1", ROOT.RooArgList(nsig, efficiency))

		# Create index category and join samples
		sample = ROOT.RooCategory("sample", "sample")
		sample.defineType("Fail")
		sample.defineType("Pass")

		#allData = ROOT.RooFit.RooDataset("allData","all data", ROOT.RooArgSet(mass), ROOT.RooFit.Import("Fail", data1))

		jointData = ROOT.RooDataSet("jointData", "joint data", ROOT.RooArgSet(mass), ROOT.RooFit.Index(sample), ROOT.RooFit.Import("Pass", data2), ROOT.RooFit.Import("Fail", data1))

		# C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
		# -----------------------------------------------------------------------------------
		simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)

		# Create model for each dataset
		modelPass = ROOT.RooAddPdf("modelPass", "Signal and Background PDF for Pass", ROOT.RooArgList(signalPDF, backgroundPDFPass), ROOT.RooArgList(signal_yield2_formula, nbkgPass))
		modelFail = ROOT.RooAddPdf("modelFail", "Signal and Background PDF for Fail", ROOT.RooArgList(signalPDF, backgroundPDFFail), ROOT.RooArgList(signal_yield1_formula, nbkgFail))

		simPdf.addPdf(modelPass, "Pass")
		simPdf.addPdf(modelFail, "Fail")
		
		'''
		allData = ROOT.RooDataSet("allData","all data",ROOT.RooArgSet(mass),ROOT.RooFit.Index(sample), ROOT.RooFit.Import("Fail", data1))
		modelFail.fitTo(allData)

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
		frame1.SetTitle("Fail selection")
		jointData.plotOn(frame1, ROOT.RooFit.Cut("sample==sample::Fail"), ROOT.RooFit.Name("data1"))
		simPdf.plotOn(frame1, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Fail"), ROOT.RooFit.Components("backgroundPDFFail"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("curve1_bkg"))
		simPdf.plotOn(frame1, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Fail"), ROOT.RooFit.Name("curve1_total"))
		frame1.GetXaxis().SetLabelSize(0.033)
		frame1.GetXaxis().SetTitle("m_{KK} [Gev]")
		frame1.SetMaximum(1.8 * frame1.GetMaximum())

		#FRAME2
		frame2 = mass.frame(30)
		frame2.SetTitle("Fail selection and trigger")
		jointData.plotOn(frame2, ROOT.RooFit.Cut("sample==sample::Pass"), ROOT.RooFit.Name("data2"))
		simPdf.plotOn(frame2, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Pass"), ROOT.RooFit.Components("backgroundPDFPass"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("curve2_bkg"))
		simPdf.plotOn(frame2, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Pass"), ROOT.RooFit.Name("curve2_total"))
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
		legend1.AddEntry(0, "Fail", "")
		legend1.AddEntry(frame1.findObject("curve1_total"), "Fit", "l")
		legend1.AddEntry(frame1.findObject("curve1_bkg"), "Bkg-only", "l")
		legend1.AddEntry(0, "#chi^{{2}}/ndf = {:.2f}".format(chi2_ndf1), "")
		legend1.AddEntry(0,"N_{sig} = "+str(round(nsig1,1))+" #pm "+str(round(nsig1_err,1)),"brNDC")
		legend1.AddEntry(0,"N_{bkg} = "+str(round(nbkgFail.getValV(),1))+" #pm "+str(round(nbkgFail.errorVar().getValV(),1)),"brNDC")
		#legend1.AddEntry(0,"#epsilon = "+str(round(eff,3))+" #pm "+str(round(eff_err,3)),"brNDC")
		legend1.AddEntry(0,"width = "+str(round(sigma_val,5))+" #pm "+str(round(sigma_err,5)),"brNDC")

		legend2 = ROOT.TLegend(x0, y0, x1, y1)
		legend2.SetBorderSize(0)
		legend2.SetFillStyle(0)
		legend2.SetTextSize(0.04)
		legend2.AddEntry(0, "Pass", "")
		legend2.AddEntry(frame2.findObject("curve2_total"), "Fit", "l")
		legend2.AddEntry(frame2.findObject("curve2_bkg"), "Bkg-only", "l")
		legend2.AddEntry(0, "#chi^{{2}}/ndf = {:.2f}".format(chi2_ndf2), "")
		legend2.AddEntry(0,"N_{sig} = "+str(round(nsig2,1))+" #pm "+str(round(nsig2_err,1)),"brNDC")
		legend2.AddEntry(0,"N_{bkg} = "+str(round(nbkgPass.getValV(),1))+" #pm "+str(round(nbkgPass.errorVar().getValV(),1)),"brNDC")
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
			h_width_Fail_Data.Fill(binIndex + 0.5,sigma_val)
			h_width_Fail_Data.SetBinError(binIndex + 1,sigma_err)
			h_mean_Fail_Data.Fill(binIndex + 0.5,mean_val)
			h_mean_Fail_Data.SetBinError(binIndex + 1,mean_err)

			h_width_Pass_Data.Fill(binIndex + 0.5,sigma_val)
			h_width_Pass_Data.SetBinError(binIndex + 1,sigma_err)
			h_mean_Pass_Data.Fill(binIndex + 0.5,mean_val)
			h_mean_Pass_Data.SetBinError(binIndex + 1,mean_err)

			h_efficiency_Data.SetBinContent(binIndex + 1,eff)
			h_efficiency_Data.SetBinError(binIndex + 1,eff_err)

		else:
			h_width_Fail_MC.Fill(binIndex + 0.5,sigma_val)
			h_width_Fail_MC.SetBinError(binIndex + 1,sigma_err)
			h_mean_Fail_MC.Fill(binIndex + 0.5,mean_val)
			h_mean_Fail_MC.SetBinError(binIndex + 1,mean_err)

			h_width_Pass_MC.Fill(binIndex + 0.5,sigma_val)
			h_width_Pass_MC.SetBinError(binIndex + 1,sigma_err)
			h_mean_Pass_MC.Fill(binIndex + 0.5,mean_val)
			h_mean_Pass_MC.SetBinError(binIndex + 1,mean_err)

			h_efficiency_MC.SetBinContent(binIndex + 1,eff)
			h_efficiency_MC.SetBinError(binIndex + 1,eff_err)

# P l o t s
# ---------------------------------------------------

#width Fail 
widthFailCanvas = ROOT.TCanvas("c1", "c1", 1000, 1000)
widthFailCanvas.cd()

leg1 = ROOT.TLegend(0.2,0.75,0.58,0.92) #left positioning
leg1.SetHeader(" ")
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.SetNColumns(1)
leg1.AddEntry(0,"Fail dataset","")
leg1.AddEntry(h_width_Fail_Data,"Data","ep")
leg1.AddEntry(h_width_Fail_MC,"MC","ep")

h_width_Fail_MC.SetMarkerStyle(20)
h_width_Fail_MC.SetMarkerSize(2.)
h_width_Fail_MC.SetMarkerColor(ROOT.kRed)
h_width_Fail_MC.SetLineColor(ROOT.kRed)
h_width_Fail_MC.SetAxisRange(0.0005,0.007,"Y")
h_width_Fail_MC.GetXaxis().SetBinLabel(1,"bin0")
h_width_Fail_MC.GetXaxis().SetBinLabel(2,"bin1")
h_width_Fail_MC.GetXaxis().SetBinLabel(3,"bin2")
h_width_Fail_MC.GetYaxis().SetTitle("width")
h_width_Fail_MC.Draw()
h_width_Fail_Data.SetMarkerStyle(20)
h_width_Fail_Data.SetMarkerSize(2.)
h_width_Fail_Data.SetLineColor(ROOT.kBlack)
h_width_Fail_Data.Draw("SAME")

leg1.Draw()

#mean Fail 
meanFailCanvas = ROOT.TCanvas("c2", "c2", 1000, 1000)
meanFailCanvas.cd()

leg2 = ROOT.TLegend(0.2,0.75,0.58,0.92) #left positioning
leg2.SetHeader(" ")
leg2.SetFillColorAlpha(0,0.)
leg2.SetBorderSize(0)
leg2.SetLineColor(1)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillStyle(1001)
leg2.SetNColumns(1)
leg2.AddEntry(0,"Fail dataset","")
leg2.AddEntry(h_mean_Fail_Data,"Data","ep")
leg2.AddEntry(h_mean_Fail_MC,"MC","ep")

h_mean_Fail_MC.SetMarkerStyle(20)
h_mean_Fail_MC.SetMarkerSize(2.)
h_mean_Fail_MC.SetMarkerColor(ROOT.kRed)
h_mean_Fail_MC.SetLineColor(ROOT.kRed)
h_mean_Fail_MC.SetAxisRange(1.016,1.025,"Y")
h_mean_Fail_MC.GetXaxis().SetBinLabel(1,"bin0")
h_mean_Fail_MC.GetXaxis().SetBinLabel(2,"bin1")
h_mean_Fail_MC.GetXaxis().SetBinLabel(3,"bin2")
h_mean_Fail_MC.GetYaxis().SetTitle("mean")
h_mean_Fail_MC.Draw()
h_mean_Fail_Data.SetMarkerStyle(20)
h_mean_Fail_Data.SetMarkerSize(2.)
h_mean_Fail_Data.SetLineColor(ROOT.kBlack)
h_mean_Fail_Data.Draw("SAME")

leg2.Draw()

#width Pass 
widthPassCanvas = ROOT.TCanvas("c3", "c3", 1000, 1000)
widthPassCanvas.cd()

leg3 = ROOT.TLegend(0.2,0.75,0.58,0.92) #left positioning
leg3.SetHeader(" ")
leg3.SetFillColorAlpha(0,0.)
leg3.SetBorderSize(0)
leg3.SetLineColor(1)
leg3.SetLineStyle(1)
leg3.SetLineWidth(1)
leg3.SetFillStyle(1001)
leg3.SetNColumns(1)
leg3.AddEntry(0,"Pass dataset","")
leg3.AddEntry(h_width_Pass_Data,"Data","ep")
leg3.AddEntry(h_width_Pass_MC,"MC","ep")

h_width_Pass_MC.SetMarkerStyle(20)
h_width_Pass_MC.SetMarkerSize(2.)
h_width_Pass_MC.SetMarkerColor(ROOT.kRed)
h_width_Pass_MC.SetLineColor(ROOT.kRed)
h_width_Pass_MC.SetAxisRange(0.000,0.007,"Y")
h_width_Pass_MC.GetXaxis().SetBinLabel(1,"bin0")
h_width_Pass_MC.GetXaxis().SetBinLabel(2,"bin1")
h_width_Pass_MC.GetXaxis().SetBinLabel(3,"bin2")
h_width_Pass_MC.GetYaxis().SetTitle("width")
h_width_Pass_MC.Draw()
h_width_Pass_Data.SetMarkerStyle(20)
h_width_Pass_Data.SetMarkerSize(2.)
h_width_Pass_Data.SetLineColor(ROOT.kBlack)
h_width_Pass_Data.Draw("SAME")

leg3.Draw()

#mean Pass 
meanPassCanvas = ROOT.TCanvas("c4", "c4", 1000, 1000)
meanPassCanvas.cd()

leg4 = ROOT.TLegend(0.2,0.75,0.58,0.92) #left positioning
leg4.SetHeader(" ")
leg4.SetFillColorAlpha(0,0.)
leg4.SetBorderSize(0)
leg4.SetLineColor(1)
leg4.SetLineStyle(1)
leg4.SetLineWidth(1)
leg4.SetFillStyle(1001)
leg4.SetNColumns(1)
leg4.AddEntry(0,"Pass dataset","")
leg4.AddEntry(h_mean_Pass_Data,"Data","ep")
leg4.AddEntry(h_mean_Pass_MC,"MC","ep")

h_mean_Pass_MC.SetMarkerStyle(20)
h_mean_Pass_MC.SetMarkerSize(2.)
h_mean_Pass_MC.SetMarkerColor(ROOT.kRed)
h_mean_Pass_MC.SetLineColor(ROOT.kRed)
h_mean_Pass_MC.SetAxisRange(1.016,1.025,"Y")
h_mean_Pass_MC.GetXaxis().SetBinLabel(1,"bin0")
h_mean_Pass_MC.GetXaxis().SetBinLabel(2,"bin1")
h_mean_Pass_MC.GetXaxis().SetBinLabel(3,"bin2")
h_mean_Pass_MC.GetYaxis().SetTitle("mean")
h_mean_Pass_MC.Draw()
h_mean_Pass_Data.SetMarkerStyle(20)
h_mean_Pass_Data.SetMarkerSize(2.)
h_mean_Pass_Data.SetLineColor(ROOT.kBlack)
h_mean_Pass_Data.Draw("SAME")

leg4.Draw()


widthFailCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/width.png")
meanFailCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/mean.png")
#widthPassCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/widthPass.png")
#meanPassCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/meanPass.png")
widthFailCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/width.pdf")
meanFailCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/mean.pdf")
#widthPassCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/widthPass.pdf")
#meanPassCanvas.SaveAs("~/cernbox/www/latest_production/trigger_efficiency_TwoProngs_SimultaneousFit_latest_production/meanPass.pdf")

fOut = ROOT.TFile("histos/histos_trigger_efficiency/histos_TwoProngs/histos_TwoProngsSimFit.root","RECREATE")
fOut.cd()

h_efficiency_MC.Write()
h_efficiency_Data.Write()

fOut.Close()

fOut2 = ROOT.TFile("scale_factors/TwoProngsTriggerSF.root","RECREATE")
fOut2.cd()

h_efficiency_Data.Divide(h_efficiency_MC)
h_TwoProngsTrigger_SF = h_efficiency_Data
h_TwoProngsTrigger_SF.Write()

fOut2.Close()