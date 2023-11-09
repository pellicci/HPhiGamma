import ROOT
import argparse
import math
import numpy as np
import sys
import copy
from array import array
import tdrstyle, CMS_lumi
from ROOT import gROOT, RooFit
import os, sys

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True) 

p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('eta', help='Type the eta range')
args = p.parse_args()
ETA = args.eta

# B i n   L o o p
# -------------------------------------------------------------
pT_list = [38.,40.,42.,44.,46.,50.,55.,80.]
binList = [0,1,2,3,4,5,6]
#binList = [4]

# H i s t o s
# -------------------------------------------------------------
h_efficiency_Data = ROOT.TH1F("h_efficiency_Data","eff data",len(pT_list) - 1, array('d', pT_list))
h_efficiency_MC   = ROOT.TH1F("h_efficiency_MC","eff MC",len(pT_list) - 1, array('d', pT_list))

#Import the Doube Crystal Ball PDF -----------------------------------
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+")


for SAMPLE in ["Data", "MC"]:

	for binIndex in binList:

		# I n p u t   f i l e s
		# -------------------------------------------------------------
		fileInput = ROOT.TFile("histos/IsoEfficiency/histos_IsoNeuEfficiency_"+SAMPLE+"_"+ETA+".root")
		fileInput.cd()
		tree_fail = fileInput.Get("tree_output_bin"+str(binIndex)+"_fail")
		tree_pass = fileInput.Get("tree_output_bin"+str(binIndex)+"_pass")


		# C r e a t e   m o d e l  
		# -------------------------------------------------------------

		# Create observables
		mass = ROOT.RooRealVar("mumuMass","mumu invariant mass",80.,110.)

		# Construct signal pdf (Double-sided Crystal Ball)
		mean_central  ,  mean_min ,  mean_max  = 91., 89., 94.
		#sigma_central , sigma_min , sigma_max  = 0.8,0.5,1.
		alphaL_central, alphaL_min, alphaL_max = 3., 0.1, 50.
		alphaR_central, alphaR_min, alphaR_max = 1., 0.1, 50.
		nL_central    ,     nL_min,     nL_max = 3., 0.1, 50.
		nR_central    ,     nR_min,     nR_max = 1., 0.1, 50.

		#in case of Voigtian
		if binIndex == 0 or 2 or 3:
			sigma_central , sigma_min , sigma_max  = 0.84,0.5,2.
			width_central , width_min , width_max  = 2.9,2.5,5.

		else:
			sigma_central , sigma_min , sigma_max  = 0.8,0.5,1.
			width_central , width_min , width_max  = 3.55,2.5,5.

		#in case of Breit Wigner
		Zmean, Zwidth = 89., 2.495
		BWmean  = ROOT.RooRealVar("mean","The mean of the double Crystal Ball", Zmean, 85., 93.)
		BWwidth = ROOT.RooRealVar("width","The width of the Voigtian", Zwidth, 1., 10.)

		#BWmean.setConstant(1)
		#BWwidth.setConstant(1)
		
		mean   = ROOT.RooRealVar("mean","The mean of the double Crystal Ball", mean_central, mean_min, mean_max)
		sigma  = ROOT.RooRealVar("sigma","The width of the double Crystal Ball", sigma_central, sigma_min, sigma_max)
		alphaL = ROOT.RooRealVar("alphaL","The alpha left of the double Crystal Ball", alphaL_central, alphaL_min, alphaL_max)
		alphaR = ROOT.RooRealVar("alphaR","The alpha right of the double Crystal Ball", alphaR_central, alphaR_min, alphaR_max)
		nL     = ROOT.RooRealVar("nL","The n left of the double Crystal Ball", nL_central, nL_min, nL_max)
		nR     = ROOT.RooRealVar("nR","The n right of the double Crystal Ball", nR_central, nR_min, nR_max)
		
		#in case of Voigtian
		width  = ROOT.RooRealVar("width","The width of the Voigtian", width_central , width_min , width_max)

		#signalPDF = ROOT.RooDoubleCBFast("signalPDF","CB pdf",mass,mean,sigma,alphaL,alphaR,nL,nR)
		#gaussian = ROOT.RooGaussian("gaussian","CB pdf",mass,mean,sigma)
		signalPDF = ROOT.RooVoigtian("signalPDF","CB pdf",mass,mean,width,sigma)
		#BWpdf    = ROOT.RooBreitWigner("BWpdf","BW pdf",mass,BWmean,BWwidth)

		#signalPDF = gaussian
		#signalPDF = ROOT.RooFFTConvPdf("signalPDF", "BW + gauss", mass, BWpdf, gaussian)

		# Construct background pdf
		if SAMPLE == "Data":
			a1_central, a1_min, a1_max = 0.,-5.5,3.5
			a2_central, a2_min, a2_max = 0.,-5.5,0.5
			a3_central, a3_min, a3_max = 0.3,-2.,2.
		else:
			a1_central, a1_min, a1_max = 0.,-1.5,1.5
			a2_central, a2_min, a2_max = 0.,-1.5,0.5
			a3_central, a3_min, a3_max = 0.3,-2.,2.

		a1_fail = ROOT.RooRealVar("a1_fail","The a1 of background", a1_central, a1_min, a1_max)
		a2_fail = ROOT.RooRealVar("a2_fail","The a2 of background", a2_central, a2_min, a2_max)
		a3_fail = ROOT.RooRealVar("a3_fail","The a3 of background", a3_central, a3_min, a3_max)

		a1_pass = ROOT.RooRealVar("a1_pass","The a1 of background", a1_central, a1_min, a1_max)
		a2_pass = ROOT.RooRealVar("a2_pass","The a2 of background", a2_central, a2_min, a2_max)
		a3_pass = ROOT.RooRealVar("a3_pass","The a3 of background", a3_central, a3_min, a3_max)

		backgroundPDFFail = ROOT.RooChebychev("backgroundPDFFail","The background PDF for fail dataset",mass,ROOT.RooArgList(a1_fail, a2_fail))
		backgroundPDFPass = ROOT.RooChebychev("backgroundPDFPass","The background PDF for pass dataset",mass,ROOT.RooArgList(a1_pass,a2_pass))

		# N   e v e n t s
		# ---------------------------------------------------------------
		nEntries = tree_pass.GetEntriesFast()

		nsig = ROOT.RooRealVar("nsig", "signal yield 1", nEntries, nEntries/3, nEntries*2)
		if binIndex == 4: nsig = ROOT.RooRealVar("nsig", "signal yield 1", nEntries, nEntries/2, nEntries*1.5)
		#nsig = ROOT.RooRealVar("nsig", "signal yield 1", 1000000, 1000000/3, 1000000*2)
		#nsig2 = ROOT.RooRealVar("nsig2", "signal yield 2", 300., 0., 5000.)
		if SAMPLE == "MC":
			nbkgFail   = ROOT.RooRealVar("nbkgFail", "background yield", nEntries/100, 0., nEntries)
			nbkgPass   = ROOT.RooRealVar("nbkgPass", "background yield", nEntries/100, 0., nEntries)
		else:
			nbkgFail   = ROOT.RooRealVar("nbkgFail", "background yield", nEntries/8, nEntries/1000, nEntries/2)
			nbkgPass   = ROOT.RooRealVar("nbkgPass", "background yield", nEntries/8, nEntries/1000, nEntries/2)
		
		efficiency = ROOT.RooRealVar("efficiency", "efficiency", 0.8, 0.4, 1.)

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
		sample.defineType("Pass")
		sample.defineType("Fail")

		jointData = ROOT.RooDataSet("jointData", "joint data", ROOT.RooArgSet(mass), ROOT.RooFit.Index(sample), ROOT.RooFit.Import("Pass", data2), ROOT.RooFit.Import("Fail", data1))

		# C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
		# -----------------------------------------------------------------------------------
		simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)

		# Create model for each dataset
		#modelPass = ROOT.RooExtendPdf("modelPass", "Signal PDF for pass", signalPDF, signal_yield2_formula)
		#modelFail = ROOT.RooExtendPdf("modelFail", "Signal PDF for fail", signalPDF, signal_yield1_formula)
		
		modelPass = ROOT.RooAddPdf("modelPass", "Signal and Background PDF for pass", ROOT.RooArgList(signalPDF, backgroundPDFPass), ROOT.RooArgList(signal_yield2_formula, nbkgPass))
		modelFail = ROOT.RooAddPdf("modelFail", "Signal and Background PDF for fail", ROOT.RooArgList(signalPDF, backgroundPDFFail), ROOT.RooArgList(signal_yield1_formula, nbkgFail))

		simPdf.addPdf(modelPass, "Pass")
		simPdf.addPdf(modelFail, "Fail")

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
		frame1 = mass.frame(120)
		frame1.SetTitle("Fail selection")
		jointData.plotOn(frame1, ROOT.RooFit.Cut("sample==sample::Fail"), ROOT.RooFit.Name("data1"))
		simPdf.plotOn(frame1, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Fail"), ROOT.RooFit.Components("backgroundPDFFail"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("curve1_bkg"))
		simPdf.plotOn(frame1, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Fail"), ROOT.RooFit.Name("curve1_total"))
		frame1.GetXaxis().SetLabelSize(0.033)
		frame1.GetXaxis().SetTitle("m_{#mu#mu} [Gev]")
		frame1.SetMaximum(1.8 * frame1.GetMaximum())

		#FRAME2
		frame2 = mass.frame(120)
		frame2.SetTitle("Fail selection and trigger")
		jointData.plotOn(frame2, ROOT.RooFit.Cut("sample==sample::Pass"), ROOT.RooFit.Name("data2"))
		simPdf.plotOn(frame2, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Pass"), ROOT.RooFit.Components("backgroundPDFPass"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("curve2_bkg"))
		simPdf.plotOn(frame2, ROOT.RooFit.ProjWData(ROOT.RooArgSet(sample), jointData), ROOT.RooFit.Slice(sample, "Pass"), ROOT.RooFit.Name("curve2_total"))
		frame2.GetXaxis().SetLabelSize(0.033)
		frame2.GetXaxis().SetTitle("m_{#mu#mu} [Gev]")
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
		legend1.AddEntry(0, "Fail, "+ETA, "")
		legend1.AddEntry(frame1.findObject("curve1_total"), "Fit", "l")
		legend1.AddEntry(frame1.findObject("curve1_bkg"), "Bkg-only", "l")
		legend1.AddEntry(0, "#chi^{{2}}/ndf = {:.2f}".format(chi2_ndf1), "")
		legend1.AddEntry(0,"N_{sig} = "+str(round(nsig1,1))+" #pm "+str(round(nsig1_err,1)),"brNDC")
		legend1.AddEntry(0,"N_{bkg} = "+str(round(nbkgFail.getValV(),1))+" #pm "+str(round(nbkgFail.errorVar().getValV(),1)),"brNDC")
		legend1.AddEntry(0,"1 - #epsilon = "+ str(round(1-eff,3))+" #pm "+str(round(eff_err,3)),"brNDC")
		legend1.AddEntry(0,"width = "+str(round(sigma_val,2))+" #pm "+str(round(sigma_err,2)),"brNDC")

		legend2 = ROOT.TLegend(x0, y0, x1, y1)
		legend2.SetBorderSize(0)
		legend2.SetFillStyle(0)
		legend2.SetTextSize(0.04)
		legend2.AddEntry(0, "Pass, "+ETA, "")
		legend2.AddEntry(frame2.findObject("curve2_total"), "Fit", "l")
		legend2.AddEntry(frame2.findObject("curve2_bkg"), "Bkg-only", "l")
		legend2.AddEntry(0, "#chi^{{2}}/ndf = {:.2f}".format(chi2_ndf2), "")
		legend2.AddEntry(0,"N_{sig} = "+str(round(nsig2,1))+" #pm "+str(round(nsig2_err,1)),"brNDC")
		legend2.AddEntry(0,"N_{bkg} = "+str(round(nbkgPass.getValV(),1))+" #pm "+str(round(nbkgPass.errorVar().getValV(),1)),"brNDC")
		legend2.AddEntry(0,"#epsilon = "+str(round(eff,3))+" #pm "+str(round(eff_err,3)),"brNDC")
		legend2.AddEntry(0,"width = "+str(round(sigma_val,2))+" #pm "+str(round(sigma_err,2)),"brNDC")

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

		if SAMPLE == "Data":
			#h_efficiency_Data.Fill(binIndex,eff)
			#h_efficiency_Data.SetBinError(binIndex,eff_err)

			h_efficiency_Data.SetBinContent(binIndex+1, eff)
			h_efficiency_Data.SetBinError(binIndex+1, eff_err)
		
		else:
			#h_efficiency_MC.Fill(binIndex,eff)
			#h_efficiency_MC.SetBinError(binIndex,eff_err)
			h_efficiency_MC.SetBinContent(binIndex+1, eff)
			h_efficiency_MC.SetBinError(binIndex+1, eff_err)
	
		canvas.SaveAs("~/cernbox/www/latest_production/IsoNeuEfficiency_Fit_"+ETA+"_latest_production/simultaneous_fit_bin_"+str(binIndex)+"_"+SAMPLE+".png")
		canvas.SaveAs("~/cernbox/www/latest_production/IsoNeuEfficiency_Fit_"+ETA+"_latest_production/simultaneous_fit_bin_"+str(binIndex)+"_"+SAMPLE+".pdf")


fOut = ROOT.TFile("histos/IsoEfficiency/histos_IsoNeuEfficiencySimFit_"+ETA+".root","RECREATE")
fOut.cd()

h_efficiency_MC.Write()
h_efficiency_Data.Write()

fOut.Close()

