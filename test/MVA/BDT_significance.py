import ROOT
import math
import numpy as np
import argparse

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)  

Nsig_passed = 0.52596985 # Number of signal and background events from the sum of the weights (before applying BDT cuts), take this number running the generate_histos for signal
Nbkg_passed = 11926

BDT_file = ROOT.TFile("outputs/Nominal_training.root")

h_BDT_effB_effS = BDT_file.Get("default/Method_BDT/BDT/MVA_BDT_effBvsS")

canvas1 = ROOT.TCanvas()
sig_eff = []
bkg_eff = []
signif = []
_effS = 0

for jbin in range(1,h_BDT_effB_effS.GetNbinsX()+1):
    if h_BDT_effB_effS.GetBinCenter(jbin) > 0.3:
        sig_eff.append(h_BDT_effB_effS.GetBinCenter(jbin))
        if h_BDT_effB_effS.GetBinContent(jbin) < 0.:
            bkg_eff.append(0.)
            signif.append(0)
        else:
            bkg_eff.append(h_BDT_effB_effS.GetBinContent(jbin))
            signif.append(Nsig_passed*h_BDT_effB_effS.GetBinCenter(jbin)/math.sqrt(Nbkg_passed*h_BDT_effB_effS.GetBinContent(jbin)))


sig_eff_array = np.array(sig_eff)
bkg_eff_array = np.array(bkg_eff)
signif_array = np.array(signif)
#print "signif_len: ", len(signif_array)
sign = ROOT.TGraph(70,sig_eff_array,signif_array)
sign.SetTitle("")
sign.GetXaxis().SetTitle("#varepsilon_{S}^{BDT}")
sign.GetYaxis().SetTitle("Significance")
sign.SetMaximum(0.04)
sign.SetMarkerStyle(8)
sign.SetMarkerColor(4)
sign.Draw("AP")

canvas2 = ROOT.TCanvas()
sign_vs_bkg = ROOT.TGraph(70,bkg_eff_array,signif_array)
sign_vs_bkg.SetTitle("")
sign_vs_bkg.GetXaxis().SetTitle("#varepsilon_{B}^{BDT}")
sign_vs_bkg.GetYaxis().SetTitle("Significance")
sign_vs_bkg.GetXaxis().SetRangeUser(0.,0.02)
sign_vs_bkg.SetMaximum(0.04)
sign_vs_bkg.SetMarkerStyle(8)
sign_vs_bkg.SetMarkerColor(4)
sign_vs_bkg.Draw("AP")


canvas1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MVA/latest_production/signif_vs_effS.pdf")
canvas1.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MVA/latest_production/signif_vs_effS.png")
canvas2.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MVA/latest_production/signif_vs_effB.pdf")
canvas2.SaveAs("~/cernbox/www/MyAnalysis/HPhiGamma/MVA/latest_production/signif_vs_effB.png")

#----Now find the BDT output corresponding to the highest significance

h_BDT_effS = BDT_file.Get("default/Method_BDT/BDT/MVA_BDT_effS")
signif_maximizing_eff = sig_eff_array[np.argmax(signif_array)]
BDT_output = 0. 

for entry in xrange(h_BDT_effS.GetNbinsX()):

    effS = h_BDT_effS.GetBinContent(entry)
    effS = float(format(effS, '.3f'))
    signif_maximizing_eff = float(format(signif_maximizing_eff, '.3f'))
    #print "effS: ", effS#, "signif_max_eff: ", signif_maximizing_eff
    if effS == signif_maximizing_eff:
    #if effS == 0.55:
        BDT_output =  h_BDT_effS.GetBinCenter(entry)
        _effS = effS

print "For a signal efficiency of ", _effS, "the BDT output is :", BDT_output

raw_input()
