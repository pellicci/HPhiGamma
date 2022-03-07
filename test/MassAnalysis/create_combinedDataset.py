import ROOT
import argparse
import os

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+") #import the Doube Crystal Ball PDF

debug = True

p = argparse.ArgumentParser(description='Select directory')
p.add_argument('dir_path', help='Type directory path')
args = p.parse_args()

#INPUT SETTINGS-------------------------------------------------------------
fInput1 = ROOT.TFile("workspaces/ws_signal.root")
fInput2 = ROOT.TFile("workspaces/ws_sidebands.root")
fInput3 = ROOT.TFile(args.dir_path + "histos_SR_Signal.root")

#TREE RETRIEVING------------------------------------------------------------
Signal_tree = fInput3.Get("tree_output")
Signal_entries = Signal_tree.GetEntriesFast()

#GET WORKSPACES
fInput1.cd()
workspaceSignal = fInput1.Get("myworkspace")
fInput2.cd()
workspaceData = fInput2.Get("myworkspace")

#workspaceSignal.Print()
#workspaceData.Print()

#RETRIEVE PDFs
sigPDF = workspaceSignal.pdf("signalPDF")
backgroundPDF = workspaceData.pdf("bkgPDF")

#N. OF EVENTS
signal_amplifier = 1.
signalWeight = 0.227113240117*signal_amplifier
Ndata = 4149
Nsig = ROOT.RooRealVar("Signal_events","N of signal events",signalWeight,-0.1,signalWeight+0.5*signalWeight) 
Nbkg = ROOT.RooRealVar("Bkg events", "N of bkg events",Ndata,Ndata-0.5*Ndata,Ndata+0.5*Ndata)

#ARGLISTs and VARs
argListPDFs = ROOT.RooArgList(sigPDF,backgroundPDF)
argListN = ROOT.RooArgList(Nsig,Nbkg)
mass_KKg = ROOT.RooRealVar("mass_KKg","The invariant mass_KKg",125.,100.,150.,"GeV/c^2")

#ADD PDFs
totPDF = ROOT.RooAddPdf("totPDF","The total PDF",argListPDFs,argListN)
print "PDFs added!"

#GENERATE DATA
dataCombined = totPDF.generate(ROOT.RooArgSet(mass_KKg),Ndata)

#PLOT
mass_KKgplot = mass_KKg.frame(40)
print "mass_KKg plot created!"
dataCombined.plotOn(mass_KKgplot)
print "dataCombined plotted!"
totPDF.plotOn(mass_KKgplot)
print "totpdf plotted!"

canvas = ROOT.TCanvas()
canvas.cd()
mass_KKgplot.Draw()
canvas.SaveAs("/eos/user/g/gumoret/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/PDFtot.pdf")
canvas.SaveAs("/eos/user/g/gumoret/www/MyAnalysis/HPhiGamma/MassAnalysis/latest_production/PDFtot.png")

#SAVE DATACOMBINED
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(dataCombined)

fOut = ROOT.TFile("workspaces/ws_dataCombined.root","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()

del workspace
