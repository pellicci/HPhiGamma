import ROOT
import argparse
import os

debug = True

p = argparse.ArgumentParser(description='Select directory')
p.add_argument('dir_path', help='Type directory path')
args = p.parse_args()

#INPUT SETTINGS-------------------------------------------------------------
fInput1 = ROOT.TFile("workspaces/ws_signal.root")
fInput2 = ROOT.TFile("workspaces/ws_data.root")
fInput3 = ROOT.TFile(args.dir_path + "histos_Signal.root")

#TREE RETRIEVING------------------------------------------------------------
Signal_tree = fInput3.Get("tree_output")
Signal_entries = Signal_tree.GetEntriesFast()

#GET WORKSPACES
fInput1.cd()
workspaceSignal = fInput1.Get("myworkspace")
fInput2.cd()
workspaceData = fInput2.Get("myworkspace")

workspaceSignal.Print()
workspaceData.Print()

#RETRIEVE PDFs
sigPDF = workspaceSignal.pdf("signalPDF")
backgroundPDF = workspaceData.pdf("bkgPDF")

#N. OF EVENTS
signalWeight = 0.692480078053
Ndata = 789
Nsig = ROOT.RooRealVar("Signal_events","N of signal events",signalWeight,signalWeight-0.5*signalWeight,signalWeight+0.5*signalWeight) 
Nbkg = ROOT.RooRealVar("Bkg events", "N of bkg events",Ndata,Ndata-0.5*Ndata,Ndata+0.5*Ndata)

#ARGLISTs and VARs
argListPDFs = ROOT.RooArgList(sigPDF,backgroundPDF)
argListN = ROOT.RooArgList(Nsig,Nbkg)
mass_KKg = ROOT.RooRealVar("mass_KKg","The invariant mass_KKg",125.,100.,150.,"GeV/c^2")

#ADD PDFs
totPDF = ROOT.RooAddPdf("totPDF","The total PDF",argListPDFs,argListN)
print "PDFs added!"

#GENERATE DATA
dataset = totPDF.generate(ROOT.RooArgSet(mass_KKg),1000)
#totPDF.fitTo(dataset)

#PLOT
mass_KKgplot = mass_KKg.frame()
print "mass_KKg plot created!"
totPDF.plotOn(mass_KKgplot)
print "totpdf plotted!"
dataset.plotOn(mass_KKgplot)
print "dataset plotted!"

canvas = ROOT.TCanvas()
canvas.cd()
mass_KKgplot.Draw()
canvas.SaveAs("PDFtot.pdf")
