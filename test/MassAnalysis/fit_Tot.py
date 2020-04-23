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

#EVENTS COUNTING------------------------------------------------
Signal_events = 0.

#SIGNAL events loop
for jentry in xrange(Signal_entries):
    ientry = Signal_tree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = Signal_tree.GetEntry(jentry )
    if nb <= 0:
        continue

    Signal_events += Signal_tree._nEvents #s   

#GET WORKSPACES
fInput1.cd()
workspaceSignal = fInput1.Get("myworkspace")
fInput2.cd()
workspaceData = fInput2.Get("myworkspace")

#RETRIEVE PDFs
sigPDF = workspaceSignal.pdf("signalPDF")
backgroundPDF = workspaceData.pdf("bkgPDF")

#N. OF EVENTS
Nsig = ROOT.RooRealVar("Signal_events","N of signal events",Signal_events) 
Nbkg = ROOT.RooRealVar("Bkg events", "N of bkg events",789)

#ARGLISTs and VARs
argListPDFs = ROOT.RooArgList(sigPDF,backgroundPDF)
argListN = ROOT.RooArgList(Nsig,Nbkg)
mass = ROOT.RooRealVar("mass","The invariant mass",125.,100.,150.,"GeV/c^2")

#ADD PDFs
totPDF = ROOT.RooAddPdf("totPDF","The total PDF",argListPDFs,argListN)
print "PDFs added!"

#GENERATE DATA
data = totPDF.generate(ROOT.RooArgSet(mass),1000)
totPDF.fitTo(data)
print "Data fitted!"

#PLOT
massplot = mass.frame()
print "mass plot created!"
totPDF.plotOn(massplot)
print "totpdf plotted!"
data.plotOn(massplot)
print "data plotted!"

canvas = ROOT.TCanvas()
canvas.cd()
massplot.Draw()
canvas.SaveAs("PDFtot.pdf")
