import ROOT
import argparse
import math
import numpy as np
import sys
import copy
import tdrstyle, CMS_lumi


#bools
debug = True #Bool for verbose

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
mytree = fInput.Get("HPhiGammaTriggerAnalysis/mytree")

#Output
output_filename = args.outputfile_option
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_mass_mumu","h_gamma_eT","h_nPhotonTriggered_eT","h_nPhotonOffline_eT","h_triggerEff_eT"]

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{#mu#mu}",100,20.,120.)
histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"E_{T}^{#gamma}", 100, 30.,160.)
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"n. events #gamma triggered as function of E_{#gamma}", 11, 35.,90.)
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"n. events #gamma over offline selection as function of E_{#gamma}", 11, 35.,90.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"Trigger efficiency as function of E_{#gamma}", 11, 35.,90.)

#VARIABLES INITIALIZATION ##################################################################################################
nEventsIsoMuTrigger  = 0
nEventsPhotonTrigger = 0
nEventsOfflinePhoton = 0
triggerFraction      = 0

#EVENTS LOOP ##################################################################################################################### 
print "This sample has ", mytree.GetEntriesFast(), " events"
nentries = mytree.GetEntriesFast()

for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry )
    if nb <= 0:
        continue

    #Retrieve variables from the tree and fill the histos
    MuMuMass              = mytree.MuMuMass
    gammaEt               = mytree.photon_eT
    isIsoMuTrigger        = mytree.isIsoMuTrigger
    isPhoton30Mu17Trigger = mytree.isPhotonTrigger      
    isbestPhotonFound     = mytree.cand_photon_found
    
    if isbestPhotonFound == 0: continue

    if debug:
        print ""
        print "Processing EVENT n.",jentry+1,"================"
        print "isIsoMuTrigger = ",isIsoMuTrigger
        print "isPhoton30Mu17Trigger = ",isPhoton30Mu17Trigger
        print "isbestPhotonFound = ",isbestPhotonFound
    
    if isIsoMuTrigger:
        nEventsIsoMuTrigger  = nEventsIsoMuTrigger  + 1
    if isPhoton30Mu17Trigger:
        nEventsPhotonTrigger = nEventsPhotonTrigger + 1
    if isbestPhotonFound:
        nEventsOfflinePhoton = nEventsOfflinePhoton + 1


    histo_map["h_mass_mumu"].Fill(MuMuMass)        
    histo_map["h_gamma_eT"].Fill(gammaEt)
    if isPhoton30Mu17Trigger: histo_map["h_nPhotonTriggered_eT"].Fill(gammaEt)
    if isbestPhotonFound: histo_map["h_nPhotonOffline_eT"].Fill(gammaEt)

#TRIGGER EFFICIENCY CALCULATION ###################################################################################################
triggerFraction = float(nEventsPhotonTrigger)/float(nEventsOfflinePhoton)

#trigger efficiency as function of photon transverse energy -------------------------------------------------------------
histo_map["h_triggerEff_eT"] = histo_map["h_nPhotonTriggered_eT"] #copy the histo of the numerator
histo_map["h_triggerEff_eT"].Divide(histo_map["h_nPhotonOffline_eT"]) #divide it by the denominator
histo_map["h_triggerEff_eT"].Scale(100.) #multiply for the percentage

if debug:
    print "########################################"
    print "nEventsPhotonTrigger = ",nEventsPhotonTrigger
    print "nEventsOfflinePhoton = ",nEventsOfflinePhoton
    print "Trigger fraction = ",triggerFraction
    print "########################################"

#HISTO LABELS #####################################################################################################################
histo_map["h_mass_mumu"].GetXaxis().SetTitle("m_{#mu^{+}#mu^{-}} [GeV]")
histo_map["h_gamma_eT"].GetXaxis().SetTitle("E^{#gamma}_{T} [GeV]")

histo_map["h_nPhotonTriggered_eT"].GetXaxis().SetTitle("E^{#gamma}_{T} [GeV]")
histo_map["h_nPhotonTriggered_eT"].GetYaxis().SetTitle("n. #gamma triggered")

histo_map["h_nPhotonOffline_eT"].GetXaxis().SetTitle("E^{#gamma}_{T} [GeV]")
histo_map["h_nPhotonOffline_eT"].GetYaxis().SetTitle("n. offline #gamma")

histo_map["h_triggerEff_eT"].GetXaxis().SetTitle("E^{#gamma}_{T} [GeV]")
histo_map["h_triggerEff_eT"].GetYaxis().SetTitle("Trigger efficiency")
histo_map["h_triggerEff_eT"].SetTitle("Trigger efficiency as function of E_{#gamma}")

# HISTOS WRITING #########################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()
fOut.Close()


#c1 = ROOT.TCanvas()
#c1.cd()
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
#ROOT.gStyle.SetOptStat(0)
#histo_map["h_mass_mumu"].SetMarkerSize(1.4)
#histo_map["h_mass_mumu"].GetXaxis().SetRangeUser(20.,120.)
#histo_map["h_mass_mumu"].GetYaxis().SetRangeUser(0.,30.)
#histo_map["h_mass_mumu"].Draw("")
#c1.SaveAs("~/cernbox/www/testDir/h_mass_mumu.pdf")
#c1.SaveAs("~/cernbox/www/testDir/h_mass_mumu.png")  

#c2 = ROOT.TCanvas()
#c2.cd()
##ROOT.gStyle.SetPaintTextFormat("4.2f %")
#ROOT.gStyle.SetOptStat(0)
#histo_map["h_gamma_eT"].SetMarkerSize(1.4)
#histo_map["h_gamma_eT"].GetXaxis().SetRangeUser(20.,120.)
#histo_map["h_gamma_eT"].GetYaxis().SetRangeUser(0.,30.)
#histo_map["h_gamma_eT"].Draw("")
#c2.SaveAs("~/cernbox/www/testDir/h_gamma_eT.png") 
#c2.SaveAs("~/cernbox/www/testDir/h_gamma_eT.pdf") 

#c3 = ROOT.TCanvas()
#c3.cd()
#ROOT.gStyle.SetOptStat(0)
#histo_map["h_nPhotonTriggered_eT"].Draw("")
#c3.SaveAs("~/cernbox/www/testDir/h_nPhotonTriggered_eT.pdf")
#c3.SaveAs("~/cernbox/www/testDir/h_nPhotonTriggered_eT.png") 

#c4 = ROOT.TCanvas()
#c4.cd()
#ROOT.gStyle.SetOptStat(0)
#histo_map["h_nPhotonOffline_eT"].Draw("")
#c4.SaveAs("~/cernbox/www/testDir/h_nPhotonOffline_eT.pdf")
#c4.SaveAs("~/cernbox/www/testDir/h_nPhotonOffline_eT.png") 

#c3 = ROOT.TCanvas()
#c3.cd()
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPaintTextFormat("4.2f %")
#histo_map["h_triggerEff_eT"].GetYaxis().SetRangeUser(50.,100.)
#histo_map["h_triggerEff_eT"].Draw("HIST TEXT0")
#c3.SaveAs("~/cernbox/www/testDir/h_triggerEff_eT.pdf")
#c3.SaveAs("~/cernbox/www/testDir/h_triggerEff_eT.png") 