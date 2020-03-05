import ROOT
import sys
import copy

debug = False
#debug = True

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

signal_magnify = int(sys.argv[1])

list_inputfiles = []
for filename in sys.argv[2:]:
    list_inputfiles.append(filename)

output_dir = "plots/latest_production/"

hstack  = dict()
hsignal = dict()
hdata   = dict()
canvas  = dict()

histo_container = [] #just for memory management

#rename histos
list_histos = ["h_InvMass_TwoTrk_Photon","h_phi_InvMass_TwoTrk","h_InvMass_TwoTrk_Photon_NoPhiMassCut","h_firstKCand_pT","h_secondKCand_pT","h_secondKCand_Eta","h_firstKCand_Phi","h_secondKCand_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestJetPt","h_bestJetEta","h_K1_Iso","h_K1_Iso_ch","h_K2_Iso","h_K2_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons"] 

for hname in list_histos:
    hstack[hname] = ROOT.THStack("hstack_" + hname,"")

# Color mask must have the same number of entries as non-QCD backgrounds + 1 (that is the cumulative QCD background)
colors_mask = dict()
colors_mask["ttbarToSemiLeptonic"] = ROOT.kAzure+5
colors_mask["ttbarlnu"]            = ROOT.kAzure+7
colors_mask["ttbarToHadronic"]     = ROOT.kAzure+9
colors_mask["DY10to50"]            = ROOT.kGreen-4
colors_mask["DY50"]                = ROOT.kGreen-6
colors_mask["GammaJetsEnriched"]   = ROOT.kOrange+1
colors_mask["GammaJets20to40"]     = ROOT.kOrange+1
colors_mask["GammaJets20toInf"]    = ROOT.kOrange+1
colors_mask["GammaJets40toInf"]    = ROOT.kOrange+1
colors_mask["GammaJets"]           = ROOT.kYellow+1
colors_mask["GammaJetsHT40to100"]  = ROOT.kYellow+1
colors_mask["GammaJetsHT100to200"] = ROOT.kYellow+1
colors_mask["GammaJetsHT200to400"] = ROOT.kYellow+1
colors_mask["GammaJetsHT400to600"] = ROOT.kYellow+1
colors_mask["GammaJetsHT600toInf"] = ROOT.kYellow+1
colors_mask["WZ"]                  = ROOT.kPink+1
colors_mask["WW"]                  = ROOT.kPink+3
colors_mask["WJetsToLNu0J"]        = ROOT.kCyan-7
colors_mask["WJetsToLNu1J"]        = ROOT.kCyan-8
colors_mask["WJetsToLNu2J"]        = ROOT.kCyan-9
colors_mask["QCD"]                 = ROOT.kBlue+1
colors_mask["QCDHT100to200"]       = ROOT.kBlue+1
colors_mask["QCDHT200to300"]       = ROOT.kBlue+1
colors_mask["QCDHT300to500"]       = ROOT.kBlue+1
colors_mask["QCDHT500to700"]       = ROOT.kBlue+1
colors_mask["QCDHT700to1000"]      = ROOT.kBlue+1
colors_mask["QCDHT1000to1500"]     = ROOT.kBlue+1
colors_mask["QCDHT1500to2000"]     = ROOT.kBlue+1
colors_mask["QCDHT2000toInf"]      = ROOT.kBlue+1
colors_mask["DiPhotonJets"]        = ROOT.kRed+1


# leg1 = ROOT.TLegend(0.15,0.6120093,0.34,0.9491917) #left positioning
leg1 = ROOT.TLegend(0.6868687,0.6520093,0.9511784,0.9091917) #right positioning
leg1.SetHeader("") 
leg1.SetFillColor(0)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)

#LOOP ON INPUT FILES (MC backgrounds, MC Signal, DATA)
for filename in list_inputfiles:
    fileIn = ROOT.TFile(filename)

    sample_name = (filename.split("_")[2])[:-5] #retireve the samplename (QCD, DY, Data, etc...)
    if(debug):
        print "========="+sample_name+"========"

    #LOOP ON HISTOS (for each sample)
    for histo_name in list_histos:
        histo = fileIn.Get(histo_name)
       
        #REBIN
        if histo_name == "h_phi_InvMass_TwoTrk":
            histo.Rebin(5) 
        elif histo_name == "h_couple_Iso":
            histo.Rebin(5) 
        elif histo_name == "h_couple_Iso_ch":
            histo.Rebin(5) 
        elif histo_name == "h_nJets_25" or histo_name == "h_nMuons" or histo_name == "h_nElectrons" or histo_name == "h_nPhotons" :
            histo.Rebin(1) 
        else:
            histo.Rebin(2) 
        
        if(debug):        
            print histo_name

        # Set to 0 the bins containing negative values, due to negative weights
        hsize = histo.GetSize() - 2 # GetSize() returns the number of bins +2 (that is + overflow + underflow) 
        for bin in range(1,hsize+1): # The +1 is in order to get the last bin
            bincontent = histo.GetBinContent(bin)
            if bincontent < 0.:
                histo.SetBinContent(bin,0.)
        
        histo_container.append(copy.copy(histo))
        
        if "Signal" in sample_name:
            histo_container[-1].SetLineStyle(2)   #dashed
            histo_container[-1].SetLineColor(2)   #red
            histo_container[-1].SetLineWidth(4)   #kind of thick
            hsignal[histo_name] = histo_container[-1]
        if "Data" in sample_name:
            histo_container[-1].SetMarkerStyle(20)   #point
            hdata[histo_name] = histo_container[-1]
        if not sample_name == "Data" and not sample_name == "Signal":
            histo_container[-1].SetFillColor(colors_mask[sample_name])
            hstack[histo_name].Add(histo_container[-1])
            
    #fill legend
    #if histo_name == "h_nMuons": #Add the legend only once (h_nMuons is just a random variable)

    #if histo.Integral() > float(signal_magnify)/12. or sample_name == "Signal": #Only plot in the legend those samples which have some contribution
    if not sample_name == "Data" and not sample_name == "Signal":
        leg1.AddEntry(histo_container[-1],sample_name,"f")
    elif sample_name == "Data":
        leg1.AddEntry(histo_container[-1],sample_name,"lep")
    elif sample_name == "Signal":
        leg1.AddEntry(histo_container[-1],sample_name + " x " + str(signal_magnify),"f")

    fileIn.Close()


for histo_name in list_histos:
    if(debug):
        print "in second loop: "+histo_name
    canvas[histo_name] = ROOT.TCanvas("Canvas_" + histo_name,"",200,106,600,600)
    canvas[histo_name].cd()

    ##########################################
    pad1 = ROOT.TPad("pad_" + histo_name,"",0,0.28,1,1)
    pad2 = ROOT.TPad("pad_" + histo_name,'',0,0,1,0.25)
    pad1.SetBottomMargin(0.02)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    pad1.Draw()
    pad2.Draw()
    if histo_name == "h_nJets_25" or histo_name == "h_nMuons" or histo_name == "h_nElectrons": # or histo_name == "h_nPhotons"
        pad1.SetLogy()
        
    ##########################################

    ##########################################
    pad1.cd()
    ##########################################

    hstack[histo_name].SetTitle(histo_name)
    hstack[histo_name].Draw("histo")

    ##########################################
    #hstack[histo_name].GetXaxis().SetTickLength(0)
    hstack[histo_name].GetXaxis().SetLabelOffset(999)
    ##########################################

    if histo_name == "h_InvMass_TwoTrk_Photon" :
        hstack[histo_name].GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^2]")
        hstack[histo_name].GetXaxis().SetRangeUser(100.,150.)
        hstack[histo_name].SetTitle("Tracks + #gamma invariant mass (Cut on #phi inv. mass)")
        pad1.SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^2]")

    if histo_name == "h_nJets_25" :
        hstack[histo_name].GetXaxis().SetTitle("nJets_25")
        hstack[histo_name].SetTitle("n Jets over 25 GeV")

    if histo_name == "h_InvMass_TwoTrk_Photon_NoPhiMassCut" :
        hstack[histo_name].GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^{2}]")
        hstack[histo_name].GetXaxis().SetRangeUser(100.,150.)
        hstack[histo_name].SetTitle("Tracks + #gamma invariant mass (No cuts on #phi inv mass)")
    
    if histo_name == "h_phi_InvMass_TwoTrk" :
        hstack[histo_name].GetXaxis().SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]")
        hstack[histo_name].GetXaxis().SetRangeUser(1.,1.05)
        hstack[histo_name].SetTitle("Two tracks invariant mass")

    if histo_name == "h_firstKCand_pT" :
        hstack[histo_name].GetXaxis().SetTitle("pT_{K} [GeV/c]")
        hstack[histo_name].GetXaxis().SetRangeUser(0.,150.)
        hstack[histo_name].SetTitle("Transverse momentum of the first charged particle (pT_{max} of the couple)")

    if histo_name == "h_secondKCand_pT" :
        hstack[histo_name].GetXaxis().SetTitle("pT_{K} [GeV/c]")
        hstack[histo_name].GetXaxis().SetRangeUser(0.,150.)
        hstack[histo_name].SetTitle("Transverse momentum of the second charged particle (pT_{min} of the couple)")
    
    if histo_name == "h_firstKCand_Eta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta")
        hstack[histo_name].GetXaxis().SetRangeUser(-2.5,2.5)
        hstack[histo_name].SetTitle("Pseudorapidity of the first charged particle (pT_{max} of the couple)")

    if histo_name == "h_secondKCand_Eta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta")
        hstack[histo_name].GetXaxis().SetRangeUser(-2.5,2.5)
        hstack[histo_name].SetTitle("Pseudorapidity of the second charged particle (pT_{min} of the couple)")

    if histo_name == "h_firstKCand_Phi" :
        hstack[histo_name].GetXaxis().SetTitle("#phi [rad]")
        hstack[histo_name].SetTitle("Azimuthal angle of the first charged particle (pT_{max} of the couple)")

    if histo_name == "h_secondKCand_Phi" :
        hstack[histo_name].GetXaxis().SetTitle("#phi [rad]")
        hstack[histo_name].SetTitle("Azimuthal angle of the second charged particle (pT_{min} of the couple)")

    if histo_name == "h_bestCouplePt" :
        hstack[histo_name].GetXaxis().SetTitle("pT_{K^{+}K^{-}} [GeV]")
        hstack[histo_name].SetTitle("Transverse momentum of the couple")

    if histo_name == "h_bestCoupleEta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta")
        hstack[histo_name].SetTitle("Pseudorapidity of the couple")

    if histo_name == "h_bestJetPt" :
        hstack[histo_name].GetXaxis().SetTitle("pT_{jet} [GeV]")
        hstack[histo_name].SetTitle("Transverse momentum of the jet")

    if histo_name == "h_bestJetEta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{jet}")
        hstack[histo_name].SetTitle("Pseudorapidity of the jet")

    if histo_name == "h_K1_Iso" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{1}}/pT_{K_{1}}")
        hstack[histo_name].SetTitle("Isolation of the K_{1} candidate")

    if histo_name == "h_K2_Iso" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")
        hstack[histo_name].SetTitle("Isolation of the K_{2} candidate")

    if histo_name == "h_K1_Iso_ch" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{1}}/pT_{K_{1}}")
        hstack[histo_name].SetTitle("Isolation of the K_{1} candidate")

    if histo_name == "h_K2_Iso_ch" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")
        hstack[histo_name].SetTitle("Isolation of the K_{2} candidate")

    if histo_name == "h_couple_Iso" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")
        hstack[histo_name].SetTitle("Isolation of the couple candidate")

    if histo_name == "h_couple_Iso_ch" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")
        hstack[histo_name].SetTitle("Isolation of the couple candidate")

    if histo_name == "h_bestCoupleDeltaR" :
        hstack[histo_name].GetXaxis().SetTitle("#DeltaR_{K^{+}K^{-}}")
        hstack[histo_name].SetTitle("#DeltaR of the couple")

    if histo_name == "h_nPhotons" :
        hstack[histo_name].GetXaxis().SetTitle("n.#gamma")
        hstack[histo_name].SetTitle("n.#gamma over selections")


    if signal_magnify != 1:
        hsignal[histo_name].Scale(signal_magnify)      

    hstack[histo_name].Draw("histo")
    hsignal[histo_name].Draw("SAME,hist")
    hdata[histo_name].Draw("SAME,E1")
    
    hMCErr = copy.deepcopy(hstack[histo_name].GetStack().Last())
    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(1)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.Draw("sameE2")

    leg1.Draw()

    ################################################
    pad2.cd()
    pad2.SetTopMargin(0)
    pad2.SetFillColor(0)
    pad2.SetFillStyle(0)
    ROOT.gStyle.SetOptStat(0)
    totalMC = copy.deepcopy(hstack[histo_name].GetStack().Last())
    totalData = copy.deepcopy(hdata[histo_name])
    totalData.Divide(totalMC)
    totalData.SetTitle("")
    totalData.SetMarkerStyle(8)
    totalData.SetMarkerColor(1)
    totalData.SetLineColor(1)
    totalData.GetYaxis().SetLabelSize(0.10)
    totalData.GetYaxis().SetTitle("Data/MC")
    totalData.GetYaxis().SetTitleSize(0.08)
    totalData.GetYaxis().SetTitleOffset(0.5)
    #totalData.GetYaxis().SetRangeUser(0.4,1.6)

    totalData.GetXaxis().SetLabelSize(0.10)
    totalData.GetXaxis().SetTitleSize(0.12)
    totalData.GetXaxis().SetTitleOffset(1.0)

    if histo_name == "h_InvMass_TwoTrk_Photon":
        totalData.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV/c^{2}]")

    line_on_one = ROOT.TLine(totalData.GetXaxis().GetXmin(),1.,totalData.GetXaxis().GetXmax(),1.)
    line_on_one.SetLineColor(38)

    totalData.Draw()
    line_on_one.Draw("SAME")
        
    ################################################

    canvas[histo_name].SaveAs(output_dir + histo_name + ".pdf")
