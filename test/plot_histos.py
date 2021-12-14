import ROOT
import sys
import copy

debug = False
debug = True
pad2Flag = False

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

signal_magnify = int(sys.argv[1])

list_inputfiles = []
for filename in sys.argv[2:]:
    list_inputfiles.append(filename)

output_dir = "~/cernbox/Analysis_plots/latest_production/"
output_dir = "~/cernbox/www/latest_production/"

hstack  = dict()
hsignal = dict()
hdata   = dict()
canvas  = dict()

histo_container = [] #just for memory management

#rename histos
list_histos = ["h_InvMass_TwoTrk_Photon","h_InvMass_TwoTrk_Photon_NoPhiMassCut","h_phi_InvMass_TwoTrk","h_firstKCand_pT","h_secondKCand_pT","h_firstKCand_Eta","h_secondKCand_Eta","h_firstKCand_Phi","h_secondKCand_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestCoupleDeltaR","h_bestJetPt","h_bestJetEta","h_K1_Iso","h_K1_Iso_ch","h_K2_Iso","h_K2_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons","h_photonWP90"]

for hname in list_histos:
    hstack[hname] = ROOT.THStack("hstack_" + hname,"")

# Color mask must have the same number of entries as non-QCD backgrounds + 1 (that is the cumulative QCD background)
colors_mask = dict()
colors_mask["ttbarToSemiLeptonic"] = ROOT.kAzure+5
colors_mask["ttbarlnu"]            = ROOT.kAzure+7
colors_mask["ttbarToHadronic"]     = ROOT.kAzure+9
colors_mask["DY10to50"]            = ROOT.kGreen-4
colors_mask["DY50"]                = ROOT.kGreen-6
colors_mask["GammaJetsHT100to200"] = ROOT.kOrange+1
colors_mask["GammaJetsHT200to400"] = ROOT.kOrange+1
colors_mask["GammaJetsHT400to600"] = ROOT.kOrange+1
colors_mask["GammaJetsHT600toInf"] = ROOT.kOrange+1
colors_mask["GammaJets"]           = ROOT.kOrange+1
colors_mask["WZ"]                  = ROOT.kPink+1
colors_mask["WW"]                  = ROOT.kPink+3
colors_mask["WJetsToLNu0J"]        = ROOT.kCyan-7
colors_mask["WJetsToLNu1J"]        = ROOT.kCyan-8
colors_mask["WJetsToLNu2J"]        = ROOT.kCyan-9
colors_mask["QCD"]                 = ROOT.kRed+1
colors_mask["QCDpT30to50"]         = ROOT.kRed+1
colors_mask["QCDpT50to80"]         = ROOT.kRed+1
colors_mask["QCDpT80to120"]        = ROOT.kRed+1
colors_mask["QCDpT120to170"]       = ROOT.kRed+1
colors_mask["QCDpT170to300"]       = ROOT.kRed+1
colors_mask["QCDpT300toInf"]       = ROOT.kRed+1
colors_mask["DiPhotonJets"]        = ROOT.kGreen+1
colors_mask["ZGammaToLLGamma"]     = ROOT.kBlue+1


# leg1 = ROOT.TLegend(0.15,0.6120093,0.34,0.9491917) #left positioning
leg1 = ROOT.TLegend(0.5868687,0.5520093,0.9411784,0.9091917) #right positioning
leg1.SetHeader("") 
leg1.SetFillColor(0)
leg1.SetBorderSize(1)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)

#LOOP ON INPUT FILES (MC backgrounds, MC Signal, DATA)
for filename in list_inputfiles:
    fileIn = ROOT.TFile(filename)

    sample_name = (filename.split("_")[3])[:-5] #retireve the samplename (QCD, DY, Data, etc...)
    if(debug):
        print "========="+sample_name+"========"

    #LOOP ON HISTOS (for each sample)
    for histo_name in list_histos:
        histo = fileIn.Get(histo_name)

        #REBIN
        if histo_name   == "h_phi_InvMass_TwoTrk":
            histo.Rebin(5) 
        elif histo_name == "h_couple_Iso":
            histo.Rebin(5) 
        elif histo_name == "h_couple_Iso_ch":
            histo.Rebin(5) 
        elif histo_name == "h_nJets_25" or histo_name == "h_nMuons" or histo_name == "h_nElectrons" or histo_name == "h_nPhotons" :
            histo.Rebin(1) 
        elif histo_name == "h_photon_eta" or histo_name == "h_firstKCand_Eta" or histo_name == "h_secondKCand_Eta" or histo_name == "h_bestJetEta" or histo_name == "h_bestCoupleEta" :
            histo.Rebin(4)
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
            histo_container[-1].SetLineColor(4)   #blue
            histo_container[-1].SetLineWidth(4)   #kind of thick
            hsignal[histo_name] = histo_container[-1]
        if "Data" in sample_name:
            histo_container[-1].SetMarkerStyle(20)   #point
            hdata[histo_name] = histo_container[-1]
        if not sample_name == "Data" and not sample_name == "Signal":
            histo_container[-1].SetFillColor(colors_mask[sample_name])
            hstack[histo_name].Add(histo_container[-1])
            
    #fill legend
    #if histo_name == "nMuons": #Add the legend only once (nMuons is just a random variable)

    if histo.Integral() > float(signal_magnify)/500. or sample_name == "Signal": #Only plot in the legend those samples which have some contribution
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
    if pad2Flag:
        pad1 = ROOT.TPad("pad_" + histo_name,"",0.,0.28,1,1)
        pad2 = ROOT.TPad("pad_" + histo_name,'',0,0,1,0.25)
        pad1.SetBottomMargin(0.02)
        pad1.SetBorderMode(0)
        pad1.SetTicks(2,1) #ticks on the right and the upper axis are drawn inside
        pad2.SetTopMargin(0.01)
        pad2.SetBottomMargin(0.3)
        pad2.SetBorderMode(0)
        pad1.Draw()
        pad2.Draw()
        if histo_name == "h_nJets_25" or histo_name == "h_nMuons" or histo_name == "h_nElectrons": # or histo_name == "nPhotons"
            pad1.SetLogy()
        
    ##########################################

    ##########################################
    if pad2Flag:
        pad1.cd()
    ##########################################

    hstack[histo_name].SetTitle(histo_name)
    hstack[histo_name].Draw("histo")

    ##########################################
    #hstack[histo_name].GetXaxis().SetTickLength(0)
    if pad2Flag:
        hstack[histo_name].GetXaxis().SetLabelOffset(999)
    ##########################################

    if histo_name == "h_InvMass_TwoTrk_Photon":
        xMin_mKKg = 100.
        xMax_mKKg = 170.
        yMin_mKKg = 0.
        if pad2Flag:
            yMax_mKKg = 9000.
        else:
            yMax_mKKg = 400.

        hstack[histo_name].GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_mKKg,xMax_mKKg)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_mKKg))

    if histo_name == "h_nJets_25":
        hstack[histo_name].GetXaxis().SetTitle("nJets")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),2500.))
        hstack[histo_name].GetXaxis().SetRangeUser(0.,6.)

    if histo_name == "h_InvMass_TwoTrk_Photon_NoPhiMassCut" :
        hstack[histo_name].GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
        hstack[histo_name].GetXaxis().SetRangeUser(100.,170.)
    
    if histo_name == "h_phi_InvMass_TwoTrk" :
        hstack[histo_name].GetXaxis().SetTitle("m_{K^{+}K^{-}} [GeV]")
        xMin_mKK = 1.
        xMax_mKK = 1.04
        yMin_mKK = 0.
        if pad2Flag:
            yMax_mKK = 16000.
        else:
            yMax_mKK = 600.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_mKK,xMax_mKK)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_mKK))

    if histo_name == "h_firstKCand_pT" :
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{K_{1}} [GeV]")
        xMin_pTK1 = 20.
        xMax_pTK1 = 65.
        yMin_pTK1 = 0.
        if pad2Flag:
            yMax_pTK1 = 80000.
        else:
            yMax_pTK1 = 1000.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_pTK1,xMax_pTK1)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_pTK1))

    if histo_name == "h_secondKCand_pT" :
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{K_2} [GeV]")
        xMin_pTK2 = 12.
        xMax_pTK2 = 58.
        yMin_pTK2 = 0.
        if pad2Flag:
            yMax_pTK2 = 80000.
        else:
            yMax_pTK2 = 1000.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_pTK2,xMax_pTK2)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_pTK2))

    if histo_name == "h_firstKCand_Eta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{K_1}")
        xMin_etaK1 = -2.5
        xMax_etaK1 = 2.5
        yMin_etaK1 = 0.
        if pad2Flag:
            yMax_etaK1 = 27000.
        else:
            yMax_etaK1 = 600.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_etaK1,xMax_etaK1)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_etaK1))


    if histo_name == "h_secondKCand_Eta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{K_2}")
        xMin_etaK2 = -2.5
        xMax_etaK2 = 2.5
        yMin_etaK2 = 0.
        if pad2Flag:
            yMax_etaK2 = 27000.
        else:
            yMax_etaK2 = 600.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_etaK2,xMax_etaK2)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_etaK2))

    if histo_name == "h_firstKCand_Phi" :
        hstack[histo_name].GetXaxis().SetTitle("#phi [rad]")

    if histo_name == "h_secondKCand_Phi" :
        hstack[histo_name].GetXaxis().SetTitle("#phi [rad]")

    if histo_name == "h_bestCouplePt" :
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{K^{+}K^{-}} [GeV]")
        yMin_pTKK = 0.
        if pad2Flag:
            xMin_pTKK = 30.
            xMax_pTKK = 130.
            yMax_pTKK = 50000.        
        else:
            xMin_pTKK = 45.
            xMax_pTKK = 120.
            yMax_pTKK = 450.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_pTKK,xMax_pTKK)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_pTKK))

    if histo_name == "h_bestCoupleEta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{KK}")
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        xMin_etaKK = -2.5
        xMax_etaKK = 2.5
        yMin_etaKK = 0.
        if pad2Flag:
            yMax_etaKK = 27000.
        else:
            yMax_etaKK = 600.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_etaKK,xMax_etaKK)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_etaKK))


    if histo_name == "h_bestJetPt" :
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{jet} [GeV]")
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        xMin_jetPt = 45.
        xMax_jetPt = 250.
        yMin_jetPt = 0.
        if pad2Flag:
            yMax_jetPt = 42000.
        else:
            yMax_jetPt = 600.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_jetPt,xMax_jetPt)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_jetPt))

    if histo_name == "h_bestJetEta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{jet}")
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        xMin_etaJet = -2.5
        xMax_etaJet = 2.5
        yMin_etaJet = 0.
        if pad2Flag:
            yMax_etaJet = 27000.
        else:
            yMax_etaJet = 600.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_etaJet,xMax_etaJet)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_etaJet))


    if histo_name == "h_K1_Iso" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{1}}/pT_{K_{1}}")

    if histo_name == "h_K2_Iso" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")

    if histo_name == "h_K1_Iso_ch" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{1}}/pT_{K_{1}}")

    if histo_name == "h_K2_Iso_ch" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")

    if histo_name == "h_couple_Iso" :
        hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")

    if histo_name == "h_couple_Iso_ch" :
        hstack[histo_name].GetXaxis().SetTitle("#Sigmap_{T}^{ch}/p_{T}^{K^{+}K^{-}}")
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        xMin_isoKK = 0.
        xMax_isoKK = 1.
        yMin_isoKK = 0.
        if pad2Flag:
            yMax_isoKK = 25000.
        else:
            yMax_isoKK = 1000.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_isoKK,xMax_isoKK)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_isoKK))

    if histo_name == "h_bestCoupleDeltaR" :
        hstack[histo_name].GetXaxis().SetTitle("#DeltaR_{K^{+}K^{-}}")

    if histo_name == "h_nPhotons" :
        hstack[histo_name].GetXaxis().SetTitle("n.#gamma")
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        yMin_etaJet = 0.
        if pad2Flag:
            yMax_etaJet = 220000.
        else:
            yMax_etaJet = 2500.
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_etaJet))

    if histo_name == "h_photon_energy" :
        hstack[histo_name].GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        xMin_eTgamma = 30.
        xMax_eTgamma = 160.
        yMin_eTgamma = 0.
        if pad2Flag:
            yMax_eTgamma = 45000.
        else:
            yMax_eTgamma = 350.
        hstack[histo_name].GetXaxis().SetRangeUser(xMin_eTgamma,xMax_eTgamma)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_eTgamma))

    if histo_name == "h_photon_eta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{#gamma}")
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        yMin_etagamma = 0.
        if pad2Flag:
            yMax_etagamma = 30000.
        else:
            yMax_etagamma = 650.
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),yMax_etagamma))

    if histo_name == "h_nElectrons" :
        hstack[histo_name].GetXaxis().SetTitle("n.electrons")
#        hstack[histo_name].GetXaxis().SetRangeUser(0.,3.)
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),1000000.))

    if histo_name == "h_nMuons" :
        hstack[histo_name].GetXaxis().SetTitle("n.muons")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),10000000.))



    if signal_magnify != 1:
        hsignal[histo_name].Scale(signal_magnify)      


    hstack[histo_name].Draw("histo")
    hstack[histo_name].GetYaxis().SetTitleOffset(1.5);
    hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
    hstack[histo_name].GetYaxis().SetTitle("Events")
    hsignal[histo_name].Draw("SAME,hist")
    hdata[histo_name].Draw("SAME,E1")

    hMCErr = copy.deepcopy(hstack[histo_name].GetStack().Last())
    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(1)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.Draw("sameE2")

    leg1.Draw()

    ################################################
    if pad2Flag:
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
        #    totalData.GetYaxis().SetRangeUser(0.4,1.6)
        
        totalData.GetXaxis().SetLabelSize(0.10)
        totalData.GetXaxis().SetTitleSize(0.12)
        totalData.GetXaxis().SetTitleOffset(1.0)
        totalData.SetMaximum(2.2)
        totalData.SetMinimum(-0.2)


        line_on_one = ROOT.TLine(hstack[histo_name].GetXaxis().GetXmin(),1.,hstack[histo_name].GetXaxis().GetXmax(),1.)
        line_on_one.SetLineColor(38)

        if histo_name == "h_InvMass_TwoTrk_Photon":
            totalData.GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
            totalData.GetXaxis().SetRangeUser(xMin_mKKg,xMax_mKKg)
            line_on_one = ROOT.TLine(xMin_mKKg,1.,xMax_mKKg,1.)
            line_on_one.SetLineColor(38)
        if histo_name == "h_phi_InvMass_TwoTrk":
            totalData.GetXaxis().SetTitle("m_{K^{+}K^{-}} [GeV]")
            totalData.GetXaxis().SetRangeUser(xMin_mKK,xMax_mKK)
            line_on_one = ROOT.TLine(xMin_mKK,1.,xMax_mKK,1.)
            line_on_one.SetLineColor(38)
        if histo_name == "h_nJets_25" :
            totalData.GetXaxis().SetTitle("nJets")
        if histo_name == "h_photon_energy":
            totalData.GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
            totalData.GetXaxis().SetRangeUser(xMin_eTgamma,xMax_eTgamma)
            line_on_one = ROOT.TLine(xMin_eTgamma,1.,xMax_eTgamma,1.)
            line_on_one.SetLineColor(38)
        if histo_name == "h_bestCouplePt" :
            totalData.GetXaxis().SetTitle("p_{T}^{K^{+}K^{-}} [GeV]")
            totalData.GetXaxis().SetRangeUser(xMin_pTKK,xMax_pTKK)
            line_on_one = ROOT.TLine(xMin_pTKK,1.,xMax_pTKK,1.)
            line_on_one.SetLineColor(38)
        if histo_name == "h_bestJetPt" :
            totalData.GetXaxis().SetTitle("p_{T}^{jet} [GeV]")
            totalData.GetXaxis().SetRangeUser(xMin_jetPt,xMax_jetPt)
            line_on_one = ROOT.TLine(xMin_jetPt,1.,xMax_jetPt,1.)
            line_on_one.SetLineColor(38)
        if histo_name == "h_couple_Iso_ch" :
            totalData.GetXaxis().SetTitle("#Sigmap_{T}^{ch}/p_{T}^{K^{+}K^{-}}")
        if histo_name == "h_firstKCand_pT" :
            totalData.GetXaxis().SetTitle("p_{T}^{K_{1}} [GeV]")
            totalData.GetXaxis().SetRangeUser(xMin_pTK1,xMax_pTK1)
            line_on_one = ROOT.TLine(xMin_pTK1,1.,xMax_pTK1,1.)
            line_on_one.SetLineColor(38)
        if histo_name == "h_secondKCand_pT" :
            totalData.GetXaxis().SetTitle("p_{T}^{K_{2}} [GeV]")
            totalData.GetXaxis().SetRangeUser(xMin_pTK2,xMax_pTK2)
            line_on_one = ROOT.TLine(xMin_pTK2,1.,xMax_pTK2,1.)
            line_on_one.SetLineColor(38)
        if histo_name == "h_nElectrons" :
            totalData.GetXaxis().SetTitle("n.electrons")
        if histo_name == "h_nMuons" :
            totalData.GetXaxis().SetTitle("n.muons")



        
        totalData.Draw()
        line_on_one.Draw("SAME")
        
        ################################################

    canvas[histo_name].SaveAs(output_dir + histo_name + ".pdf")
    canvas[histo_name].SaveAs(output_dir + histo_name + ".png")
