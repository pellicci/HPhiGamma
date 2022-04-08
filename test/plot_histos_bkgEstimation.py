import ROOT
import math
import copy
import sys
import tdrstyle, CMS_lumi

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

signal_magnify = int(sys.argv[1])
CR_magnify = 1. #2079./1179.

plotOnlyData = False
isTightSelection = int(sys.argv[2])

inputnames = ["Signal","Data","SidebandsNorm"]

list_inputfiles = []
for filename in sys.argv[3:]:
    list_inputfiles.append(filename)

#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.9
CMS_lumi.cmsTextSize = 1.5
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}"

hstack  = dict()
hsignal = dict()
hdata   = dict()
canvas  = dict()
histo_container = [] #just for memory management

#Get the list of histograms
list_histos = []
signalfile = ROOT.TFile("histos/latest_production/histos_SR_Signal.root")
keylist = signalfile.GetListOfKeys()
key = ROOT.TKey()
for key in keylist :
    obj_class = ROOT.gROOT.GetClass(key.GetClassName())
    if not obj_class.InheritsFrom("TH1") :
        continue
    if not key.ReadObj().GetName() == "h_efficiency": #h_efficiency is a plot plotted in other way   
        list_histos.append( key.ReadObj().GetName() )

for hname in list_histos:
    hstack[hname] = ROOT.THStack("hstack_" + hname,"")

#COLOR MASK
colors_mask = dict()
#colors_mask["bkgEstimationCR"]   = ROOT.kRed-7
#colors_mask["Sidebands"]          = ROOT.kRed-7
colors_mask["SidebandsNorm"]      = ROOT.kTeal-9


#LEGEND
#leg1 = ROOT.TLegend(0.6868687,0.6120093,0.9511784,0.9491917) #right positioning
#leg1.SetHeader(" ")
#leg1.SetFillColor(0)
#leg1.SetBorderSize(0)
#leg1.SetLineColor(1)
#leg1.SetLineStyle(1)
#leg1.SetLineWidth(1)
#leg1.SetFillStyle(1001)


leg1 = ROOT.TLegend(0.45,0.62,0.85,0.95) #right positioning
#leg1 = ROOT.TLegend(0.321,0.58,0.981,0.95) #right positioning
leg1.SetHeader(" ")
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.SetNColumns(2)


for filename in list_inputfiles:
    fileIn = ROOT.TFile(filename)

    sample_name = (filename.split("_")[3])[:-5] 
    for histo_name in list_histos:
        histo = fileIn.Get(histo_name)

        # Set to 0 the bins containing negative values, due to negative weights
        hsize = histo.GetSize() - 2 # GetSize() returns the number of bins +2 (that is + overflow + underflow) 
        for bin in range(1,hsize+1): # The +1 is in order to get the last bin
            bincontent = histo.GetBinContent(bin)
            if bincontent < 0.:
                histo.SetBinContent(bin,0.)

        histo_container.append(copy.copy(histo))
        
        if not histo_name == "h_nMuons" and not histo_name == "h_nPhotons" and not histo_name == "h_nJets_25" and not histo_name == "h_nElectrons" and not histo_name == "h_photonWP90":
            print histo_name
            if isTightSelection and (histo_name == "h_K1_Iso" or histo_name == "h_K1_Iso_ch" or histo_name == "h_K2_Iso" or histo_name == "h_K2_Iso_ch" or histo_name == "h_couple_AbsIsoCh" or histo_name == "h_couple_Iso" or histo_name == "h_couple_Iso_ch"):
                histo_container[-1].Rebin(2)
            else:
                histo_container[-1].Rebin(5)

        if "Signal" in sample_name:
            histo_container[-1].SetLineStyle(2)   #dashed
            histo_container[-1].SetLineColor(4)   #blue, 2 for red
            histo_container[-1].SetLineWidth(4)   #kind of thick
            hsignal[histo_name] = histo_container[-1]
        elif "Data" in sample_name:
            histo_container[-1].SetMarkerStyle(20)   #point
            hdata[histo_name] = histo_container[-1]
        else:
            histo_container[-1].SetFillColor(colors_mask[sample_name])
            histo_container[-1].SetLineColor(colors_mask[sample_name])
            hstack[histo_name].Add(histo_container[-1])

        if plotOnlyData :
            hstack[histo_name].Add(histo_container[-1])


        if histo_name == "h_InvMass_TwoTrk_Photon" : #Add the legend only once (InvMass_TwoTrk_Photon is just a random variable)

            if not sample_name == "Data" and not sample_name == "Signal":
                leg1.AddEntry(histo_container[-1],sample_name,"f")
            elif sample_name == "Data":
                leg1.AddEntry(histo_container[-1],sample_name,"ep")
            elif sample_name == "Signal":
                leg1.AddEntry(histo_container[-1],sample_name + " x " + str(signal_magnify),"f")

    fileIn.Close()

for histo_name in list_histos:

    canvas[histo_name] = ROOT.TCanvas("Canvas_" + histo_name,"",200,106,600,600)
    canvas[histo_name].cd()
 
    if not plotOnlyData :   ##########################################
        pad1 = ROOT.TPad("pad_" + histo_name,"",0,0.28,1,1.)
        pad2 = ROOT.TPad("pad_" + histo_name,'',0,0.01,1,0.27)
        #pad1.SetTopMargin(0.047)
        pad1.SetBottomMargin(0.02)
        pad1.SetBorderMode(0)
        pad1.SetBorderSize(0)
        pad1.SetFrameBorderSize(0)
        pad1.SetTicks(2,1) #ticks on the right and the upper axis are drawn inside
        pad2.SetBorderSize(0)
        pad2.SetFrameBorderSize(0)
        pad2.SetBottomMargin(0.3)
        pad2.SetBorderMode(0)
        pad1.Draw()
        pad2.Draw()
        if histo_name == "h_nJets_25" or histo_name == "h_nMuons" or histo_name == "h_nElectrons" or histo_name == "h_nPhotons":
            pad1.SetLogy()
        
        pad1.cd()

    hstack[histo_name].SetTitle("")
    hdata[histo_name].SetTitle("")
    hsignal[histo_name].SetTitle("")

    if not plotOnlyData :

        hstack[histo_name].Draw("histo")
        hstack[histo_name].GetYaxis().SetTitleSize(0.07)
        hstack[histo_name].GetYaxis().SetTitleOffset(0.7)
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        hstack[histo_name].GetXaxis().SetLabelOffset(999)
        hstack[histo_name].GetYaxis().SetMaxDigits(3)

        if histo_name == "h_nJets_25" or histo_name == "h_nMuons" or histo_name == "h_nElectrons" or histo_name == "h_nPhotons":
            hstack[histo_name].SetMaximum(10 * hdata[histo_name].GetMaximum())
        else:
            hstack[histo_name].SetMaximum(2.1 * hdata[histo_name].GetMaximum())


        if histo_name == "h_InvMass_TwoTrk_Photon":            
            hstack[histo_name].GetXaxis().SetTitle("m_{K^{+}K^{-}#gamma} [GeV]")
            hstack[histo_name].GetXaxis().SetLimits(100.,150.)

        if histo_name == "h_nJets_25":
            hstack[histo_name].GetXaxis().SetTitle("nJets")
            hstack[histo_name].GetXaxis().SetLimits(-0.5,6.5)
    
        if histo_name == "h_phi_InvMass_TwoTrk" :
            hstack[histo_name].GetXaxis().SetTitle("m_{K^{+}K^{-}} [GeV]")
            hstack[histo_name].GetXaxis().SetLimits(1.006,1.034)

        if histo_name == "h_firstKCand_pT" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{K_{1}} [GeV]")
            hstack[histo_name].GetXaxis().SetLimits(20.,60.)

        if histo_name == "h_secondKCand_pT" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{K_{2}} [GeV]")
            hstack[histo_name].GetXaxis().SetLimits(12.,50.)

        if histo_name == "h_firstKCand_Eta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{K_1}")
            hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

        if histo_name == "h_secondKCand_Eta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{K_2}")
            hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

        if histo_name == "h_bestCoupleEta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{KK}")
            hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

        if histo_name == "h_firstKCand_Phi" :
            hstack[histo_name].GetXaxis().SetTitle("#phi [rad]")

        if histo_name == "h_secondKCand_Phi" :
            hstack[histo_name].GetXaxis().SetTitle("#phi [rad]")

        if histo_name == "h_bestCouplePt" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{K^{+}K^{-}} [GeV]")
            hstack[histo_name].GetXaxis().SetLimits(30.,82.)

        if histo_name == "h_bestJetPt" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{jet} [GeV]")
            hstack[histo_name].GetXaxis().SetLimits(45.,250.)
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(45.,220.)

        if histo_name == "h_bestJetEta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{jet}")
            hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

        if histo_name == "h_K1_Iso" :
            hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{1}}/pT_{K_{1}}")
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(0.,0.4)

        if histo_name == "h_K2_Iso" :
            hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(0.,0.52)

        if histo_name == "h_K1_Iso_ch" :
            hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{1}}/pT_{K_{1}}")
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(0.,0.3)        

        if histo_name == "h_K2_Iso_ch" :
            hstack[histo_name].GetXaxis().SetTitle("#SigmapT_{K_{2}}/pT_{K_{2}}")
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(0.,0.42)

        if histo_name == "h_couple_Iso" :
            hstack[histo_name].GetXaxis().SetTitle("#Sigmap_{T}/p_{T}^{K^{+}K^{-}}")
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(0.,0.25)

        if histo_name == "h_couple_Iso_ch" :
            hstack[histo_name].GetXaxis().SetTitle("#Sigmap_{T}^{ch}/p_{T}^{K^{+}K^{-}}")
            hstack[histo_name].GetXaxis().SetLimits(0.,1.)
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(0.,0.22)

        if histo_name == "h_bestCoupleDeltaR" :
            hstack[histo_name].GetXaxis().SetTitle("#DeltaR_{K^{+}K^{-}}")

        if histo_name == "h_nPhotons" :
            hstack[histo_name].GetXaxis().SetTitle("n.#gamma")
            hstack[histo_name].GetXaxis().SetLimits(-0.5,4.5)

        if histo_name == "h_photon_energy" :
            hstack[histo_name].GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
            hstack[histo_name].GetXaxis().SetLimits(30.,160.)
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(30.,150)

        if histo_name == "h_photon_eta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{#gamma}")

        if histo_name == "h_nElectrons" :
            hstack[histo_name].GetXaxis().SetTitle("n.electrons")
            hstack[histo_name].GetXaxis().SetLimits(-0.5,2.5)

        if histo_name == "h_nMuons" :
            hstack[histo_name].GetXaxis().SetTitle("n.muons")
            hstack[histo_name].GetXaxis().SetLimits(-0.5,3.5)

        if histo_name == "h_couple_AbsIsoCh" :
            hstack[histo_name].GetXaxis().SetTitle("absIso_{KK}")
            hstack[histo_name].GetXaxis().SetLimits(0.,30.)
            if isTightSelection:
                hstack[histo_name].GetXaxis().SetLimits(0.,12.)


        hstack[histo_name].Draw("SAME,histo")

    if signal_magnify != 1:
        hsignal[histo_name].Scale(signal_magnify)   
  

    hdata[histo_name].Draw("SAME,E1,X0")
    hsignal[histo_name].Draw("SAME,hist")

    hMCErr = copy.deepcopy(hstack[histo_name].GetStack().Last())
    hMCErr_size = hMCErr.GetSize() - 2
    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(0)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.SetLineColor(0)
    if not plotOnlyData :
        hMCErr.Draw("sameE2")

    if "threegamma" in histo_name and not plotOnlyData :#Add the legend only once
        leg1.AddEntry(hMCErr,"Bkg unc","f")
    leg1.Draw()

    CMS_lumi.CMS_lumi(pad1, iPeriod, iPos) #Print integrated lumi and energy information


    ################################################

    if not plotOnlyData :
        pad2.cd()
        pad2.SetTopMargin(0.03)
        pad2.SetFillColor(0)
        pad2.SetFillStyle(0)
        ROOT.gStyle.SetOptStat(0)
        totalMC = copy.deepcopy(hMCErr)
        totalData = copy.deepcopy(hdata[histo_name])
        totalData_forErrors = copy.deepcopy(hdata[histo_name])
        totalData.Divide(totalMC)

        for bin in range(1,hMCErr_size+1):
        
            #Set MC error band to MC relative uncertainty
            if not totalMC.GetBinContent(bin) == 0:
                new_MC_BinError = totalMC.GetBinError(bin)/totalMC.GetBinContent(bin)
            else:
                new_MC_BinError = 0.

            #Set data/MC ratio points error bar to data relative uncertainty
            if not totalData_forErrors.GetBinContent(bin) == 0:
                new_Data_BinError = totalData_forErrors.GetBinError(bin)/totalData_forErrors.GetBinContent(bin)
            else:
                new_Data_BinError = 0.

            totalMC.SetBinError(bin,new_MC_BinError)
            totalMC.SetBinContent(bin,1.)
            totalData.SetBinError(bin,new_Data_BinError)
    
        totalData.SetTitle("")
        totalData.SetMarkerStyle(8)
        totalData.SetMarkerColor(1)
        totalData.SetLineColor(1)
        totalData.GetYaxis().SetLabelSize(0.1)
        totalData.GetYaxis().SetTitle("Data/MC")
        totalData.GetYaxis().SetTitleSize(0.16)
        totalData.GetYaxis().SetTitleOffset(0.3)
        totalData.GetYaxis().SetRangeUser(0.,2.)
        totalData.GetYaxis().SetNdivisions(502,ROOT.kFALSE)
        totalData.GetXaxis().SetLabelSize(0.10)
        totalData.GetXaxis().SetTitleSize(0.12)
        totalData.GetXaxis().SetTitleOffset(1.0)
        totalData.GetXaxis().SetTitle(hstack[histo_name].GetXaxis().GetTitle())

        totalMC.SetTitle("")
        totalMC.SetFillStyle(3002)

        totalData.GetXaxis().SetRangeUser(hstack[histo_name].GetXaxis().GetXmin(),hstack[histo_name].GetXaxis().GetXmax())
        line_on_one = ROOT.TLine(totalData.GetXaxis().GetXmin(),1.,totalData.GetXaxis().GetXmax(),1.)
        line_on_one.SetLineColor(4)
        line_on_one.SetLineStyle(2)

        totalData.Draw("E1,X0")
        totalMC.Draw("sameE2")
        line_on_one.Draw("SAME")
    ################################################

    output_dir = "~/cernbox/www/testDir/"
    canvas[histo_name].SaveAs(output_dir + histo_name + ".pdf")
    canvas[histo_name].SaveAs(output_dir + histo_name + ".png")
