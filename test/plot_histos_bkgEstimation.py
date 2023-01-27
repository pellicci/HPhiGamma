import ROOT
import math
import copy
import sys
import tdrstyle, CMS_lumi
from ROOT import gROOT

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

signal_magnify = float(sys.argv[1])
CR_magnify = 1. #2079./1179.

plotOnlyData = False
isTightSelection = int(sys.argv[2])

isPhi = int(sys.argv[3]) #note that in python true = 1 and false = 0


if not isTightSelection:
    inputnames = ["Data","SignalggH","SignalVBF","SidebandsNorm"]
else:
    inputnames = ["Data","SignalggH","SignalVBF"]

list_inputfiles = []
for filename in sys.argv[4:]:
    list_inputfiles.append(filename)

fileIn = ROOT.TFile(list_inputfiles[0])
CAT = (filename.split("_")[3]) 

print "#############################"
print "is Phi = ",isPhi
print "Category = ", CAT
print "#############################"

#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.9
CMS_lumi.cmsTextSize = 1.5
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}"

hstack     = dict()
#hsignal    = dict()
hsignalggH = dict()
hsignalVBF = dict()
hdata      = dict()
canvas     = dict()
histo_container = [] #just for memory management

#Get the list of histograms
signalfile = ROOT.TFile("histos/latest_production/histos_SR_"+CAT+"_SignalggH.root")
list_histos = []
keylist = signalfile.GetListOfKeys()
key = ROOT.TKey()
for key in keylist :
    obj_class = ROOT.gROOT.GetClass(key.GetClassName())
    if not obj_class.InheritsFrom("TH1") :
        continue
    if not (key.ReadObj().GetName() == "h_efficiency" or key.ReadObj().GetName() == "h_cutOverflow"): #h_efficiency and h_cutOverflow is a plot plotted in other way   
        list_histos.append( key.ReadObj().GetName() )

for hname in list_histos:
    hstack[hname] = ROOT.THStack("hstack_" + hname,"")

#COLOR MASK
colors_mask = dict()
#colors_mask["bkgEstimationCR"]   = ROOT.kRed-7
#colors_mask["Sidebands"]          = ROOT.kRed-7
if isPhi:
    colors_mask["SidebandsNorm"]      = ROOT.kCyan-7
    decayChannel = "#phi#gamma "
else:
    colors_mask["SidebandsNorm"]      = ROOT.kRed-4
    decayChannel = "#rho#gamma "

colors_mask["GammaJets"] = ROOT.kOrange
colors_mask["QCD"]       = ROOT.kRed



#LEGEND
#leg1 = ROOT.TLegend(0.6868687,0.6120093,0.9511784,0.9491917) #right positioning
#leg1.SetHeader(" ")
#leg1.SetFillColor(0)
#leg1.SetBorderSize(0)
#leg1.SetLineColor(1)
#leg1.SetLineStyle(1)
#leg1.SetLineWidth(1)
#leg1.SetFillStyle(1001)


leg1 = ROOT.TLegend(0.65,0.62,0.95,0.95) #right positioning
#leg1 = ROOT.TLegend(0.321,0.58,0.981,0.95) #right positioning
leg1.SetHeader(" ")
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.SetNColumns(1)

for filename in list_inputfiles:
    fileIn = ROOT.TFile(filename)
    sample_name = (filename.split("_")[4])[:-5] 
    print "=============== ", sample_name
    for histo_name in list_histos:
        histo = fileIn.Get(histo_name)
        print "histo_name = ",histo_name
        # Set to 0 the bins containing negative values, due to negative weights
        hsize = histo.GetSize() - 2 # GetSize() returns the number of bins +2 (that is + overflow + underflow) 
        for bin in range(1,hsize+1): # The +1 is in order to get the last bin
            bincontent = histo.GetBinContent(bin)
            if bincontent < 0.:
                histo.SetBinContent(bin,0.)

        histo_container.append(copy.copy(histo))
        
        if not histo_name == "h_nMuons" and not histo_name == "h_nPhotons" and not histo_name == "h_nJets_25" and not histo_name == "h_nElectrons":
            if isTightSelection and not (histo_name == "h_firstTrk_Iso" or histo_name == "h_firstTrk_Iso_ch" or histo_name == "h_firstTrk_Iso_neutral" or histo_name == "h_secondTrk_Iso" or histo_name == "h_secondTrk_Iso_ch" or histo_name == "h_couple_AbsIsoCh" or histo_name == "h_couple_Iso" or histo_name == "h_couple_Iso_ch"):
                histo_container[-1].Rebin(10)
            elif histo_name == "h_InvMass_TwoTrk_Photon":
                histo_container[-1].Rebin(10)
                #hsignalVBF[histo_name].Rebin(1/5)
                #hsignalggH[histo_name].Rebin(1/5)

            else:
                histo_container[-1].Rebin(5)

        #if sample_name == "Signal":
         #   histo_container[-1].SetLineStyle(1)   
          #  histo_container[-1].SetLineColor(1)   #4 for blue, 2 for red
          #  histo_container[-1].SetLineWidth(4)   #kind of thick
          #  hsignal[histo_name] = histo_container[-1]
        if sample_name == "SignalggH":
            histo_container[-1].SetLineStyle(1)   
            histo_container[-1].SetLineColor(4)   #4 for blue, 2 for red
            histo_container[-1].SetLineWidth(4)   #kind of thick
            hsignalggH[histo_name] = histo_container[-1]
        elif sample_name == "SignalVBF":
            histo_container[-1].SetLineStyle(1)   
            histo_container[-1].SetLineColor(8)   #4 for blue, 2 for red
            histo_container[-1].SetLineWidth(4)   #kind of thick
            hsignalVBF[histo_name] = histo_container[-1]
        elif "Data" in sample_name:
            histo_container[-1].SetMarkerStyle(20)   #point
            histo_container[-1].SetBinErrorOption(ROOT.TH1.kPoisson)
            hdata[histo_name] = histo_container[-1]
        else:
            histo_container[-1].SetFillColor(colors_mask[sample_name])
            histo_container[-1].SetLineColor(colors_mask[sample_name])
            histo_container[-1].SetBinErrorOption(ROOT.TH1.kPoisson)
            hstack[histo_name].Add(histo_container[-1])

        if plotOnlyData :
            hstack[histo_name].Add(histo_container[-1])


        if histo_name == "h_InvMass_TwoTrk_Photon" : #Add the legend only once (InvMass_TwoTrk_Photon is just a random variable)
            
            if not sample_name == "Data" and not sample_name == "Signal" and not sample_name == "SignalVBF" and not sample_name == "SignalggH":
                leg1.AddEntry(histo_container[-1],"bkg estimation","f")
            elif sample_name == "Data":
                leg1.AddEntry(histo_container[-1],sample_name,"ep")
            elif sample_name == "SignalggH":
                leg1.AddEntry(histo_container[-1],decayChannel + "ggH" ,"f")
            elif sample_name == "SignalVBF":
                leg1.AddEntry(histo_container[-1],decayChannel + "VBF" ,"f")
            elif sample_name == "Signal":
                leg1.AddEntry(histo_container[-1],decayChannel + "ggH + VBF","f")
   
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
    #hsignal[histo_name].SetTitle("")
    hsignalggH[histo_name].SetTitle("")
    hsignalVBF[histo_name].SetTitle("")


    if not plotOnlyData :

        hstack[histo_name].Draw("histo")
        hstack[histo_name].GetYaxis().SetTitleSize(0.07)
        hstack[histo_name].GetYaxis().SetTitleOffset(0.7)
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        hstack[histo_name].GetXaxis().SetLabelOffset(999)
        hstack[histo_name].GetYaxis().SetMaxDigits(3)

        if histo_name == "h_nJets_25" or histo_name == "h_nMuons" or histo_name == "h_nElectrons" or histo_name == "h_nPhotons":
            hstack[histo_name].SetMaximum(100 * hdata[histo_name].GetMaximum())
            hstack[histo_name].SetMinimum(1.)  #cannot use values < 1 otherwise log is negative
        else:
            hstack[histo_name].SetMaximum(2.1 * hdata[histo_name].GetMaximum())


        if histo_name == "h_InvMass_TwoTrk_Photon":
            #hstack[histo_name].Rebin(2)            
            hstack[histo_name].GetXaxis().SetTitle("m_{ditrk#gamma} [GeV]")
            hstack[histo_name].GetXaxis().SetLimits(100.,170.)


        if histo_name == "h_nJets_25":
            hstack[histo_name].GetXaxis().SetTitle("nJets")
            hstack[histo_name].GetXaxis().SetLimits(-0.5,6.5)

        if histo_name == "h_meson_InvMass_TwoTrk" :
            hstack[histo_name].GetXaxis().SetTitle("m_{ditrk} [GeV]")
            if isPhi:
                hstack[histo_name].GetXaxis().SetLimits(1.00,1.042)
            else:
                hstack[histo_name].GetXaxis().SetLimits(0.501,0.999)

        if histo_name == "h_firstTrk_pT" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{Trk_{1}} [GeV]")
            #hstack[histo_name].GetXaxis().SetLimits(15.,60.)

        if histo_name == "h_secondTrk_pT" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{Trk_{2}} [GeV]")
            #hstack[histo_name].GetXaxis().SetLimits(5.,50.)

        if histo_name == "h_firstTrk_Eta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{Trk_{1}}")
            hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

        if histo_name == "h_secondTrk_Eta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{Trk_{2}}")
            hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

        if histo_name == "h_bestCoupleEta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{ditrk}")
            hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

        if histo_name == "h_firstTrk_Phi" :
            hstack[histo_name].GetXaxis().SetTitle("#phi_{Trk_{1}} [rad]")

        if histo_name == "h_secondTrk_Phi" :
            hstack[histo_name].GetXaxis().SetTitle("#phi_{Trk_{2}} [rad]")

        if histo_name == "h_bestCouplePt" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{ditrk} [GeV]")

        if histo_name == "h_bestJetPt" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{jet} [GeV]")
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(40.,140.)

        if histo_name == "h_bestJetEta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{jet}")
            hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

        if histo_name == "h_firstTrk_Iso" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{Trk_{1}}/(#Sigmap_{T} + p_{T}^{Trk_{1}})")
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(0.6,1.)

        if histo_name == "h_secondTrk_Iso" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{Trk_{2}}/(#Sigmap_{T} + p_{T}^{Trk_{2}})")
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(0.5,1.)

        if histo_name == "h_firstTrk_Iso_ch" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{Trk_{1}}/(#Sigmap_{T}^{ch} + p_{T}^{Trk_{1}})")      
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(0.94,1.)

        if histo_name == "h_secondTrk_Iso_ch" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{Trk_{2}}/(#Sigmap_{T}^{ch} + p_{T}^{Trk_{2}})")
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(0.92,1.)

        if histo_name == "h_couple_Iso" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{ditrk}/(#Sigmap_{T} + p_{T}^{ditrk})")
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(0.7,1.)

        if histo_name == "h_couple_Iso_ch" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{ditrk}/(#Sigmap_{T}^{ch} + p_{T}^{ditrk})")
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(0.97,1.)

        if histo_name == "h_couple_Iso_neutral" :
            hstack[histo_name].GetXaxis().SetTitle("p_{T}^{ditrk}/(#Sigmap_{T}^{0} + p_{T}^{ditrk})")
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(0.7,1.)

        if histo_name == "h_bestCoupleDeltaR" :
            hstack[histo_name].GetXaxis().SetTitle("#DeltaR_{ditrk}/m_{ditrk}")

        if histo_name == "h_nPhotons" :
            hstack[histo_name].GetXaxis().SetTitle("n.#gamma")
            hstack[histo_name].GetXaxis().SetLimits(-0.5,4.5)

        if histo_name == "h_photon_energy" :
            hstack[histo_name].GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
            hstack[histo_name].GetXaxis().SetLimits(38.,160.)
            if isTightSelection: hstack[histo_name].GetXaxis().SetLimits(38.,120.)

        if histo_name == "h_photon_eta" :
            hstack[histo_name].GetXaxis().SetTitle("#eta_{#gamma}")

        if histo_name == "h_nElectrons" :
            hstack[histo_name].GetXaxis().SetTitle("n.electrons")
            hstack[histo_name].GetXaxis().SetLimits(-0.5,2.5)

        if histo_name == "h_nMuons" :
            hstack[histo_name].GetXaxis().SetTitle("n.muons")
            hstack[histo_name].GetXaxis().SetLimits(-0.5,3.5)

        if histo_name == "h_decayChannel":
            hstack[histo_name].GetXaxis().SetTitle("decay channel")

        if histo_name == "h_dPhiGammaTrk":
            hstack[histo_name].GetXaxis().SetTitle("#Delta#phi_{#gamma,Trk_{1}} [rad]")
                


        hstack[histo_name].Draw("SAME,histo")

    if plotOnlyData:
        hstack[histo_name].Draw("")
        hstack[histo_name].GetYaxis().SetTitleSize(0.07)
        hstack[histo_name].GetYaxis().SetTitleOffset(0.7)
        hstack[histo_name].GetYaxis().SetTitle("Events")
        hstack[histo_name].SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
        #hstack[histo_name].GetXaxis().SetLabelOffset(999)
        hstack[histo_name].GetYaxis().SetMaxDigits(2)

        if histo_name == "h_meson_InvMass_TwoTrk" :
            hstack[histo_name].SetMaximum(2.1 * hdata[histo_name].GetMaximum())
            if isPhi:
                hstack[histo_name].GetXaxis().SetLimits(1.,1.04)
                hstack[histo_name].GetXaxis().SetTitle("m_{KK} [GeV]")
            else:
                hstack[histo_name].GetXaxis().SetLimits(0.5,1.)
                hstack[histo_name].GetXaxis().SetTitle("m_{#pi#pi} [GeV]")

    if signal_magnify != 1:
        #hsignal[histo_name].Scale(signal_magnify)   
        hsignalggH[histo_name].Scale(signal_magnify)   
        hsignalVBF[histo_name].Scale(signal_magnify)   
        if histo_name == "h_InvMass_TwoTrk_Photon":
            hsignalVBF[histo_name].Rebin(1/5) #do this to have more granularity in the signal curve
            hsignalggH[histo_name].Rebin(1/5)

    hdata[histo_name].Draw("SAME,E1,X0")
    #hsignal[histo_name].Draw("SAME,hist")
    hsignalggH[histo_name].Draw("SAME,hist")
    hsignalVBF[histo_name].Draw("SAME,hist")


    hMCErr = copy.deepcopy(hstack[histo_name].GetStack().Last())
    hMCErr_size = hMCErr.GetSize() - 2
    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(0)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.SetLineColor(0)
    if not plotOnlyData :
        hMCErr.Draw("sameE2")

    leg1.Draw()

    if not plotOnlyData:
        CMS_lumi.CMS_lumi(pad1, iPeriod, iPos) #Print integrated lumi and energy information
    else: 
        CMS_lumi.CMS_lumi(canvas[histo_name], iPeriod, iPos) #Print integrated lumi and energy information


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
            #else:
             #   new_MC_BinError = 0.

            #Set data/MC ratio points error bar to data relative uncertainty
            if not totalData_forErrors.GetBinContent(bin) == 0:
                new_Data_BinError = totalData_forErrors.GetBinError(bin)/totalData_forErrors.GetBinContent(bin)
            #else:
             #   new_Data_BinError = 0.

            totalMC.SetBinError(bin,new_MC_BinError)
            totalMC.SetBinContent(bin,1.)
            totalData.SetBinError(bin,new_Data_BinError)
    
        totalData.SetTitle("")
        totalData.SetMarkerStyle(8)
        totalData.SetMarkerColor(1)
        totalData.SetLineColor(1)
        totalData.GetYaxis().SetLabelSize(0.1)
        totalData.GetYaxis().SetTitle("Data/Bkg")
        totalData.GetYaxis().SetTitleSize(0.16)
        totalData.GetYaxis().SetTitleOffset(0.3)
        totalData.GetYaxis().SetRangeUser(0.8,1.2)
        if isPhi: totalData.GetYaxis().SetRangeUser(0.5,1.5)
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
    output_dir = "~/cernbox/www/latest_production/"+CAT+"_latest_production/"


    canvas[histo_name].SaveAs(output_dir + histo_name + ".pdf")
    canvas[histo_name].SaveAs(output_dir + histo_name + ".png")
