import ROOT

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+") #import the Doube Crystal Ball PDF


#INPUT SETTINGS
fInput = ROOT.TFile("workspaces/ws_data_blind.root") 
fInput.cd()
workspace = fInput.Get("myworkspace")
workspace.Print()
#VAR and CONSTANTS SETTINGS
mesonGammaMass = workspace.var("mesonGammaMass")
#datasetGenerated = workspace.data("totPDFData")

B_R_ = workspace.var("B_R_")

#workspace.var("a_bkg").setConstant(1)
#workspace.var("b_bkg").setConstant(1)

#workspace.var("dCB_nL").setConstant(1) 
#workspace.var("dCB_nR").setConstant(1)
#workspace.var("dCB_aL").setConstant(1)
#workspace.var("dCB_aR").setConstant(1)
#workspace.var("dCB_pole").setConstant(1)
#workspace.var("dCB_width").setConstant(1)


#RETRIEVE total PDF
totPDF_chebychev   = workspace.pdf("totPDF_chebychev")
totPDF_exponential = workspace.pdf("totPDF_exponential")

#Construct the Toy-MC machinery
#Binned(kTRUE)
print "pdf retrieved"
mcstudy = ROOT.RooMCStudy(totPDF_chebychev, ROOT.RooArgSet(mesonGammaMass), ROOT.RooFit.Silence(), ROOT.RooFit.Extended(1), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))
#mcstudy = ROOT.RooMCStudy(totPDF, ROOT.RooArgSet(mesonGammaMass), ROOT.RooFit.Silence(), ROOT.RooFit.Extended(), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))

#mcstudy = ROOT.RooMCStudy(totPDF, ROOT.RooArgSet(mesonGammaMass), ROOT.RooFit.Silence(), ROOT.RooFit.ProtoData(datasetGenerated), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))
print "mcstudy set"
mcstudy.generateAndFit(100)
print "toy generated"

#Plot the distributions of the fitted parameter, the error and the pull
BRval_frame = mcstudy.plotParam(B_R_, ROOT.RooFit.Bins(28))
print "BRval_frame set!"
BRerr_frame = mcstudy.plotError(B_R_, ROOT.RooFit.Bins(28))
print "BRerr_frame set!"
BRpull_frame = mcstudy.plotPull(B_R_, ROOT.RooFit.Bins(28), ROOT.RooFit.FitGauss(1))
print "fit done"

#print "RMS=", BRpull_frame.getHist().GetRMS()
#print "mean=", BRpull_frame.getHist().GetMean()

#Plot distribution of minimized likelihood
NLLframe = mcstudy.plotNLL(ROOT.RooFit.Bins(28))

BRval_frame.SetTitle("")
BRval_frame.SetTitleOffset(1.5,"y")
BRval_frame.SetXTitle("BR(H#rightarrow#phi#gamma)")
BRerr_frame.SetTitle("")
BRerr_frame.SetTitleOffset(1.5,"y")
BRerr_frame.SetXTitle("#sigma_{BR(H#rightarrow#phi#gamma)}")
BRpull_frame.SetTitle("")
#BRpull_frame.SetMaximum(900)
BRpull_frame.SetTitleOffset(1.5,"y")
BRpull_frame.SetXTitle("PULL_{BR(H#rightarrow#phi#gamma)}")
NLLframe.SetTitle("")

#Actually plot
c1 = ROOT.TCanvas()
c1.cd()
BRval_frame.Draw()
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRval.pdf")
c1.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRval.png")

c2 = ROOT.TCanvas()
c2.cd()
BRerr_frame.Draw()
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRerr.pdf")
c2.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRerr.png")

c3 = ROOT.TCanvas()
c3.cd()
BRpull_frame.Draw()
c3.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull_cheby.pdf")
c3.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/BRpull_cheby.png")

c4 = ROOT.TCanvas()
c4.cd()
NLLframe.Draw()
c4.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/NLL.pdf")
c4.SaveAs("/eos/user/g/gumoret/www/latest_production/massanalysis_latest_production/NLL.png")
