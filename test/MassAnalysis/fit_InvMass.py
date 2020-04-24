import ROOT

#INPUT SETTINGS
fInput = ROOT.TFile("workspaces/ws_dataCombined.root")
fInput.cd()
workspace = fInput.Get("myworkspace")
workspace.Print()

#VAR SETTINGS
mass_KKg = workspace.var("mass_KKg")
a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",0.1,-10.,10.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",0.,-10.,10.)
c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.,-10.,10.)
d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",0.,-10.,10.)
e_bkg = ROOT.RooRealVar("e_bkg","e_bkg",0.,-10.,10.)

#PDF
bkgPDF = ROOT.RooChebychev("bkgPDF","bkgPDF",mass_KKg,ROOT.RooArgList(a_bkg,b_bkg,c_bkg,d_bkg,e_bkg))
print "bkgPDF created!"

#DATASET
dataCombined = workspace.data("totPDFData")
print "dataCombined retrieved!"

#FIT
bkgPDF.fitTo(dataCombined)
print "dataCombined fitted!"

xframe = mass_KKg.frame(40)
dataCombined.plotOn(xframe)
bkgPDF.plotOn(xframe)

c1 = ROOT.TCanvas()
c1.cd()
xframe.Draw()
c1.SaveAs("fit_dataCombined.pdf")
