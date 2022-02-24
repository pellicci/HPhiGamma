import ROOT
#ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#Get the model and the data
fInputData = ROOT.TFile("histos/latest_production/histos_SR_Data.root")
fInputData.cd()
mytree = fInputData.Get("tree_output")

fInputPDF = ROOT.TFile("workspaces/ws_data.root")
workspace = fInputPDF.Get("myworkspace")
workspace.Print()

#Define the parameter of interest
B_R_ = workspace.var("B_R_")
poi = ROOT.RooArgSet(B_R_)

#Define observables
mass_KKg = workspace.var("mass_KKg")

observables = ROOT.RooArgSet()
observables.add(mass_KKg)

model = ROOT.RooStats.ModelConfig(workspace)
model.SetObservables(observables)
model.SetPdf("totPDF")

#Null Hypo
nullParams = poi.snapshot()
nullParams.setRealValue("B_R_",0.)

#Build the profile likelihood calculator: performs the scan of the likelihood ratio
plc = ROOT.RooStats.ProfileLikelihoodCalculator()

dataset = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass_KKg),ROOT.RooFit.Import(mytree))
plc.SetData(dataset)
plc.SetModel(model)
plc.SetParameters(poi)
plc.SetNullParameters(nullParams)

#The the Hypotest result
htr = plc.GetHypoTest()

print "-------------------------------------------------"
print "The p-value for the null is ", htr.NullPValue()
print "Corresponding to a signifcance of ", htr.Significance()
print "-------------------------------------------------"

#PyROOT sometimes fails cleaning memory, this helps
del plc

