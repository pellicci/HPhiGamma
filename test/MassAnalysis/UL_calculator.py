import ROOT

#Some useful things before starting --------------------------------------------------------------------
ROOT.gROOT.SetBatch(True) #Supress the opening of many Canvas's
ROOT.gROOT.ProcessLineSync(".L MassAnalysis/dCB/RooDoubleCBFast.cc+") #import the Doube Crystal Ball PDF

#Import from workspaces --------------------------------------------------------------------
fInput = ROOT.TFile("workspaces/ws_data_blind.root")
fInput.cd()
ws = fInput.Get("myworkspace")
ws.Print() #print some info of the workspace content

B_R_ = ws.var("B_R_")
poi_list = ROOT.RooArgSet(B_R_)
mass_KKg = ws.var("mesonGammaMass")
obs_list = ROOT.RooArgSet(ws.var("mesonGammaMass"))

#Retrieve the generated dataset in fit_Data_Blind.py
data = ws.data("totPDF_chebychevData")

#Set the RooModelConfig and let it know what the content of the workspace is about
model = ROOT.RooStats.ModelConfig()
model.SetWorkspace(ws)
model.SetPdf("totPDF_chebychev")
model.SetParametersOfInterest(poi_list)
model.SetObservables(obs_list)
model.SetName("S+B Model")
model.SetProtoData(data)
print "model set!"

bModel = model.Clone()
bModel.SetName("B model")
oldval = poi_list.find("B_R_").getVal()
poi_list.find("B_R_").setVal(0) #BEWARE that the range of the POI has to contain zero! -> signal is suppressed
bModel.SetSnapshot(poi_list)
poi_list.find("B_R_").setVal(oldval)

fc = ROOT.RooStats.AsymptoticCalculator(data, bModel, model)
fc.SetOneSided(1)
#Create hypotest inverter passing desired calculator
calc = ROOT.RooStats.HypoTestInverter(fc)
calc.SetConfidenceLevel(0.95)

#Use CLs
calc.UseCLs(1)
calc.SetVerbose(0)
#Configure ToyMC sampler
toymc = fc.GetTestStatSampler()
#Set profile likelihood test statistics
profl = ROOT.RooStats.ProfileLikelihoodTestStat(model.GetPdf())
#For CLs (bounded intervals) use one-sided profile likelihood
profl.SetOneSided(1)
#Set the test statistic to use
toymc.SetTestStatistic(profl)

npoints = 100 #Number of points to scan
# min and max for the scan (better to choose smaller intervals)
poimin = poi_list.find("B_R_").getMin()
poimax = poi_list.find("B_R_").getMax()

min_scan = 10**(-10)
max_scan = 10**(-2)
print "Doing a fixed scan  in interval : ",min_scan, " , ", max_scan
calc.SetFixedScan(npoints,min_scan,max_scan)
#calc.SetAutoScan()

# In order to use PROOF, one should install the test statistic on the workers
# pc = ROOT.RooStats.ProofConfig(workspace, 0, "workers=6",0)
# toymc.SetProofConfig(pc)

result = calc.GetInterval() #This is a HypoTestInveter class object

upperLimit = result.UpperLimit()

print "################"
print "The observed CLs upper limit is: ", upperLimit

##################################################################

#Compute expected limit
print "Expected upper limits, using the B (alternate) model : "
print " expected limit (median) ", result.GetExpectedUpperLimit(0)
print " expected limit (-1 sig) ", result.GetExpectedUpperLimit(-1)
print " expected limit (+1 sig) ", result.GetExpectedUpperLimit(1)
print "################"

del fc
