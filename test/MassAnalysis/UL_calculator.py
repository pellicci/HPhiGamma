import ROOT

fInputData = ROOT.TFile("histos/latest_production/histos_Data.root")
fInputData.cd()
mytree = fInputData.Get("tree_output")

fInputPDF = ROOT.TFile("workspaces/ws_data.root")
ws = fInputPDF.Get("myworkspace")
ws.Print()

B_R_ = ws.var("B_R_")
poi_list = ROOT.RooArgSet(B_R_)
mass_KKg = ws.var("mass_KKg")
obs_list = ROOT.RooArgSet(ws.var("mass_KKg"))
data = ROOT.RooDataSet("dataset","dataset",ROOT.RooArgSet(mass_KKg),ROOT.RooFit.Import(mytree))
#data_binned = data.binnedClone("data_binned","data_binned")
#getattr(ws,'import')(data_binned)

ws.Print()

#glb_list = ROOT.RooArgSet()
#glb_list.add(ws.var("global_eff"))

#Set the RooModelConfig and let it know what the content of the workspace is about
model = ROOT.RooStats.ModelConfig()
model.SetWorkspace(ws)
model.SetPdf("totPDF")
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
max_scan = 10**(-4)
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


'''
fc = ROOT.RooStats.FeldmanCousins(data_binned,model)
fc.UseAdaptiveSampling(1)
fc.FluctuateNumDataEntries(1)
fc.SetNBins(10)
fc.SetTestSize(.05)

result = fc.GetInterval() 

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



mcmc = ROOT.MCMCCalculator(data_binned,model);
mcmc.SetConfidenceLevel(0.95)
sp = ROOT.SequentialProposal(0.1)
mcmc.SetProposalFunction(sp)
mcmc.SetNumIters(1000000)  #Metropolis-Hastings algorithm iterations
mcmc.SetNumBurnInSteps(50) #first N steps to be ignored as burn-in
mcmc.SetLeftSideTailFraction(0.)

interval = mcmc.GetInterval()

fc = ROOT.RooStats.FrequentistCalculator(data_binned, bModel, model)
fc.SetToys(350,350)

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

npoints = 10 #Number of points to scan
# min and max for the scan (better to choose smaller intervals)
poimin = poi_list.find("br_emu").getMin()
poimax = poi_list.find("br_emu").getMax()

min_scan = 0.000000001
max_scan = 0.000001
print "Doing a fixed scan  in interval : ",min_scan, " , ", max_scan
#calc.SetFixedScan(npoints,min_scan,max_scan)
calc.SetAutoScan()

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

#Plot the results
freq_plot = ROOT.RooStats.HypoTestInverterPlot("HTI_Result_Plot","Frequentist scan result for the W -> pigamma BR",result)

#xPlot in a new canvas with style
canvas = ROOT.TCanvas()
canvas.cd()
#freq_plot.Draw("2CL")
freq_plot.Draw("EXP")
# freq_plot.GetYaxis().SetRangeUser(0.,0.8)
# freq_plot.GetXaxis().SetRange(0.,0.0000107)
canvas.SaveAs("plots/latest_production/2016_2017_2018/UL_CLs.pdf")
'''

del fc
