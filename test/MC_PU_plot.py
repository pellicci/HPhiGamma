from ROOT import *
import numpy as np

MC_PU_scenario_2018 = [4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05,
                       3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473,
                       0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138,
                       0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411,
                       0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554,
                       0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895,
                       0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877,
                       0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612,
                       0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551,
                       0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934,
                       0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915,
                       0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932,
                       0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885,
                       0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012,
                       0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05,
                       2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06,
                       3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07,
                       5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07,
                       1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08,
                       6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08]

print len(MC_PU_scenario_2018)

h_PU = TH1F("pileup","",99,0.,99.)

i = 0

for pu in MC_PU_scenario_2018:
    h_PU.Fill(i,pu)
    i+=1

print h_PU.Integral()
canvas = TCanvas()
h_PU.Draw("hist")
out_file = TFile("MCpileUp_2018_25ns_JuneProjectionFull18_PoissonOOTPU.root","RECREATE")
h_PU.Write()

raw_input()
