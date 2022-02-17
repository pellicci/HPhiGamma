import ROOT
import os
import argparse

#INPUT FILES
fIn_bkg  = ROOT.TFile("../histos/latest_production/histos_CR_Sidebands.root")
tree_bkg = fIn_bkg.Get("tree_output")
fIn_sig  = ROOT.TFile("../histos/latest_production/histos_SR_Signal.root")
tree_sig = fIn_sig.Get("tree_output")

#OUTPUT FILE
fOut = ROOT.TFile("outputs/Nominal_training.root","RECREATE")

#START MVA
ROOT.TMVA.Tools.Instance()

factory = ROOT.TMVA.Factory("TMVAClassification", fOut,":".join(["!V","Transformations=I;D;P;G,D","AnalysisType=Classification"]))

dataloader = ROOT.TMVA.DataLoader()

#VARIABLES FROM THE TREE -------------------------------------------------------------------------------

#first track 
#dataloader.AddVariable("_firstCandPt","F") # Both Float and Double variable types must be indicated as F
dataloader.AddVariable("_firstCandIso","F")
#dataloader.AddVariable("_firstCandEta","F")

#second track
#dataloader.AddVariable("_secondCandPt","F")
dataloader.AddVariable("_secondCandIso","F")
#dataloader.AddVariable("_secondCandEta","F")

#best pair
#dataloader.AddVariable("_bestCouplePt","F")
#dataloader.AddVariable("_coupleIso","F")  
#dataloader.AddVariable("_bestCoupleEta","F")
#dataloader.AddVariable("_bestCoupleDeltaR","F")
#dataloader.AddVariable("mass_KK","F") #using m_kk variable is not correct since you use it to define the bkg estimation, use it only to compute scatter plots to see correlation with other variables

#photon
dataloader.AddVariable("_photonEt","F")
#dataloader.AddVariable("_photonEta","F")

#jet
dataloader.AddVariable("_bestJetPt","F")

#-------------------------------------------------------------------------------------------------------

sig_weight = 1.
bkg_weight = 1.

print "before AddSignalTree"
dataloader.AddSignalTree(tree_sig, sig_weight)
print "before AddBackgroundTree"
dataloader.AddBackgroundTree(tree_bkg, bkg_weight)

dataloader.SetWeightExpression("_eventWeight") #_eventWeight is the weight variable of the tree


mycutSig = ROOT.TCut("")
mycutBkg = ROOT.TCut("")


dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"]))

#dataloader.PrepareTrainingAndTestTree(mycutBkg, ":".join(["!V","nTrain_Background=0:nTest_Background=0"]))

method_btd  = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()

#weightfile_dir = "default/weights/TMVAClassification_BDT.weights.xml"

#if evaluate_BDT_systematic:
    # weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_shifted.xml"
    # weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_shifted.xml"
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_secondPart.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_secondPart.xml"
#elif test_on_signal_shifted_up:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Pythia_up.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Pythia_up.xml"
#elif test_on_signal_shifted_down:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Pythia_down.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Pythia_down.xml"
#elif test_on_signal_sin2:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_sin2.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_sin2.xml"
#elif test_on_signal_cos:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_cos_with_flat_theta.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_cos_with_flat_theta.xml"
#else:
    #weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu.xml"
    #weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele.xml"
    #weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Wmass.xml"
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_firstPart.xml"
    #weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Wmass.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_firstPart.xml"

#if isMuon:
 #   rename_weightfile = "mv " + weightfile_dir + " " + weightfile_mu
  #  os.system(rename_weightfile)
#else:
 #   rename_weightfile = "mv " + weightfile_dir + " " + weightfile_ele
  #  os.system(rename_weightfile)