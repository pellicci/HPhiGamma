import FWCore.ParameterSet.Config as cms
process = cms.Process("USER")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff') 
process.load('Geometry.CaloEventSetup.CaloTowerConstituents_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') #Could be not the coprrect one, but should contain the one without "condDBv2"
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.parseArguments()

################################################################################################################
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,
#                       runEnergyCorrections=False, #as energy corrections are not yet availible for 2018
#                       era='2018-Prompt')  
################################################################################################################

#Input source
if options.runningOnData:
   process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018ABC_v2')
   inputFiles = {"root://cms-xrd-global.cern.ch//store/data/Run2018C/Tau/MINIAOD/17Sep2018-v1/270000/E95D5792-39AA-CE4B-B2BA-9A65FA639EAA.root"}

else:
   process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v18')
   #inputFiles = {"/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/F9947F2D-F185-4E43-9A4B-EA7FAF2CE4C2.root"}
   inputFiles = {"/store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_MINIAOD_10213V1/190507_131615/0004/HPhiGamma_signal_MINIAOD_4769.root"}

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles)
)

# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("HPhiGammaAnalysis_output.root")
)


###############################################################################################################################
#                                                                                                                             #
#-------------------------------------------------------- Jet corrections ----------------------------------------------------#
#                                                                                                                             #
#                                      https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC                                     #
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?redirectedfrom=CMS.WorkBookJetEnergyCorrections #
#                                                                                                                             #
###############################################################################################################################

#Use the latest JECs
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
if not options.runningOnData:      #This loop is for MC
   print "Jet Energy Corrections on Monte Carlo will be applied"
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')# The right corrections are taken through the chosen global tag
else:
   print "Jet Energy Corrections on Data will be applied "
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')# The right corrections are taken through the chosen global tag

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = jetCorrectionsList
)

process.load("HiggsAnalysis.HPhiGamma.HPhiGammaAnalysis_cfi")
process.HPhiGammaAnalysis.runningOnData = options.runningOnData

# Apply JEC to both MC and data
process.HPhiGammaAnalysis.slimmedJets = cms.InputTag("updatedPatJetsUpdatedJEC")

###############################################
#                                             #
#------------------ Trigger ------------------#
#                                             #
###############################################


#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

# Trigger filters for muons
process.trigger_filter = hlt.triggerResultsFilter.clone()
process.trigger_filter.triggerConditions = cms.vstring('HLT_Photon35_TwoProngs35_v*')
process.trigger_filter.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter.l1tResults = cms.InputTag("")
process.trigger_filter.throw = cms.bool( False )


###############################################
#                                             #
#----------------- Sequence ------------------#
#                                             #
###############################################

process.seq = cms.Path(process.trigger_filter * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.HPhiGammaAnalysis)

process.schedule = cms.Schedule(process.seq)
