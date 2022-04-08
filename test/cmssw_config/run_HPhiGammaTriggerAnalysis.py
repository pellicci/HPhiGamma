import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("USER")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff') 
process.load('Geometry.CaloEventSetup.CaloTowerConstituents_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') #Could be not the coprrect one, but should contain the one without "condDBv2"
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
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
   inputFiles = {"root://cms-xrd-global.cern.ch//store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/100000/B29BDCC8-5429-7244-96D9-A5BF6C04BB73.root","/store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/110000/00A4B13C-EBCD-6C46-93D7-FA721EA4DB41.root","/store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/110000/012BD380-53E6-4B4E-802B-3F96BF7DA0E5.root","/store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/110000/027C0A31-BD7E-F04A-99C0-067C4F606BD1.root"}

else:
   process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v18')
   inputFiles = {"root://cms-xrd-global.cern.ch//store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/100000/B29BDCC8-5429-7244-96D9-A5BF6C04BB73.root","/store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/110000/00A4B13C-EBCD-6C46-93D7-FA721EA4DB41.root","/store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/110000/012BD380-53E6-4B4E-802B-3F96BF7DA0E5.root","/store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/110000/027C0A31-BD7E-F04A-99C0-067C4F606BD1.root"}
   
process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles),
                             duplicateCheckMode = cms.untracked.string ('noDuplicateCheck')
)


# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("HPhiGammaTriggerAnalysis_output.root")
)

process.load("HiggsAnalysis.HPhiGamma.HPhiGammaTriggerAnalysis_cfi")
process.HPhiGammaTriggerAnalysis.runningOnData = options.runningOnData

# Apply JEC to both MC and data
#process.HPhiGammaTriggerAnalysis.slimmedJets = cms.InputTag("updatedPatJetsUpdatedJEC")

###############################################
#                                             #
#------------------ Trigger ------------------#
#                                             #
###############################################


#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

# Trigger filters for muons
#process.trigger_filter = hlt.triggerResultsFilter.clone()
#process.trigger_filter.triggerConditions = cms.vstring('HLT_Photon35_TwoProngs35_v*')
#process.trigger_filter.hltResults = cms.InputTag("TriggerResults", "", "HLT")
#process.trigger_filter.l1tResults = cms.InputTag("")
#process.trigger_filter.throw = cms.bool( False )


###############################################
#                                             #
#----------------- Sequence ------------------#
#                                             #
###############################################

process.seq = cms.Path(process.HPhiGammaTriggerAnalysis)

process.schedule = cms.Schedule(process.seq)
