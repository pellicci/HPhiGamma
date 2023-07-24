import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
            '/store/mc/RunIISummer20UL18RECO/GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v3/40000/0867E20C-6C92-DD44-B9D2-0B1ED385105C.root'
                )
                            )

#Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("HPhiGammaGenLevel_output.root")
)

process.HPhiGammaGenLevel = cms.EDAnalyzer('HPhiGammaGenLevel'
                              )

process.p = cms.Path(process.HPhiGammaGenLevel)
