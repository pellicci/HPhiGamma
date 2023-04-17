import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
            'file:/eos/user/p/pellicci/MesonGamma_root/2022//eos/user/p/pellicci/MesonGamma_root/2022/ZPhiGamma_RAW/Zphigamma_RAW2022_10.root','file:/eos/user/p/pellicci/MesonGamma_root/2022//eos/user/p/pellicci/MesonGamma_root/2022/ZPhiGamma_RAW/Zphigamma_RAW2022_11.root',
            'file:/eos/user/p/pellicci/MesonGamma_root/2022//eos/user/p/pellicci/MesonGamma_root/2022/ZPhiGamma_RAW/Zphigamma_RAW2022_100.root','file:/eos/user/p/pellicci/MesonGamma_root/2022//eos/user/p/pellicci/MesonGamma_root/2022/ZPhiGamma_RAW/Zphigamma_RAW2022_101.root'
                )
                            )

#Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("ZPhiGammaGenLevel_output.root")
)

process.ZPhiGammaGenLevel = cms.EDAnalyzer('ZPhiGammaGenLevel'
                              )

process.p = cms.Path(process.ZPhiGammaGenLevel)
