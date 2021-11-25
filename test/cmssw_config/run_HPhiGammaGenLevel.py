import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
            'root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_10.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_11.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_12.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_13.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_14.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_15.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_16.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_17.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_18.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_19.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_20.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_21.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_22.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_23.root','root://cms-xrd-global.cern.ch//store/user/pellicci/HPhiGamma_GENSIM_1026/HPhiGamma_GENSIM_1026/190409_152135/0000/HPhiGamma_signal_GENSIM_3701.root'
                )
                            )

#Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("HPhiGammaGenLevel_output.root")
)

process.HPhiGammaGenLevel = cms.EDAnalyzer('HPhiGammaGenLevel'
                              )

process.p = cms.Path(process.HPhiGammaGenLevel)
