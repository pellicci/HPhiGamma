import FWCore.ParameterSet.Config as cms

HPhiGammaAnalysis = cms.EDAnalyzer('HPhiGammaAnalysis',
                                  packedPFCandidates = cms.InputTag("packedPFCandidates"),
                                  slimmedMuons       = cms.InputTag("slimmedMuons"),
                                  prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                  slimmedPhotons     = cms.InputTag("slimmedPhotons"),
                                  slimmedElectrons   = cms.InputTag("slimmedElectrons"),
                                  slimmedJets        = cms.InputTag("slimmedJets"),
                                  slimmedMETs        = cms.InputTag("slimmedMETs"),
                                  slimmedMETsPuppi   = cms.InputTag("slimmedMETsPuppi"),
                                  runningOnData      = cms.bool(False),
                                  pvCollection       = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                  bsCollection       = cms.InputTag("offlineBeamSpot"),
                                  PileupSrc          = cms.InputTag("slimmedAddPileupInfo"),
                                  triggerbits        = cms.InputTag("TriggerResults","","HLT"),
                                  rho                = cms.InputTag("fixedGridRhoFastjetAll"),
                                  phoIdVerbose = cms.bool(False),
                                  effAreasConfigFile_el = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"),
                                  effAreasConfigFile_ph = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased.txt")
)
