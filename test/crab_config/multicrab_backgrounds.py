from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects/samples_MC'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_HPhiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MCpileUp_2018_25ns_UltraLegacy_PoissonOOTPU.root','MyDataPileupHistogram.root' ] #data files for PileUp reweighting
config.JobType.outputFiles = ['HPhiGammaAnalysis_output.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.pyCfgParams = ['runningOnData=False'] # Configure 2016 data jobs - muons

config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    #################################################
    #                                               #
    #--------------- Running 2018 MC ---------------#
    #                                               #
    #################################################
"""
config.General.requestName = 'HPhiGammaTwoProngsTriggerAnalysis_DY50' #Requires FileBased splitting
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()


config.General.requestName = 'HPhiGammaTwoProngsTriggerAnalysis_DY10to50'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()
"""
#config.General.requestName = 'HPhiGammaAnalysis_ttbarToHadronic'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM' #OLD
#config.Data.inputDataset = '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' #NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_ttbarToSemiLeptonic' #Requires FileBased splitting
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext3-v2/MINIAODSIM' #OLD
#config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_ttbarlnu'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM' #OLD
#config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_DY10to50'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM' #OLD
#config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' #NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_DY50'#Requires FileBased job splitting
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM'#OLD
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_QCDpT30to50'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_QCDpT50to80'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_QCDpT80to120'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/QCD_Pt-80to120_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_QCDpT120to170'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/QCD_Pt-120to170_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_QCDpT170to300'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_QCDpT300toInf'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/QCD_Pt-300toInf_EMEnriched_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()


#config.General.requestName = 'HPhiGammaAnalysis_WW'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'#OLD
#config.Data.inputDataset = '/WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT100to200'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' #NEW
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT200to400'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_DR-0p4_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' #NEW
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT400to600'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_DR-0p4_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' #NEW
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT600toInf'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_DR-0p4_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' #NEW
p = Process(target=submit, args=(config,))
p.start()
p.join()


#config.General.requestName = 'HPhiGammaAnalysis_WZ'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM'#OLD
#config.Data.inputDataset = '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' #NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_WJetsToLNu0J'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'#OLD
#config.Data.inputDataset = '/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_WJetsToLNu1J'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'#OLD
#config.Data.inputDataset = '/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()


#config.General.requestName = 'HPhiGammaAnalysis_WJetsToLNu2J'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'#OLD
#config.Data.inputDataset = '/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_DiPhotonJets'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/DiPhotonJets_MGG-80toInf_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'#OLD
#config.Data.inputDataset = '/DiPhotonJets_MGG-80toInf_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()

#config.General.requestName = 'HPhiGammaAnalysis_ZGammaToLLGamma'
#config.Data.unitsPerJob = 5
#config.Data.inputDataset = '/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'#NEW
#p = Process(target=submit, args=(config,))
#p.start()
#p.join()
