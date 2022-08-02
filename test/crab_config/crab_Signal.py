from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HPhiGammaAnalysis_Signal'
config.General.workArea = 'crab_projects/samples_signal'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_HPhiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MCpileUp_2018_25ns_UltraLegacy_PoissonOOTPU.root','MyDataPileupHistogram.root' ] #data files for PileUp reweighting
config.JobType.outputFiles = ['HPhiGammaAnalysis_output.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.pyCfgParams = ['runningOnData=False'] # Configure 2016 data jobs - muons


config.section_('Data')
config.Data.inputDataset = '/HPhiGamma_GENSIM_1026/pellicci-HPhiGamma_MINIAOD_10213V1-b4c7160d7e500816c9f483e0e4630d18/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 15
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
