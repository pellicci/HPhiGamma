from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HPhiGammaAnalysis_Signal'
config.General.workArea = 'crab_projects/samples'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_HPhiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MCpileUp_2018_25ns_JuneProjectionFull18_PoissonOOTPU.root','MyDataPileupHistogram.root' ] #data files for PileUp reweighting
config.JobType.outputFiles = ['HPhiGammaAnalysis_output.root']

config.section_('Data')
config.Data.inputDataset = '/HPhiGamma_GENSIM_1026/pellicci-HPhiGamma_MINIAOD_10213V1-b4c7160d7e500816c9f483e0e4630d18/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
