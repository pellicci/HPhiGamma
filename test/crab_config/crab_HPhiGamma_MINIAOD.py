from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
 
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HPhiGamma_Pythia8_MINIAOD_10213V2'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/GluGluHPhiGamma_Signal_MINIAOD.py'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
config.Data.inputDataset = '/HPhiGamma_GENSIM_1026/pellicci-HPhiGamma_RECO_10213V1-889dc1528a7a51924818b24d55ffac27/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

config.Data.outputDatasetTag = 'HPhiGamma_MINIAOD_10213V2'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
