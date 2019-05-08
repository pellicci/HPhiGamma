from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
 
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HPhiGamma_Pythia8_HLT_10213V1'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/GluGluHPhiGamma_Signal_HLT.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDataset = '/HPhiGamma_GENSIM_1026/pellicci-HPhiGamma_GENSIM_1026-3ce49215944f840647bbc4e19518f7a7/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

config.Data.outputDatasetTag = 'HPhiGamma_HLT_10213V1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
