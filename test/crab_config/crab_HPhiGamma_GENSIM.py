from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HPhiGamma_Pythia8_GENSIM_1026'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/GluGluHPhiGamma_Signal_GENSIM.py'
config.JobType.pluginName = 'PrivateMC'
config.JobType.outputFiles = ['HPhiGamma_signal_GENSIM.root']
config.JobType.numCores = 8

config.section_('Data')
config.Data.outputPrimaryDataset = 'HPhiGamma_GENSIM_1026'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 6
NJOBS = 8000 #Do not increase: maximum number of jobs per task is 10k
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'HPhiGamma_GENSIM_1026'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
