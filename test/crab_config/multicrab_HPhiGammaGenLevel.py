from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_HPhiGammaGenLevel.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['HPhiGammaGenLevel_output.root']
config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(9_4_10). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

config.General.workArea = 'crab_projects/samples_HPhiGammaGenLevel/'

config.section_('Data')
config.Data.inputDBS = 'phys03'
#config.Data.splitting = 'Automatic'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 6
NJOBS = 8000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFNDirBase = '/store/user/gumoret/' 
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


    #config.JobType.pyCfgParams = ['runningOnData=False','runningEra=2'] # Configure 2018 MC signal jobs 

    config.General.requestName = 'HPhiGammaGenLevel'
    config.Data.inputDataset = '/HPhiGamma_GENSIM_1026/pellicci-HPhiGamma_GENSIM_1026-3ce49215944f840647bbc4e19518f7a7/USER'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
        
