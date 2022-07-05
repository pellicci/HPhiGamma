from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects/samples_data'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_HPhiGammaTwoProngsTriggerAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['HPhiGammaTwoProngsTriggerAnalysis_output.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.pyCfgParams = ['runningOnData=True']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = False
config.Data.lumiMask = 'json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

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


    config.General.requestName = '2018_HPhiGammaTwoProngsTriggerAnalysis_SingleMuon'
    config.Data.inputDataset = '/SingleMuon/Run2018B-UL2018_MiniAODv2-v2/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

