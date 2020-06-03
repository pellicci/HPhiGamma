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
config.JobType.inputFiles = ['MCpileUp_2018_25ns_JuneProjectionFull18_PoissonOOTPU.root','MyDataPileupHistogram.root' ] #data files for PileUp reweighting
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

config.General.requestName = 'HPhiGammaAnalysis_ttbarToHadronic'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM' #There is also a version without ext2 (and with v1), it has a few less events but still more than 100M. This one has 200M, I think it is enough
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_ttbarToSemiLeptonic' #Requires FileBased splitting
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext3-v2/MINIAODSIM' #Same story as TTToHadronic
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_ttbarlnu'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_DY10to50'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_DY50'#Requires FileBased job splitting
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_QCDHT100to200'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_QCDHT200to300'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join() 

config.General.requestName = 'HPhiGammaAnalysis_QCDHT300to500'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_QCDHT500to700'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_QCDHT700to1000'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_QCDHT1000to1500'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_QCDHT1500to2000'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_QCDHT2000toInf'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_WW'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_WZ'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM' 
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJets20to40'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJets20toInf'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJets40toInf'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()


config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT600toInf'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT40to100'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT400to600'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT100to200'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-4cores5k_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_GammaJetsHT200to400'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_WJetsToLNu0J'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_WJetsToLNu1J'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()


config.General.requestName = 'HPhiGammaAnalysis_WJetsToLNu2J'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'HPhiGammaAnalysis_DiPhotonJets'
config.Data.unitsPerJob = 5
config.Data.inputDataset = '/DiPhotonJets_MGG-80toInf_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()
