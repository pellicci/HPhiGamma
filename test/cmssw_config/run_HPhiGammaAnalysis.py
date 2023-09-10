import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("USER")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff') 
process.load('Geometry.CaloEventSetup.CaloTowerConstituents_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') #Could be not the coprrect one, but should contain the one without "condDBv2"
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.parseArguments()

################################################################################################################
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,
#                       runEnergyCorrections=False, #as energy corrections are not yet availible for 2018
#                       era='2018-Prompt')  
################################################################################################################

#INPUT FILE LIST      

'''                                                                                                                                                                                                   
    For the given path, get the List of all files in the directory tree                                                                                                                               
'''                                                                                                                                                                                                   
def getListOfFiles(dirName):                                                                                                                                                                          
    # create a list of file and sub directories                                                                                                                                                       
    # names in the given directory                                                                                                                                                                    
    listOfFile = os.listdir(dirName)                                                                                                                                                                  
    allFiles = list()                                                                                                                                                                                 
    # Iterate over all the entries                                                                                                                                                                    
    for entry in listOfFile:                                                                                                                                                                          
        # Create full path                                                                                                                                                                            
        fullPath = os.path.join(dirName, entry)                                                                                                                                                       
        # If entry is a directory then get the list of files in this directory                                                                                                                        
        if os.path.isdir(fullPath):                                                                                                                                                                   
            allFiles = allFiles + getListOfFiles(fullPath)                                                                                                                                            
        else:                                                                                                                                                                                         
            allFiles.append("file:" + fullPath)                                                                                                                                                                 
                                                                                                                                                                                                      
    return allFiles     

# Get the list of all files in directory tree at given path                                                                                                                                           
#listOfFiles = getListOfFiles(input_path)                                                                                                                                                              
#print(listOfFiles)


#Input source
if options.runningOnData:
    process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v33') # OLD ONE : 102X_dataRun2_Sep2018ABC_v2
    inputFiles = { '/store/data/Run2018B/Tau/MINIAOD/UL2018_MiniAODv2-v2/70000/09FD3540-D36B-6249-9817-D3BAC2F02E74.root'}
    #inputFiles = listOfFiles

else:
   process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v15_L1v1')  # OLD ONE : 102X_upgrade2018_realistic_v18
   #inputFiles = {"/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/F9947F2D-F185-4E43-9A4B-EA7FAF2CE4C2.root"}
   #inputFiles = listOfFiles
   #Choose this for Phi ggH signal
   inputFiles = {'/store/mc/RunIISummer20UL18MiniAODv2/GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/40000/2EBC4FBF-89AB-5946-A9A5-45A005D5EB0C.root','/store/mc/RunIISummer20UL18MiniAODv2/GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/40000/466C2A54-EAA8-D441-AECE-92FBA662611A.root','/store/mc/RunIISummer20UL18MiniAODv2/GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/40000/10358A5B-CBF1-554A-983D-37E0015C2DCC.root'} #'/store/mc/RunIISummer20UL18MiniAODv2/GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/40000/436E2054-BAEE-AA45-B654-5262771C26C7.root'
   #Choose this for Phi VBF signal   
   #inputFiles = {'/store/mc/RunIISummer20UL18MiniAODv2/VBF_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2560000/00B7EB0B-2697-F241-BE5F-1735C9656C30.root'}
   #Choose this for Rho VBF signal
   #inputFiles = {'/store/mc/RunIISummer20UL18MiniAODv2/GluGlu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/40000/0A1840F8-EB89-CC46-B448-A4DF868F13CF.root','/store/mc/RunIISummer20UL18MiniAODv2/VBF_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/40000/43213ED1-8A32-9A41-B818-382412E0B6AD.root','/store/mc/RunIISummer20UL18MiniAODv2/VBF_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/40000/1E9AFC8F-9D47-3F49-808D-33F4FF9086DE.root'}
   #Choose fot Rho ggH signal
   #inputFiles = {'/store/mc/RunIISummer20UL18MiniAODv2/GluGlu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/40000/1D0FF4E1-E390-2E40-B24A-3B196B7E8DA0.root'}
   #Choose for K0star ggH signal
   #inputFiles = {'/store/mc/RunIISummer20UL18MiniAODv2/GluGluHtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2820000/019F13F4-1C5C-4649-B4BF-0124156D43A1.root'}

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles),
                             duplicateCheckMode = cms.untracked.string ('noDuplicateCheck')
)


# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("HPhiGammaAnalysis_output.root")
)


###############################################################################################################################
#                                                                                                                               #
#-------------------------------------------------------- Jet corrections ----------------------------------------------------#
#                                                                                                                             #
#                                      https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC                                     #
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?redirectedfrom=CMS.WorkBookJetEnergyCorrections #
#                                                                                                                             #
###############################################################################################################################

#Use the latest JECs
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
if not options.runningOnData:      #This loop is for MC
   print "Jet Energy Corrections on Monte Carlo will be applied"
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')# The right corrections are taken through the chosen global tag
else:
   print "Jet Energy Corrections on Data will be applied "
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')# The right corrections are taken through the chosen global tag

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = jetCorrectionsList
)

#EGamma correction
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2018-UL')    
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)


process.load("HiggsAnalysis.HPhiGamma.HPhiGammaAnalysis_cfi")
process.HPhiGammaAnalysis.runningOnData = options.runningOnData

# Apply JEC to both MC and data
process.HPhiGammaAnalysis.slimmedJets = cms.InputTag("updatedPatJetsUpdatedJEC")

###############################################
#                                             #
#------------------ Trigger ------------------#
#                                             #
###############################################


#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

# Trigger filters for muons
#process.trigger_filter = hlt.triggerResultsFilter.clone()
#process.trigger_filter.triggerConditions = cms.vstring('HLT_Photon35_TwoProngs35_v*')
#process.trigger_filter.hltResults = cms.InputTag("TriggerResults", "", "HLT")
#process.trigger_filter.l1tResults = cms.InputTag("")
#process.trigger_filter.throw = cms.bool( False )


###############################################
#                                             #
#----------------- Sequence ------------------#
#                                             #
###############################################

process.seq = cms.Path(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.egammaPostRecoSeq * process.HPhiGammaAnalysis)

process.schedule = cms.Schedule(process.seq)
