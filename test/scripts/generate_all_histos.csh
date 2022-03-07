
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_DY10to50.root             
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_DY50.root             
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT100to200.root        
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT200to300.root                
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT300to500.root                
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT500to700.root                
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT700to1000.root                
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT1000to1500.root                
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT1500to2000.root               
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT2000toInf.root       
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJets20to40.root      
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJets20toInf.root      
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJets40toInf.root      
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT40to100.root      
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT100to200.root      
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT200to400.root      
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT400to600.root      
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT600toInf.root      
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu0J.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu1J.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu2J.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WW.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_DiPhotonJets.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WZ.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarToSemiLeptonic.root  
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarToHadronic.root  
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarlnu.root

#Generate signal histos
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_Signal.root

#Generate data histos
python generate_histos.py rootfiles/latest_production/dataprocess/HPhiGammaAnalysis_Data.root

#Merge QCD histos
hadd -f histos/latest_production/histos_QCD.root histos/latest_production/histos_QCDHT*
rm histos/latest_production/histos_QCDHT*

#Merge GammaJetEnriched histos
hadd -f histos/latest_production/histos_GammaJetsEnriched.root histos/latest_production/histos_GammaJets??to*
rm histos/latest_production/histos_GammaJets??to*

#Merge GammaJet histos
hadd -f histos/latest_production/histos_GammaJets.root histos/latest_production/histos_GammaJetsHT*
rm histos/latest_production/histos_GammaJetsHT*

#Merge background samples
#hadd histos/latest_production/all_background.root histos/latest_production/histos*.root

#Merge signal + background
#hadd histos/latest_production/all_background_singal.root histos/latest_production/histos_*.root
