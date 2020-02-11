python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_DY10to50.root             
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_DY50.root             
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT100to200.root        
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDHT200to300.root                
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
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu0J.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu1J.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu2J.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WW.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WZ.root         
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarToSemiLeptonic.root  
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarToHadronic.root  
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarlnu.root

#Merge background samples
hadd histos/latest_production/all_background.root histos/latest_production/histos*.root

#Generate signal histos
python generate_histos.py rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_Signal.root

#Merge signal + background
hadd histos/latest_production/all_background_singal.root histos/latest_production/histos_*.root
