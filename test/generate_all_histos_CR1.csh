
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ZGammaToLLGamma.root histos/latest_production/histos_CR1_ZGammaToLLGamma.root            
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_DiPhotonJets.root histos/latest_production/histos_CR1_DiPhotonJets.root        
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu2J.root histos/latest_production/histos_CR1_WJetsToLNu2J.root     
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu1J.root histos/latest_production/histos_CR1_WJetsToLNu1J.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WJetsToLNu0J.root histos/latest_production/histos_CR1_WJetsToLNu0J.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WZ.root histos/latest_production/histos_CR1_WZ.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT600toInf.root histos/latest_production/histos_CR1_GammaJetsHT600toInf.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT400to600.root histos/latest_production/histos_CR1_GammaJetsHT400to600.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT200to400.root histos/latest_production/histos_CR1_GammaJetsHT200to400.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_GammaJetsHT100to200.root histos/latest_production/histos_CR1_GammaJetsHT100to200.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_WW.root histos/latest_production/histos_CR1_WW.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDpT300toInf.root histos/latest_production/histos_CR1_QCDpT300toInf.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDpT170to300.root histos/latest_production/histos_CR1_QCDpT170to300.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDpT120to170.root histos/latest_production/histos_CR1_QCDpT120to170.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDpT80to120.root histos/latest_production/histos_CR1_QCDpT80to120.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDpT50to80.root histos/latest_production/histos_CR1_QCDpT50to80.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_QCDpT30to50.root histos/latest_production/histos_CR1_QCDpT30to50.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_DY50.root histos/latest_production/histos_CR1_DY50.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_DY10to50.root histos/latest_production/histos_CR1_DY10to50.root 
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarlnu.root histos/latest_production/histos_CR1_ttbarlnu.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarToSemiLeptonic.root histos/latest_production/histos_CR1_ttbarToSemiLeptonic.root
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_ttbarToHadronic.root histos/latest_production/histos_CR1_ttbarToHadronic.root

#Generate signal histos
python generate_histos.py 1 0 rootfiles/latest_production/MC/backgrounds/HPhiGammaAnalysis_Signal.root histos/latest_production/histos_CR1_Signal.root

#Generate data histos
python generate_histos.py 1 0 rootfiles/latest_production/dataprocess/HPhiGammaAnalysis_Data.root histos/latest_production/histos_CR1_Data.root

#Merge QCD histos
hadd -f histos/latest_production/histos_CR1_QCD.root histos/latest_production/histos_CR1_QCDpT*
rm histos/latest_production/histos_CR1_QCDpT*

#Merge GammaJet histos
hadd -f histos/latest_production/histos_CR1_GammaJets.root histos/latest_production/histos_CR1_GammaJetsHT*
rm histos/latest_production/histos_CR1_GammaJetsHT*

#Merge background samples
#hadd histos/latest_production/CR1_all_background.root histos/latest_production/histos_CR1*.root

#Merge signal + background
#hadd histos/latest_production/CR1_all_background_signal.root histos/latest_production/histos_CR1*.root
