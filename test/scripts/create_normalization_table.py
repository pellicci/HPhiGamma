###All normalizations are provided to 1fb-1 of lumi in these tables
import os
import sys
import argparse

#---------------------------------#

#xsec lookup in pb
secs_table = dict()

#secs_table["QCDPt30toInf"]        = 242700.0
#secs_table["QCDPt40toInf"]        = 117400.0
#secs_table["SingleToptW"]         = 34.91
#secs_table["SingleAntiToptW"]     = 34.97
#secs_table["Wlnu"]                = 52850.0

#Signal
secs_table["Signal_Phi_ggH"]       = 48.58*10**(-5)*0.49 #cross section * B.R.(10^-5) * BR(Phi->KK)
secs_table["Signal_Rho_ggH"]       = 48.58*10**(-5) #cross section * B.R.(10^-5)
secs_table["Signal_Phi_VBF"]       =  3.78*10**(-5)*0.49
secs_table["Signal_Rho_VBF"]       =  3.78*10**(-5)
'''
#ttbar
#secs_table["ttbarToSemiLeptonic"] = 687.1 #OLD = 365.34 # accounting for the 2 possible charge signs of the W
#secs_table["ttbarToHadronic"]     = 687.1 #OLD = 377.96
#secs_table["ttbarlnu"]            = 687.1 #OLD = 88.29 #NNLO-2018

#Drell-Yan
secs_table["DY10to50"]            = 15890.0 #OLD = 18810.0
secs_table["DY50"]                =  6404.0 #OLD = 2075.14*3 #amcatnlo 2017

#Bosons
secs_table["WW"]                  = 0. #No xsec on DAS #OLD = 12.178
secs_table["WZ"]                  = 27.59 #OLD = 27.6

#GammaJets
secs_table["GammaJetsHT40to100"]  =    18700.0  * 1.26 #UL   
secs_table["GammaJetsHT100to200"] =     5034.0  * 1.26 #UL
secs_table["GammaJetsHT200to400"] =     1128.0  * 1.26 #UL
secs_table["GammaJetsHT400to600"] =      124.8  * 1.26 #UL to check
secs_table["GammaJetsHT600toInf"] =       40.72 * 1.26 #UL

#QCD
secs_table["QCDpT15to20"]         =   66590.0 * 1.26 #UL 
secs_table["QCDpT20to30"]         = 1324000.0 * 1.26 #UL 
secs_table["QCDpT30to50"]         = 6447000.0 * 1.26 #UL 
secs_table["QCDpT50to80"]         = 1988000.0 * 1.26 #UL  
secs_table["QCDpT80to120"]        =  367500.0 * 1.26 #UL
secs_table["QCDpT120to170"]       =   66590.0 * 1.26 #UL
secs_table["QCDpT170to300"]       =   16620.0 * 1.26 #UL 
secs_table["QCDpT300toInf"]       =    1104.0 * 1.26 #UL

#W jets
secs_table["WJetsToLNu0J"]        = 53330.0 #OLD = 50131.98
secs_table["WJetsToLNu1J"]        =  8875.0 #OLD = 8426.09
secs_table["WJetsToLNu2J"]        =  3338.0 #OLD = 3172.96

#Di-photon
secs_table["DiPhotonJets"]        = 126.2 #OLD = 134.3

#ZGamma
secs_table["ZGammaToLLGamma"]     =     55.48 #UL
'''
#fraction of negative-weighted events in NLO samples (2018)
frac_table = dict()

#frac_table["QCDPt30toInf"]        = 0.
#frac_table["QCDPt40toInf"]        = 0. 
#frac_table["SingleToptW"]         = 0.003758
#frac_table["SingleAntiToptW"]     = 0.0034
#frac_table["Wlnu"]                = 0.0003866

#Signal
frac_table["Signal_Phi_ggH"]       = 0.
frac_table["Signal_Rho_ggH"]       = 0.
frac_table["Signal_Phi_VBF"]       = 0.
frac_table["Signal_Rho_VBF"]       = 0.
'''
#ttbar
frac_table["ttbarToSemiLeptonic"] = 0.003957 #OLD = 0.
frac_table["ttbarToHadronic"]     = 0.00388 #OLD = 0.
frac_table["ttbarlnu"]            = 0.00379 #OLD = 0.

#Drell-Yan
frac_table["DY50"]                = 0.1643 #UL
frac_table["DY10to50"]            = 0.     #OLD = 0.1367

#Bosons
frac_table["WW"]                  = 0. #No negative fraction on DAS #OLD = 0.001755
frac_table["WZ"]                  = 0. #UL

#QCD
frac_table["QCDpT15to20"]         = 0. #UL
frac_table["QCDpT20to30"]         = 0. #UL
frac_table["QCDpT30to50"]         = 0. #UL
frac_table["QCDpT50to80"]         = 0. #UL
frac_table["QCDpT80to120"]        = 0. #UL
frac_table["QCDpT120to170"]       = 0. #UL
frac_table["QCDpT170to300"]       = 0. #UL
frac_table["QCDpT300toInf"]       = 0. #UL

#GammaJets
frac_table["GammaJetsHT40to100"]  = 9.997e-06 #UL
frac_table["GammaJetsHT100to200"] = 6.985e-05 #UL
frac_table["GammaJetsHT200to400"] = 0.0003769 #UL
frac_table["GammaJetsHT400to600"] = 0.0008598 #UL to check
frac_table["GammaJetsHT600toInf"] = 0.002077 #UL

#W jets
frac_table["WJetsToLNu0J"]        = 0.09824 #OLD = 0.09868
frac_table["WJetsToLNu1J"]        = 0.2584  #OLD = 0.269
frac_table["WJetsToLNu2J"]        = 0.3432  #OLD = 0.3460

#Di-photon
frac_table["DiPhotonJets"]        = 0.2363 #OLD = 0.2238

#ZGamma
frac_table["ZGammaToLLGamma"]     = 0.1847 #UL
'''



##Now starts the program
def main():

    dir_input = "crab_projects/samples_MC/"
    list_dirs = os.listdir(dir_input)

    if not os.path.exists("rootfiles"):
        os.makedirs("rootfiles")
    if not os.path.exists("rootfiles/latest_production"):
        os.makedirs("rootfiles/latest_production")
    if not os.path.exists("rootfiles/latest_production/MC"):
        os.makedirs("rootfiles/latest_production/MC")
    if not os.path.exists("rootfiles/latest_production/MC/normalizations"):
        os.makedirs("rootfiles/latest_production/MC/normalizations")

    output_filename = "rootfiles/latest_production/MC/normalizations/Normalizations_table.txt" 
    out_file = open(output_filename,"w")

    events_cumul = dict()

    for dirname in list_dirs:

        print "dirname = ",dirname
        samplename = dirname.split("crab_HPhiGammaAnalysis_")[1]
        print "samplename = ",samplename

        print "Processing sample dir " + dirname
        crab_command = "crab report -d " + dir_input + dirname + " | grep read | awk '{print $5}'"
        print crab_command
        
        event_string = os.popen(crab_command).read()
        number_events = int(event_string)
        print "No. of events processed = " + event_string + "\n"
        
        events_cumul[samplename] = number_events*(1-2*frac_table[samplename])
        
    for sample,event_count in events_cumul.iteritems():
        if event_count == 0:
            scale_factor = 0.
            print "NUMBER OF EVENTS RETRIEVED = 0. SCALE FACTOR SET TO 0"
        else:
            xsection = secs_table[sample]
            print "Cross section = ", xsection
            scale_factor = float(xsection*1000./events_cumul[sample])
            print sample + " scale factor = ", scale_factor
            
        write_string = sample + " " + str(scale_factor) + "\n"
        print "Output Norm = ", write_string
        out_file.write(write_string)        
        
    print "All done!"
    
if __name__ == "__main__":
    main()
