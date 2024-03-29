import ROOT
import os, sys
import subprocess
import argparse

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to download MC or data')
p.add_argument('isData_option', help='Type <<MC>> or <<data>>')
args = p.parse_args()

if args.isData_option == "MC":
    isData = False
if args.isData_option == "data":
    isData = True

if not os.path.exists("rootfiles"):
    os.makedirs("rootfiles")
if not os.path.exists("rootfiles/latest_production"):
    os.makedirs("rootfiles/latest_production")
if not os.path.exists("rootfiles/latest_production/MC"):
    os.makedirs("rootfiles/latest_production/MC")
if not os.path.exists("rootfiles/latest_production/MC/signals"):
    os.makedirs("rootfiles/latest_production/MC/signals")
if not os.path.exists("rootfiles/latest_production/MC/backgrounds"):
    os.makedirs("rootfiles/latest_production/MC/backgrounds")
if not os.path.exists("rootfiles/latest_production/dataprocess"):
    os.makedirs("rootfiles/latest_production/dataprocess")

print "Processing ", args.isData_option

if not isData :
    dir_input = "crab_projects/samples_MC/"
    dir_output_bkg = "rootfiles/latest_production/MC/backgrounds/"
    dir_output_sig = "rootfiles/latest_production/MC/signals/"  
else :
    dir_input = "crab_projects/samples_data/"
    dir_output_data = "rootfiles/latest_production/dataprocess/"

list_dirs = os.listdir(dir_input)

for dirname in list_dirs:

    print "Processing sample dir " + dirname

    n_jobs_command = "crab status -d " + dir_input + dirname + " | grep status: " + "| awk " + """'{split($0,array,"/") ; print array[2]}'""" + "| sed 's/.$//'"
    n_jobs = int(subprocess.check_output(n_jobs_command, shell=True))

    print "Number of jobs to be retrieved: ", n_jobs

    if n_jobs <= 500:
        crab_command = "crab getoutput -d " + dir_input + dirname
        os.system(crab_command)
    elif (n_jobs > 500 and n_jobs <= 1000):
        crab_command = "crab getoutput -d " + dir_input + dirname + " --jobids 1-500"
        os.system(crab_command)
        crab_command_1 = "crab getoutput -d " + dir_input + dirname + " --jobids 501-" + str(n_jobs) # Because it is impossible to concatenate str and int objects
        os.system(crab_command_1)
    else:
        crab_command = "crab getoutput -d " + dir_input + dirname + " --jobids 1-500"
        os.system(crab_command)
        crab_command_1 = "crab getoutput -d " + dir_input + dirname + " --jobids 501-1000"
        os.system(crab_command_1)
        crab_command_2 = "crab getoutput -d " + dir_input + dirname + " --jobids 1001-" + str(n_jobs) # Because it is impossible to concatenate str and int objects
        os.system(crab_command_2)
    
    if args.isData_option == "MC":
        samplename = dirname.split("crab_HPhiGammaAnalysis_")
    if args.isData_option == "data":
        samplename = dirname.split("crab_2018_HPhiGammaAnalysis_")
    if "Signal" in dirname:
        hadd_command = "hadd -f " + dir_output_sig + "HPhiGammaAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"
    elif isData:
        hadd_command = "hadd -f " + dir_output_data + "HPhiGammaAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"
    else:
        hadd_command = "hadd -f " + dir_output_bkg + "HPhiGammaAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"

    os.system(hadd_command)

print "All done!"
