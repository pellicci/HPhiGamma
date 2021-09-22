#!/bin/bash

HOMEDIR=/afs/cern.ch/user/p/pellicci/work/HPhiGamma/Prod/CMSSW_10_2_16_UL/src/HiggsAnalysis/HPhiGamma/generation
CMSSW_TO_USE=CMSSW_10_2_16_UL
INPUTDIR=/eos/user/p/pellicci/MesonGamma_root/2018/HRhoGamma_ggH/DIGI
OUTPUTDIR=/eos/user/p/pellicci/MesonGamma_root/2018/HRhoGamma_ggH/HLT
PYTHONAME=HMesonGamma_HLT_2018_cfg.py

#this is necessary only if EOS access is required
export X509_USER_PROXY=/afs/cern.ch/user/p/pellicci/voms_proxy/x509up_u28550

OFFSET=$2
if [ "$2" == "" ]; then
 OFFSET=0
fi

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

#export SCRAM_ARCH=slc7_amd64_gcc820
if [ -r $CMSSW_TO_USE/src ] ; then 
 echo release $CMSSW_TO_USE already exists
else
scram p CMSSW $CMSSW_TO_USE
fi
cd $CMSSW_TO_USE/src
eval `scram runtime -sh`
echo "check ld_library_path = $LD_LIBRARY_PATH"

cp $HOMEDIR/$PYTHONAME .
scram b

cat << EOF >> $PYTHONAME

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()
EOF

jobNumber=$(($1+$OFFSET));

xrdcp $INPUTDIR/process_${jobNumber}.root process_IN.root
cmsRun $PYTHONAME
xrdcp process.root $OUTPUTDIR/process_${jobNumber}.root
