#!/bin/bash

echo "Starting  job"

#setup environment
module() { eval `/usr/bin/modulecmd bash $*`; }
export -f module

MODULESHOME=/usr/share/Modules
export MODULESHOME

if [ "${LOADEDMODULES:-}" = "" ]; then
  LOADEDMODULES=
  export LOADEDMODULES
fi

if [ "${MODULEPATH:-}" = "" ]; then
  MODULEPATH=`sed -n 's/[       #].*$//; /./H; $ { x; s/^\n//; s/\n/:/g; p; }' ${MODULESHOME}/init/.modulespath`
  export MODULEPATH
fi

if [ ${BASH_VERSINFO:-0} -ge 3 ] && [ -r ${MODULESHOME}/init/bash_completion ]; then
 . ${MODULESHOME}/init/bash_completion
fi
source /afs/desy.de/user/c/chayanit/.profile

cd /nfs/dust/cms/user/chayanit/MSSMHBB_Myfork2017/CMSSW_9_4_4/src/
eval `scramv1 runtime -sh`
cmsenv

export LOCAL_DIR=/nfs/dust/cms/user/chayanit/MSSMHBB_Myfork2017/CMSSW_9_4_4/src/Analysis/Mssmhbb/test

$LOCAL_DIR/MssmHbbAnalysis -c examples/allhadronic/mssmhbb_onlcsv2017.cfg

echo "Job finished"
