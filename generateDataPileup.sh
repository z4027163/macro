#!/bin/bash
#version=`python -c "from DevTools.Utilities.utilities import getCMSSWVersion; print getCMSSWVersion()"`
runPeriod="Collisions16"
central=69200
#lumimask=`python -c "from DevTools.Utilities.utilities import getJson; print getJson('$runPeriod')"`
lumimask="/uscms/home/zwang4/data_lumi/BCDEF.json"
#normtag=`python -c "from DevTools.Utilities.utilities import getNormtag; print getNormtag('$runPeriod')"`
pileupjson="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/$runPeriod/13TeV/PileUp/pileup_latest.txt"
mkdir -p pileup

# new (doesnt work)
#brilcalc lumi -b "STABLE BEAMS" -i $lumimask -o pileup.txt --normtag $normtag --byls --minBiasXsec 71000

maxBins=90

# old pileup
for xsec in $central; do
    up=$(echo "$xsec*1.05" | bc)
    down=$(echo "$xsec*0.95" | bc)
    echo $xsec
    pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $xsec --maxPileupBin $maxBins --numPileupBins $maxBins pileup/PileUpData_BCDEF.root
    echo $up
    pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $up --maxPileupBin $maxBins --numPileupBins $maxBins pileup/PileUpData_up.root
    echo $down
    pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $down --maxPileupBin $maxBins --numPileupBins $maxBins pileup/PileUpData_down.root
done
#for xsec in 60000 61000 62000 63000 64000 65000 66000 67000 68000 69000 70000 71000 72000 73000 74000 75000 76000 77000 78000 79000 80000; do
#    echo $xsec
#    pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $xsec --maxPileupBin $maxBins --numPileupBins $maxBins pileup/PileUpData_$xsec.root
#done
