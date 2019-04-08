#! /usr/bin/env python

#@ CONFIGURATION                                                                                                                                                                                                                              
from optparse import OptionParser
parser = OptionParser()

parser.add_option('--files', type='string', action='store',
                  dest='files',
                  help='Input files')

parser.add_option('--maxruns', type='int', action='store',
                  default=-1,
                  dest='maxruns',
                  help='Number of runs to run. -1 is all runs')

(options, args) = parser.parse_args()
argv = []


#@ FWLITE STUFF

import ROOT
import sys
from DataFormats.FWLite import Runs, Handle
#ROOT.gROOT.Macro("rootlogon.C")
import copy
from array import array
from math import *
ROOT.gSystem.Load("libAnalysisPredictedDistribution")

#@ Labels and Handles
h_lheRun = Handle("LHERunInfoProduct")
l_lheRun = ("externalLHEProducer", "")

    
#@ RUN LOOP

#filelist = file( options.files )
#filesraw = filelist.readlines()

files = []
files.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/045AF4D5-8CC1-E611-9F87-008CFA197D18.root')
#files.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/0A341391-FFD5-E611-9BB7-0CC47A78A3E8.root')
#files.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/110000/024C1DAE-6303-E711-B126-C45444922C46.root')
#files.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/50000/16E7A136-1FBE-E611-B8DA-0025905C3DD0.root')
#('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/FE8A7852-66E4-E611-B5D0-002590E7E01A.root ')
#files.append('root://xrootd.unl.edu//store/mc/RunIISpring16MiniAODv2/RSGluonToTT_M-1000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/20000/02E562CD-9038-E611-A654-B083FED13803.root')
#files.append('root://xrootd.unl.edu//store/mc/RunIISpring16MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/001AFDCE-C33B-E611-B032-0025905D1C54.root')
#files.append('root://cmsxrootd.fnal.gov///store/user/jdolen/B2G2016/ZprimeToTT_M-3000_W-30_TuneCUETP8M1_13TeV-madgraphMLM_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_MINIAOD.root')
#files.append('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/ZprimeToTT_M-3000_W-900_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/02C3658C-8B08-E511-9624-AC853D9DACE3.root')
Nruns = 0
#for ifile in filesraw : #{ Loop over text file and find root files linked
    #if len( ifile ) > 2 :
        #s = ifile.rstrip()
        #files.append( s )
        #print 'Added ' + s
        #} End loop over txt file

# loop over files
for ifile in files : #{ Loop over root files
    print 'Processing file ' + ifile
    runs = Runs (ifile)
    if options.maxruns > 0 and Nruns > options.maxruns :
        break

    # Make sure the handles we want are in the files so we can
    # avoid leaking memory
    readFilters = True

    
    #NrunsInFile = runs.size()
    #if NrunsInFile <= 0 : 
        #continue
    # loop over runs in this file
    for run in runs: #{ Loop over runs in root files

        #Event weight errors
        #gotLHE = event.getByLabel( l_lhe, h_lhe )
        #lhe = h_lhe.product()
        #weight_id = lhe.weights()[0].wgt                                                                                                                                                                                                
            #print '***********' + str(weight_id)                                                                                                                                                                                             
            ## if not gotGenerator or not gotLHE :                                                                                                                                                                                            
            ##     continue                   
        
        gotLHErun = run.getByLabel( l_lheRun, h_lheRun )
        lheRun = h_lheRun.product() 

        

        for i_lheRun in lheRun:
            #print lheRun.headers()[i_lheRun].tag
            print i_lheRun

