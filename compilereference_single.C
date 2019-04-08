#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis_llbb.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <string>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <TDCacheFile.h>

#endif

using namespace std;

int main(int argc, char ** argv){
  

  string dataconf=argv[1];
  cout << dataconf << " dataconfiguration" <<endl;

  string mcconf=argv[2];
  cout << mcconf << " mc configuration" <<endl;
  
  Char_t nome[300];

//    sprintf(nome,"/eos/uscms/store/user/wangz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-noMuScale-v1/170718_172042/DYJetsToLL.root");
//  sprintf(nome,"/eos/uscms/store/user/wangz/DoubleEG/DoubleEG_Run2016H-03Feb2017-v1_part1/170717_204041/DoubleEG_H.root");
//    sprintf(nome,"/eos/uscms/store/user/wangz/DoubleMuon/DoubleMuon_Run2016G-03Feb2017-MuCal-v1_part1/170720_230940/DoubleMuon_G_part1.root");
    sprintf(nome,"/eos/uscms/store/user/wangz/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v1/170725_230003/TT_TuneCUETP8M1.root");
//  sprintf(nome,"/eos/uscms/store/user/wangz/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-MuCal-v1/170802_212239/DYBB.root");  

  cout << "test01" << endl;
  TFile *file3;
  file3 = TFile::Open(nome);
  
  cout << "test02" << endl;
  TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
  cout << "test03" << endl;
  HZZ4LeptonsAnalysis make3(tree3,1.,dataconf,mcconf);
  //HZZ4LeptonsAnalysis make3(tree3);

  sprintf(nome,"output_TT_TuneCUETP8M1.root");

  make3.Loop(nome);

  cout << "Create file with name: " << nome << endl;
  delete tree3;
  file3 -> Close();

  return 0; 

}

