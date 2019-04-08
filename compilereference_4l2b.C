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

//    sprintf(nome,"/eos/uscms/store/user/wangz/ZZTo4L_13TeV-amcatnloFXFX-pythia8/ZZTo4L_13TeV_amcatnloFXFX_pythia8_v1_part1/170703_233157/ZZTo4L.root");
//  sprintf(nome,"/eos/uscms/store/user/wangz/DoubleMuon/DoubleMuon_Run2016C-03Feb2017-v2_part1/170628_230635/DoubleMuon.root");
    sprintf(nome,"/eos/uscms/store/user/wangz/DoubleMuon/DoubleMuon_Run2016E-03Feb2017-v2/170714_163811/DoubleMuon_E_par2.root");

  cout << "test01" << endl;
  TFile *file3;
  file3 = TFile::Open(nome);
  
  cout << "test02" << endl;
  TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
  cout << "test03" << endl;
  HZZ4LeptonsAnalysis make3(tree3,1.,dataconf,mcconf);
  //HZZ4LeptonsAnalysis make3(tree3);

  sprintf(nome,"output_DoubleMuon_E_part2.root");
  make3.Loop(nome);

  cout << "Create file with name: " << nome << endl;
  delete tree3;
  file3 -> Close();

  return 0; 

}

