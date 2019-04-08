#if !defined(__CINT__) || defined(__MAKECINT__)

#include "MonoHiggsAnalysis4e.h"
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
  
  string site=argv[1];
  cout << site << " configuration" <<endl;

  string dataconf=argv[2];
  cout << dataconf << " dataconfiguration" <<endl;

  string mcconf=argv[3];
  cout << mcconf << " mc configuration" <<endl;
  
  Char_t nome[300];

  if (site.find("CERN")<5){
    sprintf(nome,"/castor/cern.ch/user/n/ndefilip/DAS2013/bkg/roottree_leptons_ZZTo2e2mu_8TeV-powheg-pythia6.root"); 
  }
  else if (site.find("FNAL")<5){
    sprintf(nome,"dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/cmsdas/2013/HZZ4lExercise/sig/roottree_leptons_GluGluToHToZZTo4L_M-1000_8TeV-powheg-pythia6.root");
  }
  else {
    // sprintf(nome,"/lustre/cms/store/user/defilip/MonoHiggs/Spring15_25ns_merged/roottree_leptons_MZP1000_MA0300.root");
    //sprintf(nome,"roottree_leptons_vector.root");
    sprintf(nome,"/localdata/Syncr13TeV/roottree_leptons_sync_Fall15_HiggsToZZ_76x_vector.root");
  }


  TFile *file3;
  file3 = TFile::Open(nome);
  
  //if (site.find("FNAL")<5){
  //  file3 = new  TDCacheFile (nome,"READ","ROOT file with trees",0);
  //}
  //else {
  //  file3 = TFile::Open(nome);
  //}

  cout << "Read file with name: " << nome << endl;
  TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
  cout << "Tree has " << tree3->GetEntries() << endl;
  MonoHiggsAnalysis4e make3(tree3,1.,dataconf,mcconf);
  //MonoHiggsAnalysis4mu make3(tree3);

  sprintf(nome,"output_Fall15_Higgs_vector_4e.root");
  make3.Loop(nome);

  cout << "Create file with name: " << nome << endl;
  delete tree3;
  file3 -> Close();

  return 0; 

}

