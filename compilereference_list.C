#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis_4l2b.h"
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
#include <TString.h>
#include <sstream>

#endif

using namespace std;

int main(int argc, char ** argv){
  

  string dataconf=argv[1];
  cout << dataconf << " dataconfiguration" <<endl;

  string mcconf=argv[2];
  cout << mcconf << " mc configuration" <<endl;
  
  ifstream fbkg;
  fbkg.open(argv[3]);

  Char_t nome[300];

  string in;
  string out;
  string s;

   while(getline(fbkg,s)){

    stringstream ss(s);
    ss >> in >> out;

    cout << "input= " << in << endl;
    TFile *file3;
    file3 = TFile::Open(in.c_str());
  
    TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
    HZZ4LeptonsAnalysis make3(tree3,1.,dataconf,mcconf);
    //HZZ4LeptonsAnalysis make3(tree3);
    sprintf(nome,out.c_str());
    make3.Loop(nome);

    cout << "Create file with name: " << nome << endl;
    delete tree3;
    file3 -> Close();

  }

  return 0; 

}

