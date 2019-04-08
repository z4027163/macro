#include <TMath.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>

#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>

void copyhist(){

  TString dir="";

  TString file[1]={"output_WJetsToLNu_amcatnloFXFX_012.root"};

  TString hist_h[1]={"hMZ1_5"};
  TString hist_light[1]={"hMZ1_5"};

  for(int i=0; i<1; i++){
    TFile *f1 = new TFile(dir+file[i],"UPDATE");

    for(int j=0; j<1; j++){
     TH1F *h1 = (TH1F*)f1->Get(hist_h[j]);
     TH1F *h2 = (TH1F*)h1->Clone();
     h2->Write("hMZ1_5_mc_up");
     h2->Write("hMZ1_5_mc_dow");
    }
    f1->Close();
  }

} 

