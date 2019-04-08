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
void test(){
/*
  TFile *f = new TFile("fits.root");
  TGraph *g = (TGraph*)f->Get("ratio_eff_eta3_dr030e030_corr");
  TH1F *h = g->GetHistogram();
  int ix = h->GetXaxis()->FindBin(1.4);
  double x = h->GetBinContent(ix);
  double y = g->Eval(1.4);
  std::cout << "test=" << y << " test2=" << ix << std::endl;
  TCanvas *c1 = new TCanvas("c1","A");
  c1->cd();
  g->Draw("AC*");*/

  TFile *f2 = new TFile("IDSF_GH.root");
  TH2F *h2 = (TH2F*)gDirectory->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  TCanvas *c1 = new TCanvas("c1","A");
  c1->cd();
  h2->Draw("COLZ");
  double Pt=50;
  double Eta=0.1;
  int biny = h2->GetYaxis()->FindBin(Pt);
  int binx = h2->GetXaxis()->FindBin(Eta);
  cout << "xbin= " << binx <<" ybin="<<  biny << endl;
  cout << "eff_weight test = " <<h2->GetBinContent(binx,biny)<< endl;

  TFile *f3 = new TFile("ISOSF_GH.root");
  TH2F *h3 = (TH2F*)gDirectory->Get("LooseISO_LooseID_pt_eta/abseta_pt_ratio");
  TCanvas *c2 = new TCanvas("c2","B");
  c2->cd();
  h3->Draw("COLZ");

  int biny3 = h3->GetYaxis()->FindBin(Pt);
  int binx3 = h3->GetXaxis()->FindBin(Eta);
  cout << "xbin= " << binx3 <<" ybin="<<  biny3 << endl;
  cout << "eff_weight test = " <<h3->GetBinContent(binx3,biny3)<< endl;

  TFile *f4 = new TFile("SingleMuonSF_GH.root");
  TH2F *h4 = (TH2F*)gDirectory->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
  TCanvas *c3 = new TCanvas("c3","C");
  c3->cd();
  h4->Draw("COLZ");

  int biny4 = h4->GetYaxis()->FindBin(Pt);
  int binx4 = h4->GetXaxis()->FindBin(Eta);
  cout << "xbin= " << binx4 <<" ybin="<<  biny4 << endl;
  cout << "eff_weight test = " <<h4->GetBinContent(binx4,biny4)<< endl;


}
