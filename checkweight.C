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
#define mPI 3.14159
using namespace std;
double DELTAPHI( double phi1, double phi2 ){

        if( phi1 > mPI || phi1 < -mPI || phi2 > mPI || phi2 < -mPI) {
          // cout << "Angles out of range!!! " << endl;
          // cout << " phi1 " << phi1 << endl;
          // cout << " phi2 " << phi2 << endl;
          return -999;
        }
        float dp=std::abs(phi1-phi2);
        if (dp>mPI) dp-=float(2*mPI);
        return dp;
        //return  min( fabs( phi1 - phi2 ) , 2*mPI - fabs( phi1 - phi2 ) ) ;

}

void checkweight(){

   TString infilename=
//"/eos/uscms/store/user/wangz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v8/180913_180905/345.root";
//"/eos/uscms/store/user/wangz/GluGluTo4L/GluGluToContinToZZTo4mu_13TeV_DefaultShower_MCFM701_pythia8/GluGluToContinToZZTo4mu_13TeV_DefaultShower_MCFM701_pythia8_v1/180720_094516/0.root";
//"/eos/uscms/store/user/wangz/ZZTo4L_13TeV-amcatnloFXFX-pythia8/ZZTo4L_13TeV-amcatnloFXFX-pythia8_v6/180720_093746/ZZTo4L.root";
"/eos/uscms/store/user/wangz/ZZTo4L_13TeV_powheg_pythia8_ext1/ZZTo4L_13TeV_powheg_pythia8_ext1/190207_150025/0.root";
//"/eos/uscms/store/user/wangz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_dm_mcweight_noMuCal_v4/171128_184521/012.root";
   TFile *f1 = new TFile(infilename);
   TTree *t1 = (TTree*)f1->Get("HZZ4LeptonsAnalysis");

   TH1D * hNjets_8 = new TH1D("hNjets_0", "Number of jets passing VBF", 5, 0, 5);
   hNjets_8->SetXTitle("# n-jets");

   TH1F *hpt = new TH1F("hpt", "jet pt", 100, 0, 200);
   TH1F *heta = new TH1F("hpt", "jet eta", 100, -5, 5);
   float mc_weight;
   float temp=0;

   Float_t         MC_GENJET_PT[100];
   Float_t         MC_GENJET_ETA[100];
   Float_t         MC_GENJET_PHI[100];
   bool          dm_trig;
   UChar_t         RECOMU_isPFMu[100];
   UChar_t         RECOMU_isMedium[100];
   Float_t         RECOMU_muInnertrkvalidFraction[100];
   UChar_t         RECOMU_isGlobalMu[100];
   UChar_t         RECOMU_isStandAloneMu[100];
   UChar_t         RECOMU_isTrackerMu[100];
   UChar_t         RECOMU_isCaloMu[100];
   UChar_t         RECOMU_isTrackerHighPtMu[100];
   Float_t         RECOMU_E[100];
   Float_t         RECOMU_PT[100];
   Float_t         RECOMU_P[100];
   Float_t         RECOMU_ETA[100];
   Float_t         RECOMU_THETA[100];
   Float_t         RECOMU_PHI[100];
   Float_t         RECOMU_MASS[100];
   Float_t         RECOMU_CHARGE[100];
   Float_t         RECOMU_mubesttrkDxy[100];
   Int_t           RECO_NMU;
   Float_t         RECOMU_mubesttrkDz[100];
   Int_t           RECO_NELE;
   Double_t        RECOELE_mvaTrigV0[100];
   Float_t         RECOELE_PT[100];
   Float_t         RECOELE_ETA[100];
   Float_t         RECOELE_gsftrack_dxy[100];
   Float_t         RECOELE_gsftrack_dz[100];
   Float_t         RECOELE_scl_Eta[100];

   t1->SetBranchAddress("dm_trig",&dm_trig);
   t1->SetBranchAddress("RECOMU_isPFMu", RECOMU_isPFMu);
   t1->SetBranchAddress("RECOMU_isMedium", RECOMU_isMedium);
   t1->SetBranchAddress("RECOMU_muInnertrkvalidFraction", RECOMU_muInnertrkvalidFraction);
   t1->SetBranchAddress("RECOMU_isGlobalMu", RECOMU_isGlobalMu);
   t1->SetBranchAddress("RECOMU_isStandAloneMu", RECOMU_isStandAloneMu);
   t1->SetBranchAddress("RECOMU_isTrackerMu", RECOMU_isTrackerMu);
   t1->SetBranchAddress("RECOMU_isCaloMu", RECOMU_isCaloMu);
   t1->SetBranchAddress("RECOMU_isTrackerHighPtMu", RECOMU_isTrackerHighPtMu);
   t1->SetBranchAddress("RECOMU_E", RECOMU_E);
   t1->SetBranchAddress("RECOMU_PT", RECOMU_PT);
   t1->SetBranchAddress("RECOMU_P", RECOMU_P);
   t1->SetBranchAddress("RECOMU_ETA", RECOMU_ETA);
   t1->SetBranchAddress("RECOMU_THETA", RECOMU_THETA);
   t1->SetBranchAddress("RECOMU_PHI", RECOMU_PHI);
   t1->SetBranchAddress("RECOMU_MASS", RECOMU_MASS);
   t1->SetBranchAddress("RECOMU_CHARGE", RECOMU_CHARGE);
   t1->SetBranchAddress("RECOMU_mubesttrkDxy", RECOMU_mubesttrkDxy);
// t1->SetBranchAddress("RECOMU_mubesttrkDz", RECOMU_mubesttrkDZ);
   t1->SetBranchAddress("RECOMU_mubesttrkDz", RECOMU_mubesttrkDz);
   t1->SetBranchAddress("RECO_NMU", &RECO_NMU);

   t1->SetBranchAddress("RECO_NELE", &RECO_NELE);
   t1->SetBranchAddress("RECOELE_mvaTrigV0", RECOELE_mvaTrigV0);
   t1->SetBranchAddress("RECOELE_PT", RECOELE_PT);
   t1->SetBranchAddress("RECOELE_ETA", RECOELE_ETA);
   t1->SetBranchAddress("RECOELE_gsftrack_dxy", RECOELE_gsftrack_dxy);
   t1->SetBranchAddress("RECOELE_gsftrack_dz", RECOELE_gsftrack_dz);

   t1->SetBranchAddress("MC_weighting", &mc_weight);
   t1->SetBranchAddress("MC_GENJET_PT",MC_GENJET_PT);
   t1->SetBranchAddress("MC_GENJET_ETA",MC_GENJET_ETA);

   int entries = t1->GetEntries();
   for(int i=0; i < 100000; i++){
       if(i%10000==0) cout << i << endl;
       t1->GetEntry(i);
//       if(i>0){ if(abs(mc_weight)!=abs(temp)) cout << "different " << mc_weight << " and " << temp << endl;} 
//       temp=mc_weight;

//       if(!dm_trig) continue;
      
      int nmu=0;
      for( int i = 0; i < RECO_NMU; ++i ){
        if(
              ( RECOMU_isGlobalMu[i] || RECOMU_isTrackerMu[i] ) && RECOMU_isPFMu[i]
              && RECOMU_PT[i] > 10.
              && fabs(RECOMU_ETA[i]) < 2.4
              && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1. //loose muon (qier) 
            ){
              nmu++;
            }
       }

      int nele=0;
      for( int i = 0; i < RECO_NELE; ++i ){
        if( RECOELE_PT[i] > 10. && fabs(RECOELE_ETA[i]) < 2.5 );
        else continue ;
        bool BDT_ok = 0;
        if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaTrigV0[i] > 0.837 ) BDT_ok = 1 ;
        if( fabs(RECOELE_scl_Eta[i]) >= .8 &&fabs(RECOELE_scl_Eta[i])<1.479 && RECOELE_mvaTrigV0[i] > 0.715 ) BDT_ok = 1 ;
        if(fabs(RECOELE_scl_Eta[i])>=1.479 && RECOELE_mvaTrigV0[i] > 0.357) BDT_ok=1;
        if( !BDT_ok ) continue ;
        if( fabs(RECOELE_gsftrack_dxy[i]) < .5
         && fabs(RECOELE_gsftrack_dz[i])  < 1. ) /* ok */ ;
        else continue ;
     
        nele++;
      }

       if(nmu<4) continue;
  
       double Ngenjet=0;
        for(int k=0;k<100;k++){
           if(MC_GENJET_PT[k]<-1) continue;
           if(abs(MC_GENJET_ETA[k])<2.4 && MC_GENJET_PT[k]>30) Ngenjet++;
           if(abs(MC_GENJET_ETA[k])<2.4&& MC_GENJET_PT[k]>30) {hpt->Fill(MC_GENJET_PT[k]);heta->Fill(MC_GENJET_ETA[k]);}
        }
//       cout << "GenJet number= " << Ngenjet << endl;
         //k factor from SMP-17-005 
/*       if(Ngenjet==0) newweight=newweight*1.20;
       if(Ngenjet==1) newweight=newweight*0.935;
       if(Ngenjet==2) newweight=newweight*0.75;
       else newweight=newweight*0.765;
*/
    hNjets_8->Fill(Ngenjet);
   }
  TCanvas *c1 = new TCanvas("c1","c1",600,800);
//  c1->Divide(1,2);
  c1->cd();

  hNjets_8->Draw();
//  c1->cd(2);
//  hpt->Draw();
  c1->SaveAs("tem.png");
  
} 
