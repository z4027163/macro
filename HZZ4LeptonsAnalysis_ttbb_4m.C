#define HZZ4LeptonsAnalysis_cxx
#include "HZZ4LeptonsAnalysis_4l2b.h"
#include <TH2.h>
#include <TStyle.h>
//#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TSpline.h>
#include <TRandom3.h>

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
#include <libgen.h>

#include "ZZMatrixElement/MELA/src/computeAngles.h"
#include "ZZMatrixElement/MELA/src/computeAngles.cc"
#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "RoccoR.cc"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

using namespace std;
// using namespace RooFit;
// using namespace meMCFM;
using namespace MEMNames;
using namespace JME;

const double Zmass = 91.188; // nominal Z boson mass
const double mPI = 3.141592654; 

void HZZ4LeptonsAnalysis::Loop(Char_t *output)
{

   if (fChain == 0) return;

   TString pufile, puhist,ele_leg1,ele_leg2,id_sf,iso_sf,mu17_leg,mu8_leg,mu_track;
   string btagcal;
   string era="BF";
   bool botsf=false;

   if(era=="BF"){
     pufile="pileup_BCDEF.root";
     puhist="pileup_scale_BCDEF";
     ele_leg1="Leg1_BF_EGM2D.root";
     ele_leg2="Leg2_BF_EGM2D.root";
     id_sf="IDSF_BCDEF.root";
     iso_sf="ISOSF_BCDEF.root";
     mu17_leg="Mu17Leg_SF_BF.root";
     mu8_leg="Mu8Leg_SF_BF.root";
     mu_track="track_BCDEF.root";
     btagcal="CSVv2_Moriond17_B_F.csv";
   }
   if(era=="GH"){
     pufile="pileup_GH.root";
     puhist="pileup_scale";
     ele_leg1="Leg1_GH_EGM2D.root";
     ele_leg2="Leg2_GH_EGM2D.root";
     id_sf="IDSF_GH.root";
     iso_sf="ISOSF_GH.root";
     mu17_leg="Mu17Leg_SF_GH.root";
     mu8_leg="Mu8Leg_SF_GH.root";
     mu_track="track_GH.root";
     btagcal="CSVv2_Moriond17_G_H.csv";
   }

   TString puhist_up,puhist_dow;
   puhist_up ="pileup_scale_up";
   puhist_dow="pileup_scale_down";

   RoccoR  rc("/uscms/home/zwang4/nobackup/WORKSPCACE/ntuple/CMSSW_8_0_24/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros/roccor/rcdata.2016.v3"); 

// setup calibration + reader
   BTagCalibration calib("CSVv2", "btag/"+btagcal);
   BTagCalibrationReader reader(BTagEntry::OP_LOOSE,"central",{"up","down"});      // other sys types
   reader.load(calib,BTagEntry::FLAV_B,"comb");  
   reader.load(calib,BTagEntry::FLAV_C,"comb");
   reader.load(calib,BTagEntry::FLAV_UDSG,"incl");

   // Declare MEM class
   MEMs combinedMEM(13,125,"CTEQ6L");     
   
    // JME
   JME::JetParameters jetparameters;
   JME::JetResolution jetresolution;
   JME::JetResolutionScaleFactor jetresolution_sf;

   // BNN
   Char_t datasetChar[500],bnnOUT[500],eventsOUT[500];
  
   cout << "The output file is " << output << endl;
   TString out = output;
   TString datasetName=out.ReplaceAll(".root","");
   sprintf(datasetChar,"%s",datasetName.Data());
   sprintf(bnnOUT,"%s_bnn.txt",datasetName.Data());
   sprintf(eventsOUT,"%s_bnn.root",datasetName.Data());

   // Pileup reweighting in 80x
   TFile *_filePU;
   _filePU= TFile::Open("pileup/"+pufile);
   TH1D *puweight = (TH1D*)_filePU->Get(puhist);

   TH1D *puweight_up = (TH1D*)_filePU->Get(puhist_up);
   TH1D *puweight_dow = (TH1D*)_filePU->Get(puhist_dow);
   /////////////Lepton Efficiency Scale Factrons/////////////
   // Load histograms
   //
   TFile *ele_scale_factors_v3 = new TFile("SF_ELE/egammaEffi_RECO_EGM2D.root");
   TH2F *ele_scale_factors_reco = (TH2F*)gDirectory->Get("EGamma_SF2D");
   TFile *ele_scale_factors_v4 = new TFile("SF_ELE/egammaEffi_WP90_EGM2D.root");
   TH2F *ele_scale_factors_wp90 = (TH2F*)gDirectory->Get("EGamma_SF2D");

   TFile *ele_scale_factors_v1 = new TFile("SF_ELE/"+ele_leg1);
   TH2F *ele_scale_factors_leg1 = (TH2F*)gDirectory->Get("EGamma_SF2D");

   TFile *ele_scale_factors_v2 = new TFile("SF_ELE/"+ele_leg2);
   TH2F *ele_scale_factors_leg2 = (TH2F*)gDirectory->Get("EGamma_SF2D");


  TFile *mu_scale_factors3_p2 = new TFile("SF_GH/dm2/"+mu8_leg);
  TH2F *mu_scale_factors_hlt_p2 = (TH2F*)gDirectory->Get("abseta_pt_PLOT");

  TFile *mu_scale_factors1_p1 = new TFile("SF_GH/"+id_sf);
  TH2F *mu_scale_factors_id_p1 = (TH2F*)gDirectory->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile *mu_scale_factors2_p1 = new TFile("SF_GH/"+iso_sf);
  TH2F *mu_scale_factors_iso_p1 = (TH2F*)gDirectory->Get("LooseISO_LooseID_pt_eta/abseta_pt_ratio");

  TFile *mu_scale_factors3_p1 = new TFile("SF_GH/dm2/"+mu17_leg);
  TH2F *mu_scale_factors_hlt_p1 = (TH2F*)gDirectory->Get("abseta_pt_PLOT");

  TFile *mu_scale_factors4 = new TFile("SF_GH/"+mu_track); //just for GH
  TGraph *mu_scale_factors_tk = (TGraph*)gDirectory->Get("ratio_eff_eta3_dr030e030_corr");

  TFile *b_scale_factors3_p1 = new TFile("btag/mc_4l_eff.root");
  TH1F *b_eff_p1 = (TH1F*)gDirectory->Get("hPtBot_8_b");
  TH1F *b_eff_p2 = (TH1F*)gDirectory->Get("hPtBot_8_c");
  TH1F *b_eff_p4 = (TH1F*)gDirectory->Get("hPtBot_8_o");

  TFile *b_scale_factors3_p3 = new TFile("btag/QCD_all.root");
  TH1F *b_eff_p3 = (TH1F*)gDirectory->Get("hPtBot_8_l");

  TFile *b_scale_factors3_p4 = new TFile("btag/l_jet_sf.root");
  TH1F *b_eff_zz_p5 = (TH1F*)gDirectory->Get("hPtBot_8_l_zz");
   
   // Book root file (for output):
   TFile * theFile = new TFile(output,"RECREATE");

    // Clone tree for Z1
   double eta_jer[10]={0.0,0.522,0.783,1.131,1.305,1.740,1.930,2.043,2.322,2.5};
   double un_jer[9]={0.0645,0.0652,0.0632,0.1025,0.0986,0.1079,0.1214,0.1140,0.2371};
   double scale_jer[9]={1.1595,1.1948,1.1464,1.1609,1.1278,1.1000,1.1426,1.1512,1.2963};

   double DELTAPHI( double , double ) ; //call the function  
   double invmass (float M1, float PT1, float ETA1, float PHI1, float M2, float PT2, float ETA2, float PHI2 );
   
   // Book relevant variables -- counters:

   int N_0 = 0;	  // MC truth & acceptance
   int N_01 = 0;
   int N_02 = 0;

   int N_1 = 0;	  // Skim
   int N_2 = 0;

   int N_3a = 0;
   int N_3_FSR = 0;
   int N_3b = 0;

   int N_4a = 0;
   int N_4b = 0;
   int N_4c = 0;
   int N_4d = 0;

   int N_5 = 0;
   int N_6 = 0;
   int N_7 = 0;
   int N_7_PFMET = 0;
    
   int N_8 = 0;
   int N_8_a = 0;
   int N_8_PFMET = 0;
   int N_9 = 0;

   int N_9_1FSR = 0;
   int N_9_2FSR = 0;

   int N_9PS = 0;
   int N_9GRAV = 0;
   
   int N_9a_VBF = 0;
   int N_9b_VBF = 0;
   int N_9_PFMET = 0;

   int N_VBF = 0;
   int N_bjets = 0;
   int N_bjets_cut = 0;
   int N_njets_cut = 0;

   int N_10 = 0;

   // counter weighted
   double N_0_w = 0;	  // MC truth & acceptance
   double N_01_w = 0;
   double N_02_w = 0;

   double N_1_w = 0;	  // Skim
   double N_2_w = 0;

   double N_3a_w = 0;
   double N_3_FSR_w = 0;
   double N_3b_w = 0;

   double N_4a_w = 0;
   double N_4b_w = 0;
   double N_4c_w = 0;
   double N_4d_w = 0;

   double N_5_w = 0;
   double N_6_w = 0;
   double N_7_w = 0;
   double N_7_PFMET_w = 0;


   double N_8_w = 0;
   double N_8_a_w = 0;
   double N_8_PFMET_w = 0;
   double N_9_w = 0;

   double N_9_1FSR_w = 0;
   double N_9_2FSR_w = 0;

   double N_9PS_w = 0;
   double N_9GRAV_w = 0;
   
   double N_9a_VBF_w = 0;
   double N_9b_VBF_w = 0;
   double N_9_PFMET_w = 0;
     
   double N_VBF_w = 0;
   double N_bjets_w = 0;
   double N_bjets_cut_w = 0;
   double N_njets_cut_w = 0;

   int N_10_w = 0;
   //******* SENSITIVITY ON MET
   //////////// END OF SENSITIVITY ON MET



   // Book Histos ***
   TH1D *nEvent_4l_w = new TH1D("nEvent_4l_w", "nEventComplete Weighted", 16, 0., 16.);
   TH1D *nEvent_4l = new TH1D("nEvent_4l", "nEventComplete", 16, 0., 16.);

   //SENSITIVITY
   
   TH1D *nEvent_CUT_w = new TH1D("nEvent_CUT_w", "nEventCUT Weightd", 61, 0., 61.);
   TH1D *nEvent_CUT = new TH1D("nEvent_CUT", "nEventCUT", 61, 0., 61.);
   
   //

   // Pileup reweighting
   TH1F *hPUvertices             = new TH1F("hPUvertices", "hPUvertices",70,0.,70.);  
   TH1F *hPUvertices_ReWeighted  = new TH1F("hPUvertices_ReWeighted", "hPUvertices_ReWeighted",70,0.,70.);  

   //step 3

   TH1F * hMZ1_5 = new TH1F("hMZ1_5", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5 = new TH1F("hPtZ1_5", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5 = new TH1F("hYZ1_5", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5->SetXTitle("Y_Z1");

   TH1F * hMZ1_5_pdf_up = new TH1F("hMZ1_5_pdf_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_5_pdf_up = new TH1F("hPtZ1_5_pdf_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_5_pdf_up = new TH1F("hYZ1_5_pdf_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_5_pdf_dow = new TH1F("hMZ1_5_pdf_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_5_pdf_dow = new TH1F("hPtZ1_5_pdf_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_5_pdf_dow = new TH1F("hYZ1_5_pdf_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_5_pu_up = new TH1F("hMZ1_5_pu_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_5_pu_up = new TH1F("hPtZ1_5_pu_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_5_pu_up = new TH1F("hYZ1_5_pu_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_5_pu_dow = new TH1F("hMZ1_5_pu_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_5_pu_dow = new TH1F("hPtZ1_5_pu_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_5_pu_dow = new TH1F("hYZ1_5_pu_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_5_mc_1 = new TH1F("hMZ1_5_mc_1", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5_mc_1->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5_mc_1 = new TH1F("hPtZ1_5_mc_1", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5_mc_1->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5_mc_1 = new TH1F("hYZ1_5_mc_1", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5_mc_1->SetXTitle("Y_Z1");

   TH1F * hMZ1_5_mc_2 = new TH1F("hMZ1_5_mc_2", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5_mc_2->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5_mc_2 = new TH1F("hPtZ1_5_mc_2", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5_mc_2->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5_mc_2 = new TH1F("hYZ1_5_mc_2", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5_mc_2->SetXTitle("Y_Z1");

   TH1F * hMZ1_5_mc_3 = new TH1F("hMZ1_5_mc_3", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5_mc_3->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5_mc_3 = new TH1F("hPtZ1_5_mc_3", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5_mc_3->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5_mc_3 = new TH1F("hYZ1_5_mc_3", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5_mc_3->SetXTitle("Y_Z1");

   TH1F * hMZ1_5_mc_4 = new TH1F("hMZ1_5_mc_4", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5_mc_4->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5_mc_4 = new TH1F("hPtZ1_5_mc_4", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5_mc_4->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5_mc_4 = new TH1F("hYZ1_5_mc_4", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5_mc_4->SetXTitle("Y_Z1");

   TH1F * hMZ1_5_mc_5 = new TH1F("hMZ1_5_mc_5", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5_mc_5->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5_mc_5 = new TH1F("hPtZ1_5_mc_5", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5_mc_5->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5_mc_5 = new TH1F("hYZ1_5_mc_5", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5_mc_5->SetXTitle("Y_Z1");

   TH1F * hMZ1_5_mc_6 = new TH1F("hMZ1_5_mc_6", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5_mc_6->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5_mc_6 = new TH1F("hPtZ1_5_mc_6", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5_mc_6->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5_mc_6 = new TH1F("hYZ1_5_mc_6", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5_mc_6->SetXTitle("Y_Z1");

   TH1F * hMZ1_5_mc_7 = new TH1F("hMZ1_5_mc_7", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5_mc_7->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5_mc_7 = new TH1F("hPtZ1_5_mc_7", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5_mc_7->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5_mc_7 = new TH1F("hYZ1_5_mc_7", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5_mc_7->SetXTitle("Y_Z1");

   TH1F * hMZ1_5_mc_8 = new TH1F("hMZ1_5_mc_8", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5_mc_8->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5_mc_8 = new TH1F("hPtZ1_5_mc_8", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5_mc_8->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5_mc_8 = new TH1F("hYZ1_5_mc_8", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5_mc_8->SetXTitle("Y_Z1");

   TH1F * hMZ2_5 = new TH1F("hMZ2_5", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5 = new TH1F("hPtZ2_5", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5 = new TH1F("hYZ2_5", "Y of Z2 after selection step 5", 500 , -5. , 5. );
   hYZ2_5->SetXTitle("Y_Z2");

   TH1F * hMZ2_5_pdf_up = new TH1F("hMZ2_5_pdf_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_5_pdf_up = new TH1F("hPtZ2_5_pdf_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_5_pdf_up = new TH1F("hYZ2_5_pdf_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_5_pdf_dow = new TH1F("hMZ2_5_pdf_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_5_pdf_dow = new TH1F("hPtZ2_5_pdf_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_5_pdf_dow = new TH1F("hYZ2_5_pdf_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_5_pu_up = new TH1F("hMZ2_5_pu_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_5_pu_up = new TH1F("hPtZ2_5_pu_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_5_pu_up = new TH1F("hYZ2_5_pu_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_5_pu_dow = new TH1F("hMZ2_5_pu_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_5_pu_dow = new TH1F("hPtZ2_5_pu_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_5_pu_dow = new TH1F("hYZ2_5_pu_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_5_mc_1 = new TH1F("hMZ2_5_mc_1", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5_mc_1->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5_mc_1 = new TH1F("hPtZ2_5_mc_1", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5_mc_1->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5_mc_1 = new TH1F("hYZ2_5_mc_1", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_5_mc_1->SetXTitle("Y_Z2");

   TH1F * hMZ2_5_mc_2 = new TH1F("hMZ2_5_mc_2", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5_mc_2->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5_mc_2 = new TH1F("hPtZ2_5_mc_2", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5_mc_2->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5_mc_2 = new TH1F("hYZ2_5_mc_2", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_5_mc_2->SetXTitle("Y_Z2");

   TH1F * hMZ2_5_mc_3 = new TH1F("hMZ2_5_mc_3", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5_mc_3->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5_mc_3 = new TH1F("hPtZ2_5_mc_3", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5_mc_3->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5_mc_3 = new TH1F("hYZ2_5_mc_3", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_5_mc_3->SetXTitle("Y_Z2");

   TH1F * hMZ2_5_mc_4 = new TH1F("hMZ2_5_mc_4", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5_mc_4->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5_mc_4 = new TH1F("hPtZ2_5_mc_4", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5_mc_4->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5_mc_4 = new TH1F("hYZ2_5_mc_4", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_5_mc_4->SetXTitle("Y_Z2");

   TH1F * hMZ2_5_mc_5 = new TH1F("hMZ2_5_mc_5", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5_mc_5->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5_mc_5 = new TH1F("hPtZ2_5_mc_5", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5_mc_5->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5_mc_5 = new TH1F("hYZ2_5_mc_5", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_5_mc_5->SetXTitle("Y_Z2");

   TH1F * hMZ2_5_mc_6 = new TH1F("hMZ2_5_mc_6", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5_mc_6->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5_mc_6 = new TH1F("hPtZ2_5_mc_6", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5_mc_6->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5_mc_6 = new TH1F("hYZ2_5_mc_6", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_5_mc_6->SetXTitle("Y_Z2");

   TH1F * hMZ2_5_mc_7 = new TH1F("hMZ2_5_mc_7", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5_mc_7->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5_mc_7 = new TH1F("hPtZ2_5_mc_7", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5_mc_7->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5_mc_7 = new TH1F("hYZ2_5_mc_7", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_5_mc_7->SetXTitle("Y_Z2");

   TH1F * hMZ2_5_mc_8 = new TH1F("hMZ2_5_mc_8", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5_mc_8->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5_mc_8 = new TH1F("hPtZ2_5_mc_8", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5_mc_8->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5_mc_8 = new TH1F("hYZ2_5_mc_8", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_5_mc_8->SetXTitle("Y_Z2");


   TH1F * hPtLep1_5 = new TH1F("hPtLep1_5", "Pt of Lep1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep1_5->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_5 = new TH1F("hEtaLep1_5", "Y of Lep1 after selection step 5", 500 , -5. , 5. );
   hYLep1_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_5 = new TH1F("hIsoLep1_5", "Isolation of Lep1 after selection step 5", 2000 , -10. , 10. );
   hIsoLep1_5->SetXTitle("Iso");
   TH1F * hSipLep1_5 = new TH1F("hSipLep1_5", "Sip of Lep1 after selection step 5",  1000 , -20. , 40. );
   hSipLep1_5->SetXTitle("Sip");
   TH1F * hIpLep1_5 = new TH1F("hIpLep1_5", "Ip of Lep1 after selection step 5",  1000 , -20. , 40. );
   hIpLep1_5->SetXTitle("Ip");
   TH1F * hIpErLep1_5 = new TH1F("hIpErLep1_5", "IpEr of Lep1 after selection step 5",  1000 , 0. , 10. );
   hIpErLep1_5->SetXTitle("IpEr");

   TH1F * hPtLep2_5 = new TH1F("hPtLep2_5", "Pt of Lep2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep2_5->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_5 = new TH1F("hEtaLep2_5", "Y of Lep2 after selection step 5", 500 , -5. , 5. );
   hYLep2_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_5 = new TH1F("hIsoLep2_5", "Isolation of Lep2 after selection step 5", 2000 , -10. , 10. );
   hIsoLep2_5->SetXTitle("Iso");
   TH1F * hSipLep2_5 = new TH1F("hSipLep2_5", "Sip of Lep2 after selection step 5",  1000 , -20. , 40. );
   hSipLep2_5->SetXTitle("Sip");
   TH1F * hIpLep2_5 = new TH1F("hIpLep2_5", "Ip of Lep2 after selection step 5",  1000 , -20. , 40. );
   hIpLep2_5->SetXTitle("Ip");
   TH1F * hIpErLep2_5 = new TH1F("hIpErLep2_5", "IpEr of Lep2 after selection step 5",  1000 , 0. , 10. );
   hIpErLep2_5->SetXTitle("IpEr");
   
   TH1F * hPtLep3_5 = new TH1F("hPtLep3_5", "Pt of Lep3 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep3_5->SetXTitle("pt_Lep3  (GeV)");
   TH1F * hYLep3_5 = new TH1F("hEtaLep3_5", "Y of Lep3 after selection step 5", 500 , -5. , 5. );
   hYLep3_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep3_5 = new TH1F("hIsoLep3_5", "Isolation of Lep3 after selection step 5", 2000 , -10. , 10. );
   hIsoLep3_5->SetXTitle("Iso");
   TH1F * hSipLep3_5 = new TH1F("hSipLep3_5", "Sip of Lep3 after selection step 5",  1000 , -20. , 40. );
   hSipLep3_5->SetXTitle("Sip");
   TH1F * hIpLep3_5 = new TH1F("hIpLep3_5", "Ip of Lep3 after selection step 5",  1000 , -20. , 40. );
   hIpLep3_5->SetXTitle("Ip");
   TH1F * hIpErLep3_5 = new TH1F("hIpErLep3_5", "IpEr of Lep3 after selection step 5",  1000 , 0. , 10. );
   hIpErLep3_5->SetXTitle("IpEr");

   TH1F * hPtLep4_5 = new TH1F("hPtLep4_5", "Pt of Lep4 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep4_5->SetXTitle("pt_Lep4  (GeV)");
   TH1F * hYLep4_5 = new TH1F("hEtaLep4_5", "Y of Lep4 after selection step 5", 50 , -5. , 5. );
   hYLep4_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep4_5 = new TH1F("hIsoLep4_5", "Isolation of Lep4 after selection step 5", 2000 , -10. , 10. );
   hIsoLep4_5->SetXTitle("Iso");
   TH1F * hSipLep4_5 = new TH1F("hSipLep4_5", "Sip of Lep4 after selection step 5",  1000 , -20. , 40. );
   hSipLep4_5->SetXTitle("Sip");
   TH1F * hIpLep4_5 = new TH1F("hIpLep4_5", "Ip of Lep4 after selection step 5",  1000 , -20. , 40. );
   hIpLep4_5->SetXTitle("Ip");
   TH1F * hIpErLep4_5 = new TH1F("hIpErLep4_5", "IpEr of Lep4 after selection step 5",  1000 , 0. , 10. );
   hIpErLep4_5->SetXTitle("IpEr");


   TH1F * hIso_5 = new TH1F("hIso_5", "Isolation maxima after selection step 5", 2000 , -10. , 10. );
   hIso_5->SetXTitle("Iso");
   TH1F * hSip_5 = new TH1F("hSip_5", "Sip maxima after selection step 5",  1000 , -20. , 40. );
   hSip_5->SetXTitle("Sip");
   TH1F * hIp_5 = new TH1F("hIp_5", "Ip maxima after selection step 5",  1000 , -20. , 40. );
   hIp_5->SetXTitle("Ip");
   

   //step 6
   
   TH1F * hminMll_6 = new TH1F("hminMll_6", "minMll at selection step 6", 400 , 0. , 200.);
   hminMll_6->SetXTitle("minMll  (GeV)");

   TH1F * hMZ1_6 = new TH1F("hMZ1_6", "Mass of Z1 after Z1 selection", 200 , -0.5 , 199.5 );
   hMZ1_6->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6 = new TH1F("hPtZ1_6", "Pt of Z1 after Z1 selection", 200 , -0.5 , 199.5 );
   hPtZ1_6->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6 = new TH1F("hYZ1_6", "Y of Z1 after Z1 selection", 500 , -5. , 5.);
   hYZ1_6->SetXTitle("Y_Z1");

   TH1F * hMZ1_6_pdf_up = new TH1F("hMZ1_6_pdf_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_6_pdf_up = new TH1F("hPtZ1_6_pdf_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_6_pdf_up = new TH1F("hYZ1_6_pdf_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_6_pdf_dow = new TH1F("hMZ1_6_pdf_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_6_pdf_dow = new TH1F("hPtZ1_6_pdf_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_6_pdf_dow = new TH1F("hYZ1_6_pdf_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_6_pu_up = new TH1F("hMZ1_6_pu_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_6_pu_up = new TH1F("hPtZ1_6_pu_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_6_pu_up = new TH1F("hYZ1_6_pu_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_6_pu_dow = new TH1F("hMZ1_6_pu_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_6_pu_dow = new TH1F("hPtZ1_6_pu_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_6_pu_dow = new TH1F("hYZ1_6_pu_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_6_mc_1 = new TH1F("hMZ1_6_mc_1", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_6_mc_1->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6_mc_1 = new TH1F("hPtZ1_6_mc_1", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_6_mc_1->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6_mc_1 = new TH1F("hYZ1_6_mc_1", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_6_mc_1->SetXTitle("Y_Z1");

   TH1F * hMZ1_6_mc_2 = new TH1F("hMZ1_6_mc_2", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_6_mc_2->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6_mc_2 = new TH1F("hPtZ1_6_mc_2", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_6_mc_2->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6_mc_2 = new TH1F("hYZ1_6_mc_2", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_6_mc_2->SetXTitle("Y_Z1");

   TH1F * hMZ1_6_mc_3 = new TH1F("hMZ1_6_mc_3", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_6_mc_3->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6_mc_3 = new TH1F("hPtZ1_6_mc_3", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_6_mc_3->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6_mc_3 = new TH1F("hYZ1_6_mc_3", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_6_mc_3->SetXTitle("Y_Z1");

   TH1F * hMZ1_6_mc_4 = new TH1F("hMZ1_6_mc_4", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_6_mc_4->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6_mc_4 = new TH1F("hPtZ1_6_mc_4", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_6_mc_4->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6_mc_4 = new TH1F("hYZ1_6_mc_4", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_6_mc_4->SetXTitle("Y_Z1");

   TH1F * hMZ1_6_mc_5 = new TH1F("hMZ1_6_mc_5", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_6_mc_5->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6_mc_5 = new TH1F("hPtZ1_6_mc_5", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_6_mc_5->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6_mc_5 = new TH1F("hYZ1_6_mc_5", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_6_mc_5->SetXTitle("Y_Z1");

   TH1F * hMZ1_6_mc_6 = new TH1F("hMZ1_6_mc_6", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_6_mc_6->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6_mc_6 = new TH1F("hPtZ1_6_mc_6", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_6_mc_6->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6_mc_6 = new TH1F("hYZ1_6_mc_6", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_6_mc_6->SetXTitle("Y_Z1");

   TH1F * hMZ1_6_mc_7 = new TH1F("hMZ1_6_mc_7", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_6_mc_7->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6_mc_7 = new TH1F("hPtZ1_6_mc_7", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_6_mc_7->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6_mc_7 = new TH1F("hYZ1_6_mc_7", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_6_mc_7->SetXTitle("Y_Z1");


   TH1F * hMZ1_6_mc_8 = new TH1F("hMZ1_6_mc_8", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_6_mc_8->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_6_mc_8 = new TH1F("hPtZ1_6_mc_8", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_6_mc_8->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_6_mc_8 = new TH1F("hYZ1_6_mc_8", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_6_mc_8->SetXTitle("Y_Z1");


   TH1F * hMZ2_6 = new TH1F("hMZ2_6", "Mass of Z2 after jj selection", 200 , -0.5 , 199.5 );
   hMZ2_6->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6 = new TH1F("hPtZ2_6", "Pt of Z2 after jj selection", 200 , -0.5 , 199.5 );
   hPtZ2_6->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6 = new TH1F("hYZ2_6", "Y of Z2 after jj selection", 500 , -5. , 5.);
   hYZ2_6->SetXTitle("Y_Z2");

   TH1F * hMZ2_6_pdf_up = new TH1F("hMZ2_6_pdf_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_6_pdf_up = new TH1F("hPtZ2_6_pdf_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_6_pdf_up = new TH1F("hYZ2_6_pdf_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_6_pdf_dow = new TH1F("hMZ2_6_pdf_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_6_pdf_dow = new TH1F("hPtZ2_6_pdf_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_6_pdf_dow = new TH1F("hYZ2_6_pdf_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_6_pu_up = new TH1F("hMZ2_6_pu_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_6_pu_up = new TH1F("hPtZ2_6_pu_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_6_pu_up = new TH1F("hYZ2_6_pu_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_6_pu_dow = new TH1F("hMZ2_6_pu_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_6_pu_dow = new TH1F("hPtZ2_6_pu_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_6_pu_dow = new TH1F("hYZ2_6_pu_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);
 
   TH1F * hMZ2_6_mc_1 = new TH1F("hMZ2_6_mc_1", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_6_mc_1->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6_mc_1 = new TH1F("hPtZ2_6_mc_1", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_6_mc_1->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6_mc_1 = new TH1F("hYZ2_6_mc_1", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_6_mc_1->SetXTitle("Y_Z2");

   TH1F * hMZ2_6_mc_2 = new TH1F("hMZ2_6_mc_2", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_6_mc_2->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6_mc_2 = new TH1F("hPtZ2_6_mc_2", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_6_mc_2->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6_mc_2 = new TH1F("hYZ2_6_mc_2", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_6_mc_2->SetXTitle("Y_Z2");

   TH1F * hMZ2_6_mc_3 = new TH1F("hMZ2_6_mc_3", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_6_mc_3->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6_mc_3 = new TH1F("hPtZ2_6_mc_3", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_6_mc_3->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6_mc_3 = new TH1F("hYZ2_6_mc_3", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_6_mc_3->SetXTitle("Y_Z2");

   TH1F * hMZ2_6_mc_4 = new TH1F("hMZ2_6_mc_4", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_6_mc_4->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6_mc_4 = new TH1F("hPtZ2_6_mc_4", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_6_mc_4->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6_mc_4 = new TH1F("hYZ2_6_mc_4", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_6_mc_4->SetXTitle("Y_Z2");

   TH1F * hMZ2_6_mc_5 = new TH1F("hMZ2_6_mc_5", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_6_mc_5->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6_mc_5 = new TH1F("hPtZ2_6_mc_5", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_6_mc_5->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6_mc_5 = new TH1F("hYZ2_6_mc_5", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_6_mc_5->SetXTitle("Y_Z2");

   TH1F * hMZ2_6_mc_6 = new TH1F("hMZ2_6_mc_6", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_6_mc_6->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6_mc_6 = new TH1F("hPtZ2_6_mc_6", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_6_mc_6->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6_mc_6 = new TH1F("hYZ2_6_mc_6", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_6_mc_6->SetXTitle("Y_Z2");

   TH1F * hMZ2_6_mc_7 = new TH1F("hMZ2_6_mc_7", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_6_mc_7->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6_mc_7 = new TH1F("hPtZ2_6_mc_7", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_6_mc_7->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6_mc_7 = new TH1F("hYZ2_6_mc_7", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_6_mc_7->SetXTitle("Y_Z2");

   TH1F * hMZ2_6_mc_8 = new TH1F("hMZ2_6_mc_8", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_6_mc_8->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_6_mc_8 = new TH1F("hPtZ2_6_mc_8", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_6_mc_8->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_6_mc_8 = new TH1F("hYZ2_6_mc_8", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_6_mc_8->SetXTitle("Y_Z2");

   //step 7

   TH1F * hMZ1_7 = new TH1F("hMZ1_7", "Mass of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hMZ1_7->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_7 = new TH1F("hPtZ1_7", "Pt of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ1_7->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_7 = new TH1F("hYZ1_7", "Y of Z1 after selection step 7", 500 , -5. , 5.);
   hYZ1_7->SetXTitle("Y_Z1");

   TH1F * hMZ1_7_pdf_up = new TH1F("hMZ1_7_pdf_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_pdf_up = new TH1F("hPtZ1_7_pdf_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_pdf_up = new TH1F("hYZ1_7_pdf_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_pdf_dow = new TH1F("hMZ1_7_pdf_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_pdf_dow = new TH1F("hPtZ1_7_pdf_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_pdf_dow = new TH1F("hYZ1_7_pdf_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_pu_up = new TH1F("hMZ1_7_pu_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_pu_up = new TH1F("hPtZ1_7_pu_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_pu_up = new TH1F("hYZ1_7_pu_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_pu_dow = new TH1F("hMZ1_7_pu_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_pu_dow = new TH1F("hPtZ1_7_pu_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_pu_dow = new TH1F("hYZ1_7_pu_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_mc_1 = new TH1F("hMZ1_7_mc_1", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_mc_1 = new TH1F("hPtZ1_7_mc_1", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_mc_1 = new TH1F("hYZ1_7_mc_1", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_mc_2 = new TH1F("hMZ1_7_mc_2", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_mc_2 = new TH1F("hPtZ1_7_mc_2", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_mc_2 = new TH1F("hYZ1_7_mc_2", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_mc_3 = new TH1F("hMZ1_7_mc_3", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_mc_3 = new TH1F("hPtZ1_7_mc_3", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_mc_3 = new TH1F("hYZ1_7_mc_3", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_mc_4 = new TH1F("hMZ1_7_mc_4", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_mc_4 = new TH1F("hPtZ1_7_mc_4", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_mc_4 = new TH1F("hYZ1_7_mc_4", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_mc_5 = new TH1F("hMZ1_7_mc_5", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_mc_5 = new TH1F("hPtZ1_7_mc_5", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_mc_5 = new TH1F("hYZ1_7_mc_5", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_mc_6 = new TH1F("hMZ1_7_mc_6", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_mc_6 = new TH1F("hPtZ1_7_mc_6", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_mc_6 = new TH1F("hYZ1_7_mc_6", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_mc_7 = new TH1F("hMZ1_7_mc_7", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_mc_7 = new TH1F("hPtZ1_7_mc_7", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_mc_7 = new TH1F("hYZ1_7_mc_7", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ1_7_mc_8 = new TH1F("hMZ1_7_mc_8", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ1_7_mc_8 = new TH1F("hPtZ1_7_mc_8", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ1_7_mc_8 = new TH1F("hYZ1_7_mc_8", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_7 = new TH1F("hMZ2_7", "Mass of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hMZ2_7->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7 = new TH1F("hPtZ2_7", "Pt of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ2_7->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7 = new TH1F("hYZ2_7", "Y of Z2 after selection step 7", 500 , -5. , 5.);
   hYZ2_7->SetXTitle("Y_Z2");

   TH1F * hMZ2_7_pdf_up = new TH1F("hMZ2_7_pdf_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_7_pdf_up = new TH1F("hPtZ2_7_pdf_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_7_pdf_up = new TH1F("hYZ2_7_pdf_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_7_pdf_dow = new TH1F("hMZ2_7_pdf_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_7_pdf_dow = new TH1F("hPtZ2_7_pdf_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_7_pdf_dow = new TH1F("hYZ2_7_pdf_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_7_pu_up = new TH1F("hMZ2_7_pu_up", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_7_pu_up = new TH1F("hPtZ2_7_pu_up", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_7_pu_up = new TH1F("hYZ2_7_pu_up", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_7_pu_dow = new TH1F("hMZ2_7_pu_dow", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hPtZ2_7_pu_dow = new TH1F("hPtZ2_7_pu_dow", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   TH1F * hYZ2_7_pu_dow = new TH1F("hYZ2_7_pu_dow", "Y of Z1 after selection step 5", 500 , -5. , 5.);

   TH1F * hMZ2_7_mc_1 = new TH1F("hMZ2_7_mc_1", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_7_mc_1->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7_mc_1 = new TH1F("hPtZ2_7_mc_1", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_7_mc_1->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7_mc_1 = new TH1F("hYZ2_7_mc_1", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_7_mc_1->SetXTitle("Y_Z2");
   
   TH1F * hMZ2_7_mc_2 = new TH1F("hMZ2_7_mc_2", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_7_mc_2->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7_mc_2 = new TH1F("hPtZ2_7_mc_2", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_7_mc_2->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7_mc_2 = new TH1F("hYZ2_7_mc_2", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_7_mc_2->SetXTitle("Y_Z2");
   
   TH1F * hMZ2_7_mc_3 = new TH1F("hMZ2_7_mc_3", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_7_mc_3->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7_mc_3 = new TH1F("hPtZ2_7_mc_3", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_7_mc_3->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7_mc_3 = new TH1F("hYZ2_7_mc_3", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_7_mc_3->SetXTitle("Y_Z2");
   
   TH1F * hMZ2_7_mc_4 = new TH1F("hMZ2_7_mc_4", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_7_mc_4->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7_mc_4 = new TH1F("hPtZ2_7_mc_4", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_7_mc_4->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7_mc_4 = new TH1F("hYZ2_7_mc_4", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_7_mc_4->SetXTitle("Y_Z2");
   
   TH1F * hMZ2_7_mc_5 = new TH1F("hMZ2_7_mc_5", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_7_mc_5->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7_mc_5 = new TH1F("hPtZ2_7_mc_5", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_7_mc_5->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7_mc_5 = new TH1F("hYZ2_7_mc_5", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_7_mc_5->SetXTitle("Y_Z2");
   
   TH1F * hMZ2_7_mc_6 = new TH1F("hMZ2_7_mc_6", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_7_mc_6->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7_mc_6 = new TH1F("hPtZ2_7_mc_6", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_7_mc_6->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7_mc_6 = new TH1F("hYZ2_7_mc_6", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_7_mc_6->SetXTitle("Y_Z2");
   
   TH1F * hMZ2_7_mc_7 = new TH1F("hMZ2_7_mc_7", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_7_mc_7->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7_mc_7 = new TH1F("hPtZ2_7_mc_7", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_7_mc_7->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7_mc_7 = new TH1F("hYZ2_7_mc_7", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_7_mc_7->SetXTitle("Y_Z2");

   TH1F * hMZ2_7_mc_8 = new TH1F("hMZ2_7_mc_8", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_7_mc_8->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7_mc_8 = new TH1F("hPtZ2_7_mc_8", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_7_mc_8->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7_mc_8 = new TH1F("hYZ2_7_mc_8", "Y of Z2 after selection step 5", 500 , -5. , 5.);
   hYZ2_7_mc_8->SetXTitle("Y_Z2");


   TH1F * hPtLep1_7 = new TH1F("hPtLep1_7", "Pt of Lep1 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep1_7->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_7 = new TH1F("hEtaLep1_7", "Y of Lep1 after selection step 7", 500 , -5. , 5. );
   hYLep1_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_7 = new TH1F("hIsoLep1_7", "Isolation of Lep1 after selection step 7", 2000 , -10. , 10. );
   hIsoLep1_7->SetXTitle("Iso");
   TH1F * hSipLep1_7 = new TH1F("hSipLep1_7", "Sip of Lep1 after selection step 7",  1000 , -20. , 40. );
   hSipLep1_7->SetXTitle("Sip");
   TH1F * hIpLep1_7 = new TH1F("hIpLep1_7", "Ip of Lep1 after selection step 7",  1000 , -20. , 40. );
   hIpLep1_7->SetXTitle("Ip");
   TH1F * hIpErLep1_7 = new TH1F("hIpErLep1_7", "IpEr of Lep1 after selection step 7",  1000 , 0. , 10. );
   hIpErLep1_7->SetXTitle("IpEr");


   TH1F * hPtLep2_7 = new TH1F("hPtLep2_7", "Pt of Lep2 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep2_7->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_7 = new TH1F("hEtaLep2_7", "Y of Lep2 after selection step 7", 500 , -5. , 5. );
   hYLep2_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_7 = new TH1F("hIsoLep2_7", "Isolation of Lep2 after selection step 7", 2000 , -10. , 10. );
   hIsoLep2_7->SetXTitle("Iso");
   TH1F * hSipLep2_7 = new TH1F("hSipLep2_7", "Sip of Lep2 after selection step 7",  1000 , -20. , 40. );
   hSipLep2_7->SetXTitle("Sip");
   TH1F * hIpLep2_7 = new TH1F("hIpLep2_7", "Ip of Lep2 after selection step 7",  1000 , -20. , 40. );
   hIpLep2_7->SetXTitle("Ip");
   TH1F * hIpErLep2_7 = new TH1F("hIpErLep2_7", "IpEr of Lep2 after selection step 7",  1000 , 0. , 10. );
   hIpErLep2_7->SetXTitle("IpEr");

   TH1F * hPtLep3_7 = new TH1F("hPtLep3_7", "Pt of Lep3 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep3_7->SetXTitle("pt_Lep3  (GeV)");
   TH1F * hYLep3_7 = new TH1F("hEtaLep3_7", "Y of Lep3 after selection step 7", 500 , -5. , 5. );
   hYLep3_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep3_7 = new TH1F("hIsoLep3_7", "Isolation of Lep3 after selection step 7", 2000 , -10. , 10. );
   hIsoLep3_7->SetXTitle("Iso");
   TH1F * hSipLep3_7 = new TH1F("hSipLep3_7", "Sip of Lep3 after selection step 7",  1000 , -20. , 40. );
   hSipLep3_7->SetXTitle("Sip");
   TH1F * hIpLep3_7 = new TH1F("hIpLep3_7", "Ip of Lep3 after selection step 7",  1000 , -20. , 40. );
   hIpLep3_7->SetXTitle("Ip");
   TH1F * hIpErLep3_7 = new TH1F("hIpErLep3_7", "IpEr of Lep3 after selection step 7",  1000 , 0. , 10. );
   hIpErLep3_7->SetXTitle("IpEr");

   TH1F * hPtLep4_7 = new TH1F("hPtLep4_7", "Pt of Lep4 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep4_7->SetXTitle("pt_Lep4  (GeV)");
   TH1F * hYLep4_7 = new TH1F("hEtaLep4_7", "Y of Lep4 after selection step 7", 500 , -5. , 5. );
   hYLep4_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep4_7 = new TH1F("hIsoLep4_7", "Isolation of Lep4 after selection step 7", 2000 , -10. , 10. );
   hIsoLep4_7->SetXTitle("Iso");
   TH1F * hSipLep4_7 = new TH1F("hSipLep4_7", "Sip of Lep4 after selection step 7",  1000 , -20. , 40. );
   hSipLep4_7->SetXTitle("Sip");
   TH1F * hIpLep4_7 = new TH1F("hIpLep4_7", "Ip of Lep4 after selection step 7",  1000 , -20. , 40. );
   hIpLep4_7->SetXTitle("Ip");
   TH1F * hIpErLep4_7 = new TH1F("hIpErLep4_7", "IpEr of Lep4 after selection step 7",  1000 , 0. , 10. );
   hIpErLep4_7->SetXTitle("IpEr");

   TH1F * hIso_7 = new TH1F("hIso_7", "Isolation maxima after selection step 7", 2000 , -10. , 10. );
   hIso_7->SetXTitle("Iso");
   TH1F * hSip_7 = new TH1F("hSip_7", "Sip maxima after selection step 7",  1000 , -20. , 40. );
   hSip_7->SetXTitle("Sip");
   TH1F * hIp_7 = new TH1F("hIp_7", "Ip maxima after selection step 7",  1000 , -20. , 40. );
   hIp_7->SetXTitle("Ip");
   TH1F * hPFMET_7 = new TH1F("hPFMET_7", "PF MET after selection step 7", 1000 , 0., 1000.);
   hPFMET_7->SetXTitle("PF MET (GeV)"); 

   TH1F * hDPHI_7 = new TH1F("DPHI_7", "polar angle between 4l and E_{T,miss}", 1000, 0., 5. );
   hDPHI_7->SetXTitle("#DELTA#phi(4l,E_{T,miss})");


   //step 8
   TH1F * hMZ1_61 = new TH1F("hMZ1_61", "Mass of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hMZ1_61->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_61 = new TH1F("hPtZ1_61", "Pt of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ1_61->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_61 = new TH1F("hYZ1_61", "Y of Z1 after selection step 7", 500 , -5. , 5.);
   hYZ1_61->SetXTitle("Y_Z1");   

   TH1F * hMZ2_61 = new TH1F("hMZ2_61", "Mass of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hMZ2_61->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_61 = new TH1F("hPtZ2_61", "Pt of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ2_61->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_61 = new TH1F("hYZ2_61", "Y of Z2 after selection step 7", 500 , -5. , 5.);
   hYZ2_61->SetXTitle("Y_Z2");

   TH1F * hMZ1_62 = new TH1F("hMZ1_62", "Mass of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hMZ1_62->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_62 = new TH1F("hPtZ1_62", "Pt of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ1_62->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_62 = new TH1F("hYZ1_62", "Y of Z1 after selection step 7", 500 , -5. , 5.);
   hYZ1_62->SetXTitle("Y_Z1");

   TH1F * hMZ2_62 = new TH1F("hMZ2_62", "Mass of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hMZ2_62->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_62 = new TH1F("hPtZ2_62", "Pt of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ2_62->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_62 = new TH1F("hYZ2_62", "Y of Z2 after selection step 7", 500 , -5. , 5.);
   hYZ2_62->SetXTitle("Y_Z2");


   TH1F * hPtLep1_8 = new TH1F("hPtLep1_8", "Pt of Lep1 after selection step 8", 200 , -0.5 , 199.5 );
   hPtLep1_8->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_8 = new TH1F("hEtaLep1_8", "Y of Lep1 after selection step 8", 500 , -5. , 5. );
   hYLep1_8->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_8 = new TH1F("hIsoLep1_8", "Isolation of Lep1 after selection step 8", 2000 , -10. , 10. );
   hIsoLep1_8->SetXTitle("Iso");
   TH1F * hSipLep1_8 = new TH1F("hSipLep1_8", "Sip of Lep1 after selection step 8",  1000 , -20. , 40. );
   hSipLep1_8->SetXTitle("Sip");
   TH1F * hIpLep1_8 = new TH1F("hIpLep1_8", "Ip of Lep1 after selection step 8",  1000 , -20. , 40. );
   hIpLep1_8->SetXTitle("Ip");
   TH1F * hIpErLep1_8 = new TH1F("hIpErLep1_8", "IpEr of Lep1 after selection step 8",  1000 , 0. , 10. );
   hIpErLep1_8->SetXTitle("IpEr");

   
   TH1F * hPtLep2_8 = new TH1F("hPtLep2_8", "Pt of Lep2 after selection step 8", 200 , -0.5 , 199.5 );
   hPtLep2_8->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_8 = new TH1F("hEtaLep2_8", "Y of Lep2 after selection step 8", 500 , -5. , 5. );
   hYLep2_8->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_8 = new TH1F("hIsoLep2_8", "Isolation of Lep2 after selection step 8", 2000 , -10. , 10. );
   hIsoLep2_8->SetXTitle("Iso");
   TH1F * hSipLep2_8 = new TH1F("hSipLep2_8", "Sip of Lep2 after selection step 8",  1000 , -20. , 40. );
   hSipLep2_8->SetXTitle("Sip");
   TH1F * hIpLep2_8 = new TH1F("hIpLep2_8", "Ip of Lep2 after selection step 8",  1000 , -20. , 40. );
   hIpLep2_8->SetXTitle("Ip");
   TH1F * hIpErLep2_8 = new TH1F("hIpErLep2_8", "IpEr of Lep2 after selection step 8",  1000 , 0. , 10. );
   hIpErLep2_8->SetXTitle("IpEr");

   TH1F * hPtLep3_8 = new TH1F("hPtLep3_8", "Pt of Lep3 after selection step 8", 200 , -0.5 , 199.5 );
   hPtLep3_8->SetXTitle("pt_Lep3  (GeV)");
   TH1F * hYLep3_8 = new TH1F("hEtaLep3_8", "Y of Lep3 after selection step 8", 500 , -5. , 5. );
   hYLep3_8->SetXTitle("Y of Lep2");
   TH1F * hIsoLep3_8 = new TH1F("hIsoLep3_8", "Isolation of Lep3 after selection step 8", 2000 , -10. , 10. );
   hIsoLep3_8->SetXTitle("Iso");
   TH1F * hSipLep3_8 = new TH1F("hSipLep3_8", "Sip of Lep3 after selection step 8",  1000 , -20. , 40. );
   hSipLep3_8->SetXTitle("Sip");
   TH1F * hIpLep3_8 = new TH1F("hIpLep3_8", "Ip of Lep3 after selection step 8",  1000 , -20. , 40. );
   hIpLep3_8->SetXTitle("Ip");
   TH1F * hIpErLep3_8 = new TH1F("hIpErLep3_8", "IpEr of Lep3 after selection step 8",  1000 , 0. , 10. );
   hIpErLep3_8->SetXTitle("IpEr");

   TH1F * hPtLep4_8 = new TH1F("hPtLep4_8", "Pt of Lep4 after selection step 8", 200 , -0.5 , 199.5 );
   hPtLep4_8->SetXTitle("pt_Lep4  (GeV)");
   TH1F * hYLep4_8 = new TH1F("hEtaLep4_8", "Y of Lep4 after selection step 8", 500 , -5. , 5. );
   hYLep4_8->SetXTitle("Y of Lep2");
   TH1F * hIsoLep4_8 = new TH1F("hIsoLep4_8", "Isolation of Lep4 after selection step 8", 2000 , -10. , 10. );
   hIsoLep4_8->SetXTitle("Iso");
   TH1F * hSipLep4_8 = new TH1F("hSipLep4_8", "Sip of Lep4 after selection step 8",  1000 , -20. , 40. );
   hSipLep4_8->SetXTitle("Sip");
   TH1F * hIpLep4_8 = new TH1F("hIpLep4_8", "Ip of Lep4 after selection step 8",  1000 , -20. , 40. );
   hIpLep4_8->SetXTitle("Ip");
   TH1F * hIpErLep4_8 = new TH1F("hIpErLep4_8", "IpEr of Lep4 after selection step 8",  1000 , 0. , 10. );
   hIpErLep4_8->SetXTitle("IpEr");
   

   TH1F * hIso_8 = new TH1F("hIso_8", "Isolation maxima after selection step 8", 2000 , -10. , 10. );
   hIso_8->SetXTitle("Iso");
   TH1F * hSip_8 = new TH1F("hSip_8", "Sip maxima after selection step 8",  1000 , -20. , 40. );
   hSip_8->SetXTitle("Sip");
   TH1F * hIp_8 = new TH1F("hIp_8", "Ip maxima after selection step 8",  1000 , -20. , 40. );
   hIp_8->SetXTitle("Ip");
   
   TH1F * hDjj_8 = new TH1F("hDjj_8", "Delta jets vbf selection step 8", 200, -19.5, 19.5 );
   hDjj_8->SetXTitle("Delta jets");
   TH1F * hMjj_8 = new TH1F("hMjj_8", "Mass jets vbf selection step 8", 200, -0.5, 499.5 );
   hMjj_8->SetXTitle("Mass jets");
   TH1F * hVD_8 = new TH1F("hVD_8", "Discriminant vbf selection step 8", 200, -0.5, 9.5 );
   hMjj_8->SetXTitle("Discriminant");
   TH2F * hMjj_Djj_8 = new TH2F("hMjj_Djj_8", "Mjj vs Djj selection step 8", 200, -0.5,499.5,200,-19.5, 19.5 );
   hMjj_Djj_8->SetXTitle("Mjj_vs_Djj");

   
   TH1F * hPFMET_8 = new TH1F("hPFMET_8", "PF MET after selection step 8", 1000 , 0., 1000.);
   hPFMET_8->SetXTitle("PF MET");

   TH1F * hPFMET_8_pdf_up = new TH1F("hPFMET_8_pdf_up", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_pdf_dow = new TH1F("hPFMET_8_pdf_dow", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_pu_up = new TH1F("hPFMET_8_pu_up", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_pu_dow = new TH1F("hPFMET_8_pu_dow", "PF MET after selection step 8", 1000 , 0., 1000.);

   TH1F * hPFMET_8_mc_1 = new TH1F("hPFMET_8_mc_1", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_mc_2 = new TH1F("hPFMET_8_mc_2", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_mc_3 = new TH1F("hPFMET_8_mc_3", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_mc_4 = new TH1F("hPFMET_8_mc_4", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_mc_5 = new TH1F("hPFMET_8_mc_5", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_mc_6 = new TH1F("hPFMET_8_mc_6", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_mc_7 = new TH1F("hPFMET_8_mc_7", "PF MET after selection step 8", 1000 , 0., 1000.);
   TH1F * hPFMET_8_mc_8 = new TH1F("hPFMET_8_mc_8", "PF MET after selection step 8", 1000 , 0., 1000.);
   

   TH1D * hNgood_8 = new TH1D("hNgood", "Number of good leptons", 10, -0.5, 9.5);
   hNgood_8->SetXTitle("# good leptons");

   TH1D * hNbjets_8 = new TH1D("hNbjets", "Number of b jets", 10, -0.5, 9.5);
   hNbjets_8->SetXTitle("# b-jets");

   TH1D * hNbjets_8_pdf_up = new TH1D("hNbjets_pdf_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_pdf_dow = new TH1D("hNbjets_pdf_dow", "Number of b jets", 10, -0.5, 9.5);

   TH1D * hNbjets_8_pu_up = new TH1D("hNbjets_pu_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_pu_dow = new TH1D("hNbjets_pu_dow", "Number of b jets", 10, -0.5, 9.5);

   TH1D * hNbjets_8_mc_1 = new TH1D("hNbjets_mc_1", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_mc_2 = new TH1D("hNbjets_mc_2", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_mc_3 = new TH1D("hNbjets_mc_3", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_mc_4 = new TH1D("hNbjets_mc_4", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_mc_5 = new TH1D("hNbjets_mc_5", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_mc_6 = new TH1D("hNbjets_mc_6", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_mc_7 = new TH1D("hNbjets_mc_7", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_mc_8 = new TH1D("hNbjets_mc_8", "Number of b jets", 10, -0.5, 9.5);


   TH1D * h_Nbjets_8 = new TH1D("h_Nbjets_8", "Number of b jets", 10, -0.5, 9.5);
   h_Nbjets_8->SetXTitle("# b-jets");

   TH1D * hNjets_8 = new TH1D("hNjets_8", "Number of jets passing VBF", 10, -0.5, 9.5);
   hNjets_8->SetXTitle("# n-jets");

   TH1D * hNjets_8_pdf_up = new TH1D("hNjets_8_pdf_up", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_pdf_dow = new TH1D("hNjets_8_pdf_dow", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_pu_up = new TH1D("hNjets_8_pu_up", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_pu_dow = new TH1D("hNjets_8_pu_dow", "Number of jets passing VBF", 10, -0.5, 9.5);

   TH1D * hNjets_8_mc_1 = new TH1D("hNjets_8_mc_1", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_mc_2 = new TH1D("hNjets_8_mc_2", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_mc_3 = new TH1D("hNjets_8_mc_3", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_mc_4 = new TH1D("hNjets_8_mc_4", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_mc_5 = new TH1D("hNjets_8_mc_5", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_mc_6 = new TH1D("hNjets_8_mc_6", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_mc_7 = new TH1D("hNjets_8_mc_7", "Number of jets passing VBF", 10, -0.5, 9.5);
   TH1D * hNjets_8_mc_8 = new TH1D("hNjets_8_mc_8", "Number of jets passing VBF", 10, -0.5, 9.5);


   TH1F * hPtJet_7 = new TH1F("hPtJet_7", "Pt of (no ID)jet after selection step 5", 300 ,  0 , 600 );
   hPtJet_7->SetXTitle("pt_jet  (GeV)");

   TH1F * hYJet_7 = new TH1F("hEtaJet_7", "Y of (no ID)jet after selection step 5", 500 , -5. , 5. );
   hYJet_7->SetXTitle("Y of jet7");

   TH1F * hPtJet_8 = new TH1F("hPtJet_8", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   hPtJet_8->SetXTitle("pt_jet  (GeV)");

   TH1F * hPtJet_8_up = new TH1F("hPtJet_8_up", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_dow = new TH1F("hPtJet_8_dow", "Pt of jet after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtJet_8_jer_up = new TH1F("hPtJet_8_jer_up", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_jer_dow = new TH1F("hPtJet_8_jer_dow", "Pt of jet after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtJet_8_pdf_up = new TH1F("hPtJet_8_pdf_up", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_pdf_dow = new TH1F("hPtJet_8_pdf_dow", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_pu_up = new TH1F("hPtJet_8_pu_up", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_pu_dow = new TH1F("hPtJet_8_pu_dow", "Pt of jet after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtJet_8_mc_1 = new TH1F("hPtJet_8_mc_1", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_mc_2 = new TH1F("hPtJet_8_mc_2", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_mc_3 = new TH1F("hPtJet_8_mc_3", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_mc_4 = new TH1F("hPtJet_8_mc_4", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_mc_5 = new TH1F("hPtJet_8_mc_5", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_mc_6 = new TH1F("hPtJet_8_mc_6", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_mc_7 = new TH1F("hPtJet_8_mc_7", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtJet_8_mc_8 = new TH1F("hPtJet_8_mc_8", "Pt of jet after selection step 5", 300 ,  0 , 600 );


   TH1F * hYJet_8 = new TH1F("hEtaJet_8", "Y of jet after selection step 5", 500 , -5. , 5. );
   hYJet_8->SetXTitle("Y of Jet8");

   TH1F * hYJet_8_pdf_up = new TH1F("hEtaJet_8_pdf_up", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_pdf_dow = new TH1F("hEtaJet_8_pdf_dow", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_pu_up = new TH1F("hEtaJet_8_pu_up", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_pu_dow = new TH1F("hEtaJet_8_pu_dow", "Y of jet after selection step 5", 500 , -5. , 5. );

   TH1F * hYJet_8_mc_1 = new TH1F("hEtaJet_8_mc_1", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_mc_2 = new TH1F("hEtaJet_8_mc_2", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_mc_3 = new TH1F("hEtaJet_8_mc_3", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_mc_4 = new TH1F("hEtaJet_8_mc_4", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_mc_5 = new TH1F("hEtaJet_8_mc_5", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_mc_6 = new TH1F("hEtaJet_8_mc_6", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_mc_7 = new TH1F("hEtaJet_8_mc_7", "Y of jet after selection step 5", 500 , -5. , 5. );
   TH1F * hYJet_8_mc_8 = new TH1F("hEtaJet_8_mc_8", "Y of jet after selection step 5", 500 , -5. , 5. );

   
   // Step 9 with PFMET cut
   
   TH1F * hPFMET_9 = new TH1F("hPFMET_9", "PF MET after selection step 9", 1000 , 0., 1000.);
   hPFMET_9->SetXTitle("PF MET (GeV)");   
   
   TH1F * hM4l_T_9 = new TH1F("hM4l_T_9", "Transverse Mass of four leptons after full selection + MET", 1200, 4.5, 1204.5 );
   hM4l_T_9->SetXTitle("m_{T} + PF MET (GeV)");
   
   TH1F * hDPHI_4l_9 = new TH1F("DPHI_4l_9", "polar angle between 4l and E_{T,miss}", 3150, 0., 3.15 );
   hDPHI_4l_9->SetXTitle("#Delta#phi(4l,E_{T,miss})");

   TH1F * hCVSBTag_PT_9 = new TH1F("hCVSBTag_PT_9", "...", 900, 0., 300.);
   TH1F * hCVSBTag_DISCR_9 = new TH1F("hCVSBTag_DISCR_9", "...", 500, 0., 1.);

   // Step 10 with signal region cuts
   TH1F * hPFMET_10 = new TH1F("hPFMET_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_10->SetXTitle("PF MET (GeV)");  
   
   TH1F * hM4l_T_10 = new TH1F("hM4l_T_10", "Transverse Mass of four leptons + MET after step 10", 1200, 4.5, 1204.5 );
   hM4l_T_10->SetXTitle("m_{T} + PF MET (GeV)");
   
   TH1F * hDPHI_4l_10 = new TH1F("DPHI_4l_10", "polar angle between 4l and E_{T,miss} after step 10", 3150, 0., 3.15 );
   hDPHI_4l_10->SetXTitle("#DELTA#phi(4l,E_{T,miss})");


   
   TH1F * hN_loose_mu = new TH1F("hN_loose_mu", "N_loose_mu", 30 , 0. , 30. );
   hN_loose_mu->SetXTitle("N_loose_mu");
   TH1F * hN_loose_e = new TH1F("hN_loose_e", "N_loose_e", 30 , 0. , 30. );
   hN_loose_e->SetXTitle("N_loose_e");

   TH1F * hIso_loose_mu = new TH1F("hIso_loose_mu", "Isolation maxima after loose selection ", 2000 , -10. , 10. );
   hIso_loose_mu->SetXTitle("Iso");
   TH1F * hSip_loose_mu = new TH1F("hSip_loose_mu", "Sip maxima after loose selection ",  1000 , -20. , 40. );
   hSip_loose_mu->SetXTitle("Sip");
   TH1F * hIp_loose_mu = new TH1F("hIp_loose_mu", "Ip maxima after loose selection ",  1000 , -20. , 40. );
   hIp_loose_mu->SetXTitle("Ip");

   TH1F * hIso_loose_e = new TH1F("hIso_loose_e", "Isolation maxima after loose selection ", 2000 , -10. , 10. );
   hIso_loose_e->SetXTitle("Iso");
   TH1F * hSip_loose_e = new TH1F("hSip_loose_e", "Sip maxima after loose selection ",  1000 , -20. , 40. );
   hSip_loose_e->SetXTitle("Sip");
   TH1F * hIp_loose_e = new TH1F("hIp_loose_e", "Ip maxima after loose selection ",  1000 , -20. , 40. );
   hIp_loose_e->SetXTitle("Ip");


   TH1F * hN_good_lep = new TH1F("hN_good_lep", "N_good_lep", 30 , 0. , 30. );
   hN_good_lep->SetXTitle("N_good_lep");

   TH1F * hN_good_mu = new TH1F("hN_good_mu", "N_good_mu", 30 , 0. , 30. );
   hN_good_mu->SetXTitle("N_good_mu");
   TH1F * hN_good_ele = new TH1F("hN_good_ele", "N_good_ele", 30 , 0. , 30. );
   hN_good_ele->SetXTitle("N_good_ele");
   TH1F * hN_good_phot = new TH1F("hN_good_phot", "N_good_phot", 30 , 0. , 30. );
   hN_good_phot->SetXTitle("N_good_phot");

   TH1D * hfit = new TH1D("hfit", "Number of b jets", 10, -0.5, 9.5);
   hfit->SetXTitle("# b-jets");

   TH1D * hfit_pdf_up = new TH1D("hfit_pdf_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_pdf_dow = new TH1D("hfit_pdf_dow", "Number of b jets", 10, -0.5, 9.5);

   TH1D * hfit_pu_up = new TH1D("hfit_pu_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_pu_dow = new TH1D("hfit_pu_dow", "Number of b jets", 10, -0.5, 9.5);

   TH1D * hfit_mc_1 = new TH1D("hfit_mc_1", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_mc_2 = new TH1D("hfit_mc_2", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_mc_3 = new TH1D("hfit_mc_3", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_mc_4 = new TH1D("hfit_mc_4", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_mc_5 = new TH1D("hfit_mc_5", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_mc_6 = new TH1D("hfit_mc_6", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_mc_7 = new TH1D("hfit_mc_7", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hfit_mc_8 = new TH1D("hfit_mc_8", "Number of b jets", 10, -0.5, 9.5);


   //PFJET Plots
/*    TH1I * hN_PFJET_6 = new TH1I("hN_PFJET_6", "Number of PFJets after step 6", 102, 0, 102);
   hN_PFJET_6->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_6 = new TH1F("hChg_PFJET_6", "Charge of PFJets after step 6", 50, -10., 10.);
   hChg_PFJET_6->SetXTitle("Charge");
   TH1F * hEt_PFJET_6 = new TH1F("hEt_PFJET_6", "Et of PFJets after step 6", 201, -0.5, 400.5 );
   hEt_PFJET_6->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_6 = new TH1F("hPt_PFJET_6", "Pt of PFJets after step 6", 201, -0.5, 400.5 );
   hPt_PFJET_6->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_6 = new TH1F("hEta_PFJET_6", "Eta of PFJets after step 6", 50, -10., 10.);
   hEta_PFJET_6->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_6 = new TH1F("hPhi_PFJET_6", "Phi of PFJets after step 6", 50, -10., 10.);
   hPhi_PFJET_6->SetXTitle("phi_PFJET");

   TH1I * hN_PFJET_8 = new TH1I("hN_PFJET_8", "Number of PFJets after step 8", 102, 0, 102);
   hN_PFJET_8->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_8 = new TH1F("hChg_PFJET_8", "Charge of PFJets after step 8", 50, -10., 10.);
   hChg_PFJET_8->SetXTitle("Charge");
   TH1F * hEt_PFJET_8 = new TH1F("hEt_PFJET_8", "Et of PFJets after step 8", 201, -0.5, 400.5 );
   hEt_PFJET_8->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_8 = new TH1F("hPt_PFJET_8", "Pt of PFJets after step 8", 201, -0.5, 400.5 );
   hPt_PFJET_8->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_8 = new TH1F("hEta_PFJET_8", "Eta of PFJets after step 8", 50, -10., 10.);
   hEta_PFJET_8->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_8 = new TH1F("hPhi_PFJET_8", "Phi of PFJets after step 8", 50, -10., 10.);
   hPhi_PFJET_8->SetXTitle("phi_PFJET");

   TH1I * hN_PFJET_VBF = new TH1I("hN_PFJET_VBF", "Number of PFJets after VBF selection", 102, 0, 102);
   hN_PFJET_VBF->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_VBF = new TH1F("hChg_PFJET_VBF", "Charge of PFJets after VBF selection", 50, -10., 10.);
   hChg_PFJET_VBF->SetXTitle("Charge");
   TH1F * hEt_PFJET_VBF = new TH1F("hEt_PFJET_VBF", "Et of PFJets after step VBF selection", 201, -0.5, 400.5 );
   hEt_PFJET_VBF->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_VBF = new TH1F("hPt_PFJET_VBF", "Pt of PFJets after step VBF selection", 201, -0.5, 400.5 );
   hPt_PFJET_VBF->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_VBF = new TH1F("hEta_PFJET_VBF", "Eta of PFJets after step VBF selection", 50, -10., 10.);
   hEta_PFJET_VBF->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_VBF = new TH1F("hPhi_PFJET_VBF", "Phi of PFJets after step VBF selection", 50, -10., 10.);
   hPhi_PFJET_VBF->SetXTitle("phi_PFJET");
*/   

//Bottom plot
    TH1I * hN_PFJET_6 = new TH1I("hN_PFJET_6", "Number of PFJets after step 6", 102, 0, 102);
   hN_PFJET_6->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_6 = new TH1F("hChg_PFJET_6", "Charge of PFJets after step 6", 50, -10., 10.);
   hChg_PFJET_6->SetXTitle("Charge");
   TH1F * hEt_PFJET_6 = new TH1F("hEt_PFJET_6", "Et of PFJets after step 6", 201, -0.5, 400.5 );
   hEt_PFJET_6->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_6 = new TH1F("hPt_PFJET_6", "Pt of PFJets after step 6", 201, -0.5, 400.5 );
   hPt_PFJET_6->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_6 = new TH1F("hEta_PFJET_6", "Eta of PFJets after step 6", 50, -10., 10.);
   hEta_PFJET_6->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_6 = new TH1F("hPhi_PFJET_6", "Phi of PFJets after step 6", 50, -10., 10.);
   hPhi_PFJET_6->SetXTitle("phi_PFJET");

   TH1F * hPtBot_8 = new TH1F("hPtBot_8", "Pt of bot after selection step 5", 300 ,  0 , 600 );
   hPtBot_8->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_up = new TH1F("hPtBot_8_up", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_dow = new TH1F("hPtBot_8_dow", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_pdf_up = new TH1F("hPtBot_8_pdf_up", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_pdf_dow = new TH1F("hPtBot_8_pdf_dow", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_pu_up = new TH1F("hPtBot_8_pu_up", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_pu_dow = new TH1F("hPtBot_8_pu_dow", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_jer_up = new TH1F("hPtBot_8_jer_up", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_jer_dow = new TH1F("hPtBot_8_jer_dow", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );


   TH1F * hPtBot_8_mc_1 = new TH1F("hPtBot_8_mc_1", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_mc_2 = new TH1F("hPtBot_8_mc_2", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_mc_3 = new TH1F("hPtBot_8_mc_3", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_mc_4 = new TH1F("hPtBot_8_mc_4", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_mc_5 = new TH1F("hPtBot_8_mc_5", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_mc_6 = new TH1F("hPtBot_8_mc_6", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_mc_7 = new TH1F("hPtBot_8_mc_7", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );
   TH1F * hPtBot_8_mc_8 = new TH1F("hPtBot_8_mc_8", "Pt of bot after selection step 5(JEC)", 300 ,  0 , 600 );


   TH1F * hYBot_8 = new TH1F("hEtaBot_8", "Y of bot after selection step 5", 500 , -5. , 5. );
   hYBot_8->SetXTitle("Y of Bot8");

   TH1F * hPtBot_7 = new TH1F("hPtBot_7", "Pt of bot after selection step 5", 300 ,  0 , 600 );
   hPtBot_7->SetXTitle("pt_bot  (GeV)");


   TH1F * Mbb_6 = new TH1F("Mbb_6","invariant mass of bottom pair after step 6",50,20,420);
   Mbb_6->SetXTitle("M_{bb} (GeV)");
   TH1F * Mbb_6_up = new TH1F("Mbb_6_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_dow = new TH1F("Mbb_6_dow","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_jer_up = new TH1F("Mbb_6_jer_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_jer_dow = new TH1F("Mbb_6_jer_dow","invariant mass of bottom pair after step 6",50,20,420);

   TH1F * Mbb_6_pdf_up = new TH1F("Mbb_6_pdf_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_pdf_dow = new TH1F("Mbb_6_pdf_dow","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_pu_up = new TH1F("Mbb_6_pu_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_pu_dow = new TH1F("Mbb_6_pu_dow","invariant mass of bottom pair after step 6",50,20,420);

   TH1F * Mbb_6_mc_1 = new TH1F("Mbb_6_mc_1","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbb_6_mc_2 = new TH1F("Mbb_6_mc_2","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbb_6_mc_3 = new TH1F("Mbb_6_mc_3","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbb_6_mc_4 = new TH1F("Mbb_6_mc_4","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbb_6_mc_5 = new TH1F("Mbb_6_mc_5","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbb_6_mc_6 = new TH1F("Mbb_6_mc_6","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbb_6_mc_7 = new TH1F("Mbb_6_mc_7","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbb_6_mc_8 = new TH1F("Mbb_6_mc_8","invariant mass of jet pair after step 6",50,20,420);

   TH1F * Mbj_6 = new TH1F("Mbj_6","invariant mass of bottom pair after step 6",50,20,420);
   Mbj_6->SetXTitle("M_{bj} (GeV)");
   TH1F * Mbj_6_up = new TH1F("Mbj_6_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbj_6_dow = new TH1F("Mbj_6_dow","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbj_6_jer_up = new TH1F("Mbj_6_jer_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbj_6_jer_dow = new TH1F("Mbj_6_jer_dow","invariant mass of bottom pair after step 6",50,20,420);

   TH1F * Mbj_6_pdf_up = new TH1F("Mbj_6_pdf_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbj_6_pdf_dow = new TH1F("Mbj_6_pdf_dow","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbj_6_pu_up = new TH1F("Mbj_6_pu_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbj_6_pu_dow = new TH1F("Mbj_6_pu_dow","invariant mass of bottom pair after step 6",50,20,420);

   TH1F * Mbj_6_mc_1 = new TH1F("Mbj_6_mc_1","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbj_6_mc_2 = new TH1F("Mbj_6_mc_2","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbj_6_mc_3 = new TH1F("Mbj_6_mc_3","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbj_6_mc_4 = new TH1F("Mbj_6_mc_4","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbj_6_mc_5 = new TH1F("Mbj_6_mc_5","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbj_6_mc_6 = new TH1F("Mbj_6_mc_6","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbj_6_mc_7 = new TH1F("Mbj_6_mc_7","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mbj_6_mc_8 = new TH1F("Mbj_6_mc_8","invariant mass of jet pair after step 6",50,20,420);

   TH1F * Mjj_7 = new TH1F("Mjj_7","invariant mass of jet pair after step 6",50,20,420);
   Mjj_7->SetXTitle("M_{jj} (GeV)");

   TH1F * Mjj_7_jer_up = new TH1F("Mjj_7_jer_up","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_jer_dow = new TH1F("Mjj_7_jer_dow","invariant mass of jet pair after step 6",50,20,420);

   TH1F * Mjj_7_pdf_up = new TH1F("Mjj_7_pdf_up","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_pdf_dow = new TH1F("Mjj_7_pdf_dow","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_pu_up = new TH1F("Mjj_7_pu_up","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_pu_dow = new TH1F("Mjj_7_pu_dow","invariant mass of jet pair after step 6",50,20,420);

   TH1F * Mjj_7_up = new TH1F("Mjj_7_up","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_dow = new TH1F("Mjj_7_dow","invariant mass of jet pair after step 6",50,20,420);

   TH1F * Mjj_7_mc_1 = new TH1F("Mjj_7_mc_1","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_mc_2 = new TH1F("Mjj_7_mc_2","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_mc_3 = new TH1F("Mjj_7_mc_3","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_mc_4 = new TH1F("Mjj_7_mc_4","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_mc_5 = new TH1F("Mjj_7_mc_5","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_mc_6 = new TH1F("Mjj_7_mc_6","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_mc_7 = new TH1F("Mjj_7_mc_7","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_7_mc_8 = new TH1F("Mjj_7_mc_8","invariant mass of jet pair after step 6",50,20,420);


   TH1F * ptbb_6 = new TH1F("ptbb_6","pt of bottom pair after step 6",50,20,420);
   ptbb_6->SetXTitle("pt_{bb} (GeV)");
   TH1F * ptbb_6_up = new TH1F("ptbb_6_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_dow = new TH1F("ptbb_6_dow","pt of bottom pair after step 6",50,20,420);

   TH1F * ptbb_6_pdf_up = new TH1F("ptbb_6_pdf_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_pdf_dow = new TH1F("ptbb_6_pdf_dow","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_pu_up = new TH1F("ptbb_6_pu_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_pu_dow = new TH1F("ptbb_6_pu_dow","pt of bottom pair after step 6",50,20,420);

   TH1F * ptbb_6_jer_up = new TH1F("ptbb_6_jer_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_jer_dow = new TH1F("ptbb_6_jer_dow","pt of bottom pair after step 6",50,20,420);

   TH1F * ptbb_6_mc_1 = new TH1F("ptbb_6_mc_1","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_mc_2 = new TH1F("ptbb_6_mc_2","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_mc_3 = new TH1F("ptbb_6_mc_3","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_mc_4 = new TH1F("ptbb_6_mc_4","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_mc_5 = new TH1F("ptbb_6_mc_5","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_mc_6 = new TH1F("ptbb_6_mc_6","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_mc_7 = new TH1F("ptbb_6_mc_7","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_mc_8 = new TH1F("ptbb_6_mc_8","pt of bottom pair after step 6",50,20,420);

   TH1F * bdiscr_5_lead = new TH1F("bdiscr_5_lead","b-tag discr of leading jet after step 5",20,0,1);
   bdiscr_5_lead->SetXTitle("CSV");
   TH1F * bdiscr_5_sub = new TH1F("bdiscr_5_sub","b-tag discr of sub leading jet after step 5",20,0,1);
   bdiscr_5_sub->SetXTitle("CSV");

   TH1F * Mjj_6 = new TH1F("Mjj_6","invariant mass of jet pair after step 6",50,20,420);
   Mjj_6->SetXTitle("M_{jj} (GeV)");  

   TH1F * Mjj_6_jer_up = new TH1F("Mjj_6_jer_up","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_jer_dow = new TH1F("Mjj_6_jer_dow","invariant mass of jet pair after step 6",50,20,420);

   TH1F * Mjj_6_pdf_up = new TH1F("Mjj_6_pdf_up","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_pdf_dow = new TH1F("Mjj_6_pdf_dow","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_pu_up = new TH1F("Mjj_6_pu_up","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_pu_dow = new TH1F("Mjj_6_pu_dow","invariant mass of jet pair after step 6",50,20,420);

   TH1F * Mjj_6_up = new TH1F("Mjj_6_up","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_dow = new TH1F("Mjj_6_dow","invariant mass of jet pair after step 6",50,20,420);

   TH1F * Mjj_6_mc_1 = new TH1F("Mjj_6_mc_1","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_mc_2 = new TH1F("Mjj_6_mc_2","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_mc_3 = new TH1F("Mjj_6_mc_3","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_mc_4 = new TH1F("Mjj_6_mc_4","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_mc_5 = new TH1F("Mjj_6_mc_5","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_mc_6 = new TH1F("Mjj_6_mc_6","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_mc_7 = new TH1F("Mjj_6_mc_7","invariant mass of jet pair after step 6",50,20,420);
   TH1F * Mjj_6_mc_8 = new TH1F("Mjj_6_mc_8","invariant mass of jet pair after step 6",50,20,420);

   // end book histo ***

   TH1F * hPtJet_8_b = new TH1F("hPtJet_8_b", "Pt of jet_b after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtBot_8_b = new TH1F("hPtBot_8_b", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtBot_8_b_dow = new TH1F("hPtBot_8_b_dow", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtBot_8_b_up = new TH1F("hPtBot_8_b_up", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtJet_8_c = new TH1F("hPtJet_8_c", "Pt of jet_c after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtBot_8_c = new TH1F("hPtBot_8_c", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtBot_8_c_up = new TH1F("hPtBot_8_c_up", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtBot_8_c_dow = new TH1F("hPtBot_8_c_dow", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtJet_8_l = new TH1F("hPtJet_8_l", "Pt of jet_l after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtBot_8_l = new TH1F("hPtBot_8_l", "Pt of bot_l after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtBot_8_l_up = new TH1F("hPtBot_8_l_up", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtBot_8_l_dow = new TH1F("hPtBot_8_l_dow", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtJet_8_o = new TH1F("hPtJet_8_o", "Pt of jet_o after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtBot_8_o = new TH1F("hPtBot_8_o", "Pt of bot_o after selection step 5", 300 ,  0 , 600 );

   TH1F * hPtBot_8_o_up = new TH1F("hPtBot_8_o_up", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   TH1F * hPtBot_8_o_dow = new TH1F("hPtBot_8_o_dow", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );

   TH1F * Mbb_6_hh = new TH1F("Mbb_6_hh","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hl = new TH1F("Mbb_6_hl","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_ll = new TH1F("Mbb_6_ll","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hh_up = new TH1F("Mbb_6_hh_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hl_up = new TH1F("Mbb_6_hl_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_ll_up = new TH1F("Mbb_6_ll_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hh_dow = new TH1F("Mbb_6_hh_dow","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hl_dow = new TH1F("Mbb_6_hl_dow","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_ll_dow = new TH1F("Mbb_6_ll_dow","invariant mass of bottom pair after step 6",50,20,420);

   TH1F * ptbb_6_hh = new TH1F("ptbb_6_hh","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hl = new TH1F("ptbb_6_hl","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_ll = new TH1F("ptbb_6_ll","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hh_up = new TH1F("ptbb_6_hh_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hl_up = new TH1F("ptbb_6_hl_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_ll_up = new TH1F("ptbb_6_ll_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hh_dow = new TH1F("ptbb_6_hh_dow","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hl_dow = new TH1F("ptbb_6_hl_dow","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_ll_dow = new TH1F("ptbb_6_ll_dow","pt of bottom pair after step 6",50,20,420);
   TH1D * hNbjets_8_h = new TH1D("hNbjets_h", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_h_up = new TH1D("hNbjets_h_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_h_dow = new TH1D("hNbjets_h_dow", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_l = new TH1D("hNbjets_l", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_l_up = new TH1D("hNbjets_l_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_l_dow = new TH1D("hNbjets_l_dow", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_btag_up = new TH1D("hNbjets_btag_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_btag_dow = new TH1D("hNbjets_btag_dow", "Number of b jets", 10, -0.5, 9.5);

   TH1F * Mllbb = new TH1F("Mllbb","invariant mass of 2l2b after step 6",100,0,2000);
   Mllbb->SetXTitle("M_{2l2b} (GeV)");
   TH1F * Mllbb_up = new TH1F("Mllbb_up","invariant mass of 2l2b after step 6",100,0,2000);
   TH1F * Mllbb_dow = new TH1F("Mllbb_dow","invariant mass of 2l2b after step 6",100,0,2000);
   TH1F * Mllbb_jer_up = new TH1F("Mllbb_jer_up","invariant mass of 2l2b after step 6",100,0,2000);
   TH1F * Mllbb_jer_dow = new TH1F("Mllbb_jer_dow","invariant mass of 2l2b after step 6",100,0,2000);

   TH1F * Mllbb_pdf_up = new TH1F("Mllbb_pdf_up","invariant mass of 2l2b after step 6",100,0,2000);
   TH1F * Mllbb_pdf_dow = new TH1F("Mllbb_pdf_dow","invariant mass of 2l2b after step 6",100,0,2000);
   TH1F * Mllbb_pu_up = new TH1F("Mllbb_pu_up","invariant mass of 2l2b after step 6",100,0,2000);
   TH1F * Mllbb_pu_dow = new TH1F("Mllbb_pu_dow","invariant mass of 2l2b after step 6",100,0,2000);

   TH1F * Mllbb_mc_1 = new TH1F("Mllbb_mc_1","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllbb_mc_2 = new TH1F("Mllbb_mc_2","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllbb_mc_3 = new TH1F("Mllbb_mc_3","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllbb_mc_4 = new TH1F("Mllbb_mc_4","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllbb_mc_5 = new TH1F("Mllbb_mc_5","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllbb_mc_6 = new TH1F("Mllbb_mc_6","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllbb_mc_7 = new TH1F("Mllbb_mc_7","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllbb_mc_8 = new TH1F("Mllbb_mc_8","invariant mass of jet pair after step 6",100,0,2000);

   TH1F * Mllbb_hh = new TH1F("Mllbb_hh","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllbb_hl = new TH1F("Mllbb_hl","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllbb_ll = new TH1F("Mllbb_ll","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllbb_hh_up = new TH1F("Mllbb_hh_up","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllbb_hl_up = new TH1F("Mllbb_hl_up","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllbb_ll_up = new TH1F("Mllbb_ll_up","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllbb_hh_dow = new TH1F("Mllbb_hh_dow","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllbb_hl_dow = new TH1F("Mllbb_hl_dow","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllbb_ll_dow = new TH1F("Mllbb_ll_dow","invariant mass of bottom pair after step 6",100,0,2000);

   TH1F * Mllllbb = new TH1F("Mllllbb","invariant mass of 4l2b after step 6",100,0,2000);
   Mllllbb->SetXTitle("M_{4l2b} (GeV)");
   TH1F * Mllllbb_up = new TH1F("Mllllbb_up","invariant mass of 4l2b after step 6",100,0,2000);
   TH1F * Mllllbb_dow = new TH1F("Mllllbb_dow","invariant mass of 4l2b after step 6",100,0,2000);
   TH1F * Mllllbb_jer_up = new TH1F("Mllllbb_jer_up","invariant mass of 4l2b after step 6",100,0,2000);
   TH1F * Mllllbb_jer_dow = new TH1F("Mllllbb_jer_dow","invariant mass of 4l2b after step 6",100,0,2000);

   TH1F * Mllllbb_pdf_up = new TH1F("Mllllbb_pdf_up","invariant mass of 4l2b after step 6",100,0,2000);
   TH1F * Mllllbb_pdf_dow = new TH1F("Mllllbb_pdf_dow","invariant mass of 4l2b after step 6",100,0,2000);
   TH1F * Mllllbb_pu_up = new TH1F("Mllllbb_pu_up","invariant mass of 4l2b after step 6",100,0,2000);
   TH1F * Mllllbb_pu_dow = new TH1F("Mllllbb_pu_dow","invariant mass of 4l2b after step 6",100,0,2000);

   TH1F * Mllllbb_mc_1 = new TH1F("Mllllbb_mc_1","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllllbb_mc_2 = new TH1F("Mllllbb_mc_2","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllllbb_mc_3 = new TH1F("Mllllbb_mc_3","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllllbb_mc_4 = new TH1F("Mllllbb_mc_4","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllllbb_mc_5 = new TH1F("Mllllbb_mc_5","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllllbb_mc_6 = new TH1F("Mllllbb_mc_6","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllllbb_mc_7 = new TH1F("Mllllbb_mc_7","invariant mass of jet pair after step 6",100,0,2000);
   TH1F * Mllllbb_mc_8 = new TH1F("Mllllbb_mc_8","invariant mass of jet pair after step 6",100,0,2000);

   TH1F * Mllllbb_hh = new TH1F("Mllllbb_hh","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllllbb_hl = new TH1F("Mllllbb_hl","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllllbb_ll = new TH1F("Mllllbb_ll","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllllbb_hh_up = new TH1F("Mllllbb_hh_up","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllllbb_hl_up = new TH1F("Mllllbb_hl_up","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllllbb_ll_up = new TH1F("Mllllbb_ll_up","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllllbb_hh_dow = new TH1F("Mllllbb_hh_dow","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllllbb_hl_dow = new TH1F("Mllllbb_hl_dow","invariant mass of bottom pair after step 6",100,0,2000);
   TH1F * Mllllbb_ll_dow = new TH1F("Mllllbb_ll_dow","invariant mass of bottom pair after step 6",100,0,2000);

   TH2F * mass2D = new TH2F("mass2D","mass of mA and mH",100,0,1000,120,0,1200);
   mass2D->SetXTitle("M_{A} (GeV)");
   mass2D->SetYTitle("M_{H} (GeV)");
   TH2F * mass2D_up = new TH2F("mass2D_up","invariant mass of 4l2b after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_dow = new TH2F("mass2D_dow","invariant mass of 4l2b after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_jer_up = new TH2F("mass2D_jer_up","invariant mass of 4l2b after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_jer_dow = new TH2F("mass2D_jer_dow","invariant mass of 4l2b after step 6",100,0,1000,120,0,1200);

   TH2F * mass2D_pdf_up = new TH2F("mass2D_pdf_up","invariant mass of 4l2b after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_pdf_dow = new TH2F("mass2D_pdf_dow","invariant mass of 4l2b after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_pu_up = new TH2F("mass2D_pu_up","invariant mass of 4l2b after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_pu_dow = new TH2F("mass2D_pu_dow","invariant mass of 4l2b after step 6",100,0,1000,120,0,1200);

   TH2F * mass2D_mc_1 = new TH2F("mass2D_mc_1","invariant mass of jet pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_mc_2 = new TH2F("mass2D_mc_2","invariant mass of jet pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_mc_3 = new TH2F("mass2D_mc_3","invariant mass of jet pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_mc_4 = new TH2F("mass2D_mc_4","invariant mass of jet pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_mc_5 = new TH2F("mass2D_mc_5","invariant mass of jet pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_mc_6 = new TH2F("mass2D_mc_6","invariant mass of jet pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_mc_7 = new TH2F("mass2D_mc_7","invariant mass of jet pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_mc_8 = new TH2F("mass2D_mc_8","invariant mass of jet pair after step 6",100,0,1000,120,0,1200);

   TH2F * mass2D_hh = new TH2F("mass2D_hh","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_hl = new TH2F("mass2D_hl","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_ll = new TH2F("mass2D_ll","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_hh_up = new TH2F("mass2D_hh_up","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_hl_up = new TH2F("mass2D_hl_up","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_ll_up = new TH2F("mass2D_ll_up","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_hh_dow = new TH2F("mass2D_hh_dow","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_hl_dow = new TH2F("mass2D_hl_dow","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);
   TH2F * mass2D_ll_dow = new TH2F("mass2D_ll_dow","invariant mass of bottom pair after step 6",100,0,1000,120,0,1200);


   TTree *newtree  = new TTree("HZZ4LeptonsAnalysisReduced", "reduced ttree");
  
   // Add branches to output rootuple 
   Float_t f_weight, f_int_weight, f_pu_weight, f_eff_weight, f_lept1_pt, f_lept1_eta, f_lept1_phi, f_lept1_charge, f_lept1_pfx, f_lept1_sip, f_lept1_mvaid, f_lept2_pt, f_lept2_eta, f_lept2_phi, f_lept2_charge, f_lept2_pfx, f_lept2_sip, f_lept2_mvaid, f_lept3_pt, f_lept3_eta, f_lept3_phi, f_lept3_charge, f_lept3_pfx, f_lept3_sip, f_lept3_mvaid, f_lept4_pt, f_lept4_eta, f_lept4_phi, f_lept4_charge, f_lept4_pfx, f_lept4_sip, f_lept4_mvaid, f_iso_max, f_sip_max, f_Z1mass, f_Z2mass, f_angle_costhetastar, f_angle_costheta1, f_angle_costheta2, f_angle_phi, f_angle_phistar1, f_eta4l, f_pt4l, f_mass4l, f_mass4lErr, f_njets_pass, f_deltajj, f_massjj, f_D_jet, f_jet1_pt, f_jet1_eta, f_jet1_phi, f_jet1_e, f_jet2_pt, f_jet2_eta, f_jet2_phi, f_jet2_e;
   Float_t f_D_bkg_kin,f_D_bkg,f_D_gg,f_D_g4,f_Djet_VAJHU; 
   Float_t f_genmet, f_pfmet,f_mT,f_dphi,f_min_dphi_jet_met,f_max_dphi_jet_met,f_dphi_jet_met,f_mbb,f_m4l2b,f_m2l2b;
   Int_t f_lept1_pdgid,f_lept2_pdgid,f_lept3_pdgid,f_lept4_pdgid;
   Int_t f_category,f_Ngood,f_Nbjets,f_Njets,f_NVBFjets,f_NMatchbjets,f_NVHjets;
   Int_t f_run, f_lumi, f_event;
   
   TBranch *b_run= newtree->Branch("f_run", &f_run,"f_run/I");
   TBranch *b_lumi= newtree->Branch("f_lumi", &f_lumi,"f_lumi/I");    
   TBranch *b_event= newtree->Branch("f_event", &f_event,"f_event/I");    
   
   TBranch *b_weight= newtree->Branch("f_weight", &f_weight,"f_weight/F");
   TBranch *b_int_weight= newtree->Branch("f_int_weight", &f_int_weight,"f_int_weight/F");
   TBranch *b_pu_weight= newtree->Branch("f_pu_weight", &f_pu_weight,"f_pu_weight/F");
   TBranch *b_eff_weight= newtree->Branch("f_eff_weight", &f_eff_weight,"f_eff_weight/F");
   TBranch *b_lept1_pt= newtree->Branch("f_lept1_pt", &f_lept1_pt,"f_lept1_pt/F");
   TBranch *b_lept1_eta= newtree->Branch("f_lept1_eta", &f_lept1_eta,"f_lept1_eta/F");
   TBranch *b_lept1_phi= newtree->Branch("f_lept1_phi", &f_lept1_phi,"f_lept1_phi/F");
   TBranch *b_lept1_charge= newtree->Branch("f_lept1_charge", &f_lept1_charge,"f_lept1_charge/F");
   TBranch *b_lept1_pfx= newtree->Branch("f_lept1_pfx", &f_lept1_pfx,"f_lept1_pfx/F");
   TBranch *b_lept1_sip= newtree->Branch("f_lept1_sip", &f_lept1_sip,"f_lept1_sip/F");
   TBranch *b_lept1_pdgid= newtree->Branch("f_lept1_pdgid", &f_lept1_pdgid,"f_lept1_pdgid/I");
   TBranch *b_lept2_pt= newtree->Branch("f_lept2_pt", &f_lept2_pt,"f_lept2_pt/F");
   TBranch *b_lept2_eta= newtree->Branch("f_lept2_eta", &f_lept2_eta,"f_lept2_eta/F");
   TBranch *b_lept2_phi= newtree->Branch("f_lept2_phi", &f_lept2_phi,"f_lept2_phi/F");
   TBranch *b_lept2_charge= newtree->Branch("f_lept2_charge", &f_lept2_charge,"f_lept2_charge/F");
   TBranch *b_lept2_pfx= newtree->Branch("f_lept2_pfx", &f_lept2_pfx,"f_lept2_pfx/F");
   TBranch *b_lept2_sip= newtree->Branch("f_lept2_sip", &f_lept2_sip,"f_lept2_sip/F");
   TBranch *b_lept2_pdgid= newtree->Branch("f_lept2_pdgid", &f_lept2_pdgid,"f_lept2_pdgid/I");
   TBranch *b_lept3_pt= newtree->Branch("f_lept3_pt", &f_lept3_pt,"f_lept3_pt/F");
   TBranch *b_lept3_eta= newtree->Branch("f_lept3_eta", &f_lept3_eta,"f_lept3_eta/F");
   TBranch *b_lept3_phi= newtree->Branch("f_lept3_phi", &f_lept3_phi,"f_lept3_phi/F");
   TBranch *b_lept3_charge= newtree->Branch("f_lept3_charge", &f_lept3_charge,"f_lept3_charge/F");
   TBranch *b_lept3_pfx= newtree->Branch("f_lept3_pfx", &f_lept3_pfx,"f_lept3_pfx/F");
   TBranch *b_lept3_sip= newtree->Branch("f_lept3_sip", &f_lept3_sip,"f_lept3_sip/F");
   TBranch *b_lept3_pdgid= newtree->Branch("f_lept3_pdgid", &f_lept3_pdgid,"f_lept3_pdgid/I");
   TBranch *b_lept4_pt= newtree->Branch("f_lept4_pt", &f_lept4_pt,"f_lept4_pt/F");
   TBranch *b_lept4_eta= newtree->Branch("f_lept4_eta", &f_lept4_eta,"f_lept4_eta/F");
   TBranch *b_lept4_phi= newtree->Branch("f_lept4_phi", &f_lept4_phi,"f_lept4_phi/F");
   TBranch *b_lept4_charge= newtree->Branch("f_lept4_charge", &f_lept4_charge,"f_lept4_charge/F");
   TBranch *b_lept4_pfx= newtree->Branch("f_lept4_pfx", &f_lept4_pfx,"f_lept4_pfx/F");
   TBranch *b_lept4_sip= newtree->Branch("f_lept4_sip", &f_lept4_sip,"f_lept4_sip/F");
   TBranch *b_lept4_pdgid= newtree->Branch("f_lept4_pdgid", &f_lept4_pdgid,"f_lept4_pdgid/I");
   TBranch *b_Z1mass= newtree->Branch("f_Z1mass", &f_Z1mass,"f_Z1mass/F");
   TBranch *b_Z2mass= newtree->Branch("f_Z2mass", &f_Z2mass,"f_Z2mass/F");
   TBranch *b_angle_phi= newtree->Branch("f_angle_phi", &f_angle_phi,"f_angle_phi/F");
   TBranch *b_njets_pass= newtree->Branch("f_njets_pass", &f_njets_pass,"f_njets_pass/F");
   TBranch *b_deltajj= newtree->Branch("f_deltajj", &f_deltajj,"f_deltajj/F");
   TBranch *b_massjj= newtree->Branch("f_massjj", &f_massjj,"f_massjj/F");
   TBranch *b_massbb= newtree->Branch("f_mbb", &f_mbb,"f_mbb/F");
   TBranch *b_mass4l2b= newtree->Branch("f_m4l2b", &f_m4l2b,"f_m4l2b/F");
   TBranch *b_mass2l2b= newtree->Branch("f_m2l2b", &f_m2l2b,"f_m2l2b/F");
   TBranch *b_jet1_pt= newtree->Branch("f_jet1_pt", &f_jet1_pt,"f_jet1_pt/F");
   TBranch *b_jet1_eta= newtree->Branch("f_jet1_eta", &f_jet1_eta,"f_jet1_eta/F");
   TBranch *b_jet1_phi= newtree->Branch("f_jet1_phi", &f_jet1_phi,"f_jet1_phi/F");
   TBranch *b_jet1_e= newtree->Branch("f_jet1_e", &f_jet1_e,"f_jet1_e/F");
   TBranch *b_jet2_pt= newtree->Branch("f_jet2_pt", &f_jet2_pt,"f_jet2_pt/F");
   TBranch *b_jet2_eta= newtree->Branch("f_jet2_eta", &f_jet2_eta,"f_jet2_eta/F");
   TBranch *b_jet2_phi= newtree->Branch("f_jet2_phi", &f_jet2_phi,"f_jet2_phi/F");
   TBranch *b_jet2_e= newtree->Branch("f_jet2_e", &f_jet2_e,"f_jet2_e/F");
   TBranch *b_D_bkg_kin= newtree->Branch("f_D_bkg_kin", &f_D_bkg_kin,"f_D_bkg_kin/F");
   TBranch *b_D_bkg= newtree->Branch("f_D_bkg", &f_D_bkg,"f_D_bkg/F");
   TBranch *b_D_gg= newtree->Branch("f_D_gg", &f_D_gg,"f_D_gg/F");
   TBranch *b_D_g4= newtree->Branch("f_D_g4", &f_D_g4,"f_D_g4/F");
   TBranch *b_Djet_VAJHU= newtree->Branch("f_Djet_VAJHU", &f_Djet_VAJHU,"f_Djet_VAJHU/F");
   TBranch *b_pfmet= newtree->Branch("f_pfmet", &f_pfmet,"f_pfmet/F");
   TBranch *b_genmet= newtree->Branch("f_genmet", &f_genmet,"f_genmet/F");
   TBranch *b_f_mT= newtree->Branch("f_mT", &f_mT,"f_mT/F");
   TBranch *b_f_dphi= newtree->Branch("f_dphi", &f_dphi,"f_dphi/F");
   TBranch *b_f_category=newtree->Branch("f_category", &f_category, "f_category/I");
   TBranch *b_f_Ngood=newtree->Branch("f_Ngood", &f_Ngood, "f_Ngood/I");
   TBranch *b_f_Nbjets=newtree->Branch("f_Nbjets", &f_Nbjets, "f_Nbjets/I");
   TBranch *b_f_Njets=newtree->Branch("f_Njets", &f_Njets, "f_Njets/I");
   TBranch *b_f_NVBFjets=newtree->Branch("f_NVBFjets", &f_NVBFjets, "f_NVBFjets/I");
   TBranch *b_f_NMatchbjets=newtree->Branch("f_NMatchbjets", &f_NMatchbjets, "f_NMatchbjets/I");
   TBranch *b_f_dphi_jet_met= newtree->Branch("f_dphi_jet_met", &f_dphi_jet_met,"f_dphi_jet_met/F");
   TBranch *b_f_min_dphi_jet_met= newtree->Branch("f_min_dphi_jet_met", &f_min_dphi_jet_met,"f_min_dphi_jet_met/F");
   TBranch *b_f_max_dphi_jet_met= newtree->Branch("f_max_dphi_jet_met", &f_max_dphi_jet_met,"f_max_dphi_jet_met/F");

   float newweight=1.;
   
   // New tree with clone of events passing the final selection
   // Clone tree for final events
//   TTree *finaltree = fChain->CloneTree(0);

   // loop on entries
   
   Long64_t nentries = fChain->GetEntries();

   cout << "\n****************************"  <<endl;
   cout << "Analyzing " << nentries << " entries"  <<endl;     
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

 // Initialize reduced tree variables
      f_weight=1., f_int_weight=1., f_pu_weight=1., f_eff_weight=1., f_lept1_pt=-999., f_lept1_eta=-999., f_lept1_phi=-999., f_lept1_charge=-999., f_lept1_pfx=-999., f_lept1_sip=-999., f_lept1_mvaid=-999., f_lept2_pt=-999., f_lept2_eta=-999., f_lept2_phi=-999., f_lept2_charge=-999., f_lept2_pfx=-999., f_lept2_sip=-999., f_lept2_mvaid=-999., f_lept3_pt=-999., f_lept3_eta=-999., f_lept3_phi=-999., f_lept3_charge=-999., f_lept3_pfx=-999., f_lept3_sip=-999., f_lept3_mvaid=-999., f_lept4_pt=-999., f_lept4_eta=-999., f_lept4_phi=-999., f_lept4_charge=-999., f_lept4_pfx=-999., f_lept4_sip=-999., f_lept4_mvaid=-999., f_iso_max=-999., f_sip_max=-999., f_Z1mass=-999., f_Z2mass=-999., f_angle_costhetastar=-999., f_angle_costheta1=-999., f_angle_costheta2=-999., f_angle_phi=-999., f_angle_phistar1=-999., f_eta4l=-999., f_pt4l=-999., f_mass4l=-999., f_mass4lErr=-999., f_njets_pass=-999., f_deltajj=-999., f_massjj=-999., f_D_jet=-999., f_jet1_pt=-999., f_jet1_eta=-999., f_jet1_phi=-999., f_jet1_e=-999., f_jet2_pt=-999., f_jet2_eta=-999., f_jet2_phi=-999., f_jet2_e=-999.,f_D_bkg_kin=-999.,f_D_bkg=-999.,f_D_gg=-999.,f_D_g4=-999.,f_Djet_VAJHU=-999.,f_genmet=-999.,f_pfmet=-999.,f_mT=-999.,f_dphi=-999.,f_lept1_pdgid=-999,f_lept2_pdgid=-999,f_lept3_pdgid=-999,f_lept4_pdgid=-999,f_run=-999, f_lumi=-999, f_event=-999, f_category=-999, f_Ngood=-999, f_Nbjets=-999, f_NVBFjets=-999,f_NMatchbjets=-999,f_NVHjets=-999,f_dphi_jet_met=-999., f_min_dphi_jet_met=-999., f_max_dphi_jet_met=-999.; 
      f_Njets=-999;
//      f_u1=-999, f_u2=-999;
      f_m4l2b=-999;
      f_m2l2b=-999;
//      if (!(Run==274422 && LumiSection==484 && Event==843958845)) continue;

      if(jentry%5000 == 0) cout << "Analyzing entry: " << jentry << endl;
      
      double mc_weight_un[9];

      if( RECO_NMU > 100 ) RECO_NMU = 100;
      if( RECO_NELE > 100 ) RECO_NELE = 100;
      if( RECO_NPFPHOT > 20 ) RECO_NPFPHOT = 20;
      
      bool debug=false;  //debug flag  -- default false

      newweight=weight;
      if(jentry%5000 == 0) cout << "Starting weight= " << newweight << endl;

      // pileup reweighting 2012 and 2011
      if (DATA_type=="NO" && num_PU_vertices < 0) continue;                                                                                                                                              
      // pileup reweighting 2015
      hPUvertices->Fill(num_PU_vertices,weight);

      double pu_weight=1.;
      double pu_up=1;
      double pu_dow=1;
      if (MC_type == "Spring16"){
	Int_t binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
	if(debug) cout << " bin x= " << binx << " " << puweight->GetBinContent(binx) << endl;	
	pu_weight=double(puweight->GetBinContent(binx));
        pu_up=double(puweight_up->GetBinContent(binx))/pu_weight;
        pu_dow=double(puweight_dow->GetBinContent(binx))/pu_weight;	
      }      
       
      hPUvertices_ReWeighted->Fill(num_PU_vertices,weight*pu_weight);
      if(debug) cout << "Pileup interations and weight is= " << num_PU_vertices << " " << " and weight= " << pu_weight << endl;  
      
      //if (num_PU_vertices < 0) continue;

      // Changing the weight for pileup
      newweight=weight*pu_weight;
      if(jentry%5000 == 0) cout << "Starting weight + pileup = " << newweight << endl;
           
      

      // Weight for MCNLO samples                                                                                      
      if( datasetName.Contains("amcatnlo")) {
        if(jentry%5000 == 0) cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting << endl;
        newweight=weight*pu_weight*MC_weighting;
      }

      for(int l=0; l<9; l++){
         if(MC_weighting_un[0]!=0) mc_weight_un[l]=MC_weighting_un[l];
         else mc_weight_un[l]=1;
      }     
 
      float pFill[11];for(int pf=0;pf<11;pf++)pFill[11]=-999.;

      // ** Step 0:
      // simply number of entries...
      if( debug ) cout << "\n** Step 0: \nAnalyzing entry: " << jentry << " Run: " << Run << " Event: " << Event << " LumiSection: " << LumiSection << endl ;
      ++N_0 ;  // fill counter
      N_0_w=N_0_w+newweight;
      
      // ** Step 0.1:
      //trigger requirements (qier)
     if(!de_trig) continue; 
     
      bool is2e2mu=false;

/*      
//teleport 1      
      //////////// (VBF) Line Shape Correction///////////
      if ( (isggH || isvbf) && mHgen >=400. && (useLineShape==true || useVBFLineShape==true) ){
	cout << "Applying lineshape correction for masses= " << mHgen << " " << MC_MASS[0] << endl;
	//cout << bincenters_.front() << " " <<  bincenters_.back() << endl;
	if( MC_MASS[0] < bincenters_.front() || MC_MASS[0] >  bincenters_.back() ){ // set weights to 0 if out of range
	  cout << " out of range " << endl;
	  LineShapeWeight = 0.;
	  LineShapeWeightP = 0.;
	  LineShapeWeightM = 0.;
	}
	
	std::vector<double>::iterator low;
	low=lower_bound( bincenters_.begin(), bincenters_.end(),MC_MASS[0] ); 
	int lowindex=(low-  bincenters_.begin());
	//cout << "low " << bincenters_.begin() << endl;
	if(MC_MASS[0] == *low ){//exact match
	  cout << "Exact match " << endl;
	  if (useLineShape==true) (LineShapeWeight) = weightLineShape_[lowindex];
	  if (useVBFLineShape==true) (LineShapeWeight) = weightVBFLineShape_[lowindex]; 
	}
	else{//linear interpolation
	  lowindex--; // lower_bound finds the first element not smaller than X
	  if (useLineShape==true){ 
	    cout << "use LineShape" << endl;
	    LineShapeWeight = weightLineShape_[lowindex] +(MC_MASS[0]  - bincenters_[lowindex] )*(weightLineShape_[lowindex+1]-weightLineShape_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex]);
	    cout << "LineShape weight= " << LineShapeWeight << endl;
	    LineShapeWeightP = weightLineShapeP_[lowindex] +(MC_MASS[0]  - bincenters_[lowindex] )*(weightLineShapeP_[lowindex+1]-weightLineShape_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex]);
	    LineShapeWeightM = weightLineShapeM_[lowindex] +(MC_MASS[0]  - bincenters_[lowindex] )*(weightLineShapeM_[lowindex+1]-weightLineShape_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex]);
	  }
	  if (useVBFLineShape==true){
	    cout << "use VBF LineShape" << endl;
	    LineShapeWeight = weightVBFLineShape_[lowindex] +(MC_MASS[0]  - bincenters_[lowindex] )*(weightVBFLineShape_[lowindex+1]-weightVBFLineShape_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex]);
	    cout << "VBF LineShape weight= " << LineShapeWeight << endl;
	    LineShapeWeightP = weightVBFLineShapeP_[lowindex] +(MC_MASS[0]  - bincenters_[lowindex] )*(weightVBFLineShapeP_[lowindex+1]-weightVBFLineShape_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex]);
	    LineShapeWeightM = weightVBFLineShapeM_[lowindex] +(MC_MASS[0]  - bincenters_[lowindex] )*(weightVBFLineShapeM_[lowindex+1]-weightVBFLineShape_[lowindex])/(bincenters_[lowindex+1]-bincenters_[lowindex]);
	  }
	}

	if (useLineShape==true){ 
	  cout<<"Higgs Mass==== "<< MC_MASS[0]  <<" LineShapeWeight==== "<<LineShapeWeight<<endl;       
	  Gen_H_MASS->Fill(MC_MASS[0]);
	  if(MC_MASS[0]>10.)Gen_H_MASS_ReWeighted->Fill(MC_MASS[0],LineShapeWeight);
	  if(MC_MASS[0]>10.)Gen_H_MASS_ReWeightedP->Fill(MC_MASS[0],LineShapeWeightP);
	  if(MC_MASS[0]>10.)Gen_H_MASS_ReWeightedM->Fill(MC_MASS[0],LineShapeWeightM);
	}
	if (useVBFLineShape==true){ 
	  cout<<"Higgs Mass==== "<< MC_MASS[0]  <<" VBFLineShapeWeight==== "<<LineShapeWeight<<endl;       
	  Gen_H_MASS->Fill(MC_MASS[0]);
	  if(MC_MASS[0]>10.)Gen_H_MASS_ReWeighted->Fill(MC_MASS[0],LineShapeWeight);
	  if(MC_MASS[0]>10.)Gen_H_MASS_ReWeightedP->Fill(MC_MASS[0],LineShapeWeightP);
	  if(MC_MASS[0]>10.)Gen_H_MASS_ReWeightedM->Fill(MC_MASS[0],LineShapeWeightM);
	}	

	// Changing the weight for pileup and LineShape
	newweight=weight*pu_weight*LineShapeWeight;
	cout << "Starting weight + pileup + LineShape= " << newweight << endl;
      }
      //teleport1
*/
/*
      if (isSignal){
        if( ( MC_PDGID[3] == 13 && MC_PDGID[4] == -13  &&  MC_PDGID[5] == 11 && MC_PDGID[6] == -11 )
            || ( MC_PDGID[3] == 11 && MC_PDGID[4] == -11  &&  MC_PDGID[5] == 13 && MC_PDGID[6] == -13 )){
          is2e2mu=true;
          ++N_01 ;  // fill counter                                                                                                                                             
          N_01_w=N_01_w+newweight;
        }
        else continue;
	
	
	// ** Step 0.2:
	// events in the ACCEPTANCE:
	// variables contributing: eta & pt
	// NEW only eta cut
	bool isInAcceptance = 1;
	
	for( int i = 3; i  < 7; ++i){
	  
	  if( debug ) cout << "\n ** Step 0.2 (acceptance): "
			   << "\n fabs(MC_ETA[i])" << fabs(MC_ETA[i])	
			   << "\n MC_PT[i]" << MC_PT[i]	
			   << endl ;
	  
          if( fabs(MC_PDGID[i]) == 13 && ( fabs(MC_ETA[i]) >= 2.4) ){
            isInAcceptance = 0;
            break;        
          }
          if( fabs(MC_PDGID[i]) == 11 && ( fabs(MC_ETA[i]) >= 2.5) ){
	    
	    isInAcceptance = 0;
	    break;        
	  }
	  
	}// end for leptons
	
	if( is2e2mu && isInAcceptance )       {
	  ++N_02 ;  // fill counter
	  N_02_w=N_02_w+newweight;
	}
	else continue;	
      }
    */ 
      ++N_1 ;  // fill counter
      N_1_w=N_1_w+newweight;
  
      bool tag_2011=false;
      if (DATA_type=="2010" || DATA_type=="2011" || MC_type=="Fall11"){
        tag_2011=true;
      }
    
      // Loose lepton identification
      
      int N_loose_mu = 0;
      int N_loose_e = 0;
      double max_Iso_loose_mu = -1 ;
      double max_Sip_loose_mu = -1 ;
      double max_Ip_loose_mu = -1 ;
      double max_Iso_loose_e = -1 ;
      double max_Sip_loose_e = -1 ;
      double max_Ip_loose_e = -1 ;
      
      int* arraysize_mu = new int[1];
      arraysize_mu[0] = RECO_NMU;
      int iL_loose_mu[arraysize_mu[0]];
      delete [] arraysize_mu;

      for( int i = 0; i < RECO_NMU; ++i ){
	iL_loose_mu[i]=-999.;
      }
 
      for( int i = 0; i < RECO_NMU; ++i ){

        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOMU_isGlobalMu[i] " << int(RECOMU_isGlobalMu[i])
			 << "\nRECOMU_isTrackerMu[i] " << int(RECOMU_isTrackerMu[i])
			 << "\nRECOMU_PT[i] " << RECOMU_PT[i]
			 << "\nfabs(RECOMU_ETA[i]) " << fabs(RECOMU_ETA[i])
			 << "\nfabs( RECOMU_mubestrkDxy[i] ) " << fabs( RECOMU_mubesttrkDxy[i] )
			 << "\nfabs( RECOMU_mubesttrkDz[i] ) " << fabs( RECOMU_mubesttrkDz[i] )
			 << endl ;
       	
 	if(
              ( RECOMU_isGlobalMu[i] || RECOMU_isTrackerMu[i] ) && RECOMU_isPFMu[i]
              && RECOMU_PT[i] > 5.
              && fabs(RECOMU_ETA[i]) < 2.4
              && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1. //loose muon (qier) 
	    ){ 
	  iL_loose_mu[N_loose_mu]=i;
	  ++N_loose_mu ;
	  if( RECOMU_PFX_dB[i] > max_Iso_loose_mu ) max_Iso_loose_mu = RECOMU_PFX_dB[i] ;
	  if( fabs( RECOMU_SIP[i] ) > max_Sip_loose_mu ) max_Sip_loose_mu = fabs( RECOMU_SIP[i] ) ;
	  if( fabs( RECOMU_IP[i] ) > max_Ip_loose_mu ) max_Ip_loose_mu = fabs( RECOMU_IP[i] ) ;
	}
	
      } // end loop on muons

      
      int* arraysize_e = new int[1];
      arraysize_e[0] = RECO_NELE;
      int iL_loose_e[arraysize_e[0]];
      delete [] arraysize_e;

      for( int i = 0; i < RECO_NELE; ++i ){
	iL_loose_e[i]=-999.;
      }

      for( int i = 0; i < RECO_NELE; ++i ){
	
        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOELE_PT[i] " << RECOELE_PT[i]
			 << "\nfabs(RECOELE_ETA[i]) " << fabs(RECOELE_ETA[i])
			 << "\nfabs( RECOELE_gsftrack_dxy[i] ) " << fabs( RECOELE_gsftrack_dxy[i] )
			 << "\nfabs( RECOELE_gsftrack_dz[i] ) " << fabs( RECOELE_gsftrack_dz[i] )
			 << endl ;
       	
 	if( RECOELE_PT[i] > 5. 
	    && fabs(RECOELE_ETA[i]) < 2.5 
	    // && RECOELE_gsftrack_expected_inner_hits[i]<=1  not used anymore
	    && fabs(RECOELE_gsftrack_dxy[i]) < .5 
	    && fabs(RECOELE_gsftrack_dz[i]) < 1. 
	    ) {	  
	  iL_loose_e[N_loose_e]=i;
	  ++N_loose_e ;
	  if( RECOELE_PFX_rho[i] > max_Iso_loose_e ) max_Iso_loose_e = RECOELE_PFX_rho[i] ;
	  if( fabs( RECOELE_SIP[i] ) > max_Sip_loose_e ) max_Sip_loose_e = fabs( RECOELE_SIP[i] ) ;
	  if( fabs( RECOELE_IP[i] ) > max_Ip_loose_e ) max_Ip_loose_e = fabs( RECOELE_IP[i] ) ;
	}
	
      }// end loop on electrons
      
      hN_loose_mu->Fill( N_loose_mu,newweight );
      hN_loose_e->Fill( N_loose_e,newweight );
      hIso_loose_mu->Fill( max_Iso_loose_mu,newweight);
      hSip_loose_mu->Fill( max_Sip_loose_mu,newweight );
      hIp_loose_mu->Fill( max_Ip_loose_mu,newweight );
      hIso_loose_e->Fill( max_Iso_loose_e,newweight );
      hSip_loose_e->Fill( max_Sip_loose_e,newweight );
      hIp_loose_e->Fill( max_Ip_loose_e,newweight );
      
      // Electron Cross Cleaning  -- eles separated from muons (deltaR > 0.05)
      
      for(int e = 0; e < RECO_NELE; ++e)
      	for(int mu = 0; mu < RECO_NMU; ++mu){
	  
	  if( 
              ( RECOMU_isGlobalMu[mu] || RECOMU_isTrackerMu[mu] ) && RECOMU_isPFMu[mu]
	      && RECOMU_PT[mu] > 5. 
	      && fabs(RECOMU_ETA[mu]) < 2.4 
	      && fabs(RECOMU_mubesttrkDxy[mu]) < .5 && fabs(RECOMU_mubesttrkDz[mu]) < 1. 
	      );
	  else continue;
	  
	  double deltaR = sqrt( pow( DELTAPHI( RECOMU_PHI[mu] , RECOELE_PHI[e] ),2) + pow(RECOMU_ETA[mu] - RECOELE_ETA[e],2) );
	  
	  if( deltaR <= 0.05 ){
	    
	    if( debug )cout << "Electrom not passing the cross cleaning" << endl;
	    
	    RECOELE_PT[e]  = -0.01;
	    RECOELE_ETA[e] = -99.;
	    RECOELE_PHI[e] = -99.;
	    RECOELE_SIP[e] = -99.;
	  }
	}
      
                  
      // Lepton identification -- no iso
      
      int iL[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1};
      
      int N_good = 0 ;
      
      for( int i = 0; i < RECO_NMU; ++i ){
	
        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOMU_isPFMu[i] " << int(RECOMU_isPFMu[i])
			 << "\nRECOMU_isGlobalMu[i] " << int(RECOMU_isGlobalMu[i])
			 << "\nRECOMU_isTrackerMu[i] " << int(RECOMU_isTrackerMu[i])
			 << "\nRECOMU_PT[i] " << RECOMU_PT[i]
			 << "\nfabs(RECOMU_ETA[i]) " << fabs(RECOMU_ETA[i])
			 << "\nRECOMU_PFX_dB[i] " << RECOMU_PFX_dB[i]
                         << "\nRECOMU_isMedium[i] " << RECOMU_isMedium[i]
			 << "\nfabs( RECOMU_SIP[i] ) " << fabs( RECOMU_SIP[i] )
			 << "\nfabs( RECOMU_mubesttrkDxy[i] ) " << fabs( RECOMU_mubesttrkDxy[i] )
			 << "\nfabs( RECOMU_mubesttrkDz[i] ) " << fabs( RECOMU_mubesttrkDz[i] )
			 << endl ;
	
       	// loose muons
 	if( 
              ( RECOMU_isGlobalMu[i] || RECOMU_isTrackerMu[i] ) && RECOMU_isPFMu[i]
              && RECOMU_PT[i] > 5.
              && fabs(RECOMU_ETA[i]) < 2.4
              && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1. //loose muon (qier)
	    ){

          Double_t Pt = RECOMU_PT[i];
          Double_t Eta = RECOMU_ETA[i];
          Int_t Q = int(RECOMU_CHARGE[i]);
          Double_t Phi = RECOMU_PHI[i];
          Double_t nl = RECOMU_mutrktrackerLayersWithMeasurement[i];
          if(debug) cout << "Q=" << Q << " Pt=" << Pt << " Eta=" << Eta << " Phi=" << Phi << " nl=" << nl << endl;
          if( MC_type == "Spring16" && DATA_type == "NO"){
          double u1 = gRandom->Rndm();
          double u2 = gRandom->Rndm();
          double mcSF = rc.kScaleAndSmearMC(Q,Pt,Eta,Phi,nl,u1,u2);
          if(debug) cout << "calibration SF = " << mcSF << endl;
          RECOMU_PT[i]=RECOMU_PT[i]*mcSF;
          }
        if( MC_type == "NO" && DATA_type == "2016"){
          double dataSF = rc.kScaleDT(Q,Pt,Eta,Phi);
          RECOMU_PT[i]=RECOMU_PT[i]*dataSF;
          if(debug) cout << "calibration SF = " << dataSF << endl;
        }
	  
	  iL[ N_good ] = i ;
	  ++N_good ;	  
	}
      } // end loop on muons
      
      hN_good_mu->Fill( N_good,newweight );
      
      if( debug ) cout << "\nLeptons' indeces: "
		       << "\niL[0]: " << iL[0]
		       << "\niL[1]: " << iL[1]
		       << "\niL[2]: " << iL[2]
		       << "\niL[3]: " << iL[3]
		       << "\niL[4]: " << iL[4]
		       << "\niL[5]: " << iL[5]
		       << "\niL[6]: " << iL[6]
		       << "\niL[7]: " << iL[7]
		       << "\nNumber of good muons: " << N_good
		       << endl ;
      

    
      
      /// *** FSR
      // Photon identification & cleaning
      // ele identification is also needed
      
      //FSR photon identifications, will be used with MELA later
      int FSR_Z1_photid=-1;
      int FSR_Z2_photid=-1;
      int FSR_Z1_lepid=-1;
      int FSR_Z2_lepid=-1;

      //electrons:
      
      int Ne_good = 0 ;
      int iLe[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1}; //electrons
      int Ne_loose_4=0;

      for( int i = 0; i < RECO_NELE; ++i ){

        if( debug ) cout << "\n Electron i="<< i <<" properties: "
      		  << "\nRECOELE_PT[i] " << RECOELE_PT[i]
      		  << "\nfabs(RECOELE_ETA[i]) " << fabs(RECOELE_ETA[i])
		  << "\nfabs(RECOELE_scl_Eta[i]) " << fabs(RECOELE_scl_Eta[i])
      		  << "\nRECOELE_PFX_rho[i] " << RECOELE_PFX_rho[i]
      		  << "\nfabs( RECOELE_SIP[i] ) " << fabs( RECOELE_SIP[i] )
      		  << "\nRECOELE_mvaNonTrigV0[i] " << RECOELE_mvaNonTrigV0[i]
		  << "\nfabs( RECOELE_gsftrack_dxy[i] ) " << fabs( RECOELE_gsftrack_dxy[i] )
      		  << "\nfabs( RECOELE_gsftrack_dz[i] ) " << fabs( RECOELE_gsftrack_dz[i] )
		  << endl ;
       	
 	if( RECOELE_PT[i] > 5. && fabs(RECOELE_ETA[i]) < 2.5 );
	  // && RECOELE_gsftrack_expected_inner_hits[i]<=1 ) /* ok */ ;
	else continue ;
	
        Ne_loose_4++;

	bool BDT_ok = 0; // Spring16 with CMSSW_8_0_x
/*	if( RECOELE_PT[i] > 7. &&  RECOELE_PT[i] <= 10. ){
		if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.211 ) BDT_ok = 1 ;
		if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) < 1.479 )
						 && RECOELE_mvaNonTrigV0[i] > -0.396 ) BDT_ok = 1 ;
		if( fabs(RECOELE_scl_Eta[i]) >= 1.479 && RECOELE_mvaNonTrigV0[i] > -0.215 ) BDT_ok = 1 ;
	}
	else { 
		if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.870 ) BDT_ok = 1 ;
		if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
						 && RECOELE_mvaNonTrigV0[i] > -0.838 ) BDT_ok = 1 ;
		if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > -0.763 ) BDT_ok = 1 ;
	}	*/
        //qier chenge to wp90
        if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaTrigV0[i] > 0.837 ) BDT_ok = 1 ;
        if( fabs(RECOELE_scl_Eta[i]) >= .8 &&fabs(RECOELE_scl_Eta[i])<1.479 && RECOELE_mvaTrigV0[i] > 0.715 ) BDT_ok = 1 ;
        if(fabs(RECOELE_scl_Eta[i])>=1.479 && RECOELE_mvaTrigV0[i] > 0.357) BDT_ok=1;
	if( !BDT_ok ) continue ;
	
	if( fabs(RECOELE_gsftrack_dxy[i]) < .5 
	 && fabs(RECOELE_gsftrack_dz[i])  < 1. ) /* ok */ ;
	else continue ; 
		
	iLe[ Ne_good ] = i ;
	++Ne_good ;

      }// end loop on electrons

      hN_good_ele->Fill( Ne_good,newweight );
      hN_good_lep->Fill( N_good + Ne_good,newweight );


      if( debug ) cout << "\n Electrons' indeces: "
 		  << "\niLe[0]: " << iLe[0]
  		  << "\niLe[1]: " << iLe[1]
 		  << "\niLe[2]: " << iLe[2]
		  << "\niLe[3]: " << iLe[3]
 		  << "\niLe[4]: " << iLe[4]
  		  << "\niLe[5]: " << iLe[5]
 		  << "\niLe[6]: " << iLe[6]
		  << "\niLe[7]: " << iLe[7]
		  << "\nNumber of good eles: " << Ne_good
		  << endl ;
            
      // Define a new isolation array to allocate the contribution of photons
      // float RECOMU_PFX_dB_new[100],RECOELE_PFX_rho_new[100];
      for (int i=0;i<100;i++){
	RECOMU_PFX_dB_new[i]=RECOMU_PFX_dB[i];
	RECOELE_PFX_rho_new[i]=RECOELE_PFX_rho[i];	
      }
      //


      // photon definition & cleaning:
      int iLp[30];
      	for( int i = 0 ; i < 30 ; ++i )iLp[i] = -1;
      
      int Nphotons = 0;

      for( int i = 0; i < RECO_NPFPHOT; ++i ){
	
        if( debug ) cout << "\n Photon i="<< i <<" properties: "
			 << "\n RECOPFPHOT_PT[i] " << RECOPFPHOT_PT[i]
			 << "\n fabs(RECOPFPHOT_ETA[i]) " << fabs(RECOPFPHOT_ETA[i])
			 << "\n RECOPFPHOT_PHI[i] " << RECOPFPHOT_PHI[i]
			 << "\n RECOPFPHOT_PFX_rho[i] " << RECOPFPHOT_PFX_rho[i]
			 << endl ;
	
	if ( RECOPFPHOT_PT[i] > 2. && fabs(RECOPFPHOT_ETA[i]) < 2.4 && RECOPFPHOT_PFX_rho[i]<1.8) {
	  
	  bool is_clean = 1;
	  
	  // cleaning
	  for(int e = 0; e < N_loose_e; ++e){
//  	    if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;  // loose ID + SIP cut	    
	    double deltaPhi = DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ) ;
	    double deltaEta = fabs( RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]] );
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	    if(debug) cout << "debug: " << RECOELE_PT[iL_loose_e[e]] << " " << deltaPhi << " " << deltaEta << " " << deltaR << endl;
	    if( ( fabs(deltaPhi) < 2 && fabs(deltaEta) < 0.05 ) || deltaR <= 0.15 ){		  
	      if( debug )cout << "Photon not passing the electron cleaning" << endl;	
	      is_clean = 0;	  
	      
	    }
	  } // end loop on eles		             	
	  
	  if( !is_clean ) continue ;
	  
	  
	  iLp[ Nphotons ] = i ;
	  ++Nphotons ;
	  
	}
      }// end loop on photons
      
      hN_good_phot->Fill( Nphotons,newweight );
      
      if( debug ) cout << "Photons' indeces: "
		       << "\niLp[0]: " << iLp[0]
		       << "\niLp[1]: " << iLp[1]
		       << "\niLp[2]: " << iLp[2]
		       << "\niLp[3]: " << iLp[3]
		       << "\niLp[4]: " << iLp[4]
		       << "\niLp[5]: " << iLp[5]
		       << "\niLp[6]: " << iLp[6]
		       << "\niLp[7]: " << iLp[7]
		       << "\nNumber of good photons: " << Nphotons
		       << endl ;
      
      
      // assign to each photon the closest lepton
      int iLp_l[30];
      for( int i = 0 ; i < 30 ; ++i )iLp_l[i] = -1;
      int iLp_tagEM[30];
      for( int i = 0 ; i < 30 ; ++i )iLp_tagEM[i] = -1;  // tag  0: mu  1: ele
      
      float RECOPFPHOT_DR[30];
      for( int i = 0 ; i < 30 ; ++i ) RECOPFPHOT_DR[i] = -999; 

      for( int i = 0; i < Nphotons; ++i ){
	
	double min_deltaR = 1000;
	int  l_min_deltaR = -1;
	int  tag_min_deltaR = -1;   // 0: mu  1: ele
	
	for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
//	  if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;  //loose ID + SIP cut
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	  if(debug) cout << "DeltaR= " << deltaR << " " <<  deltaR/pow(RECOPFPHOT_PT[iLp[i]],2) << endl;
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    if(debug) cout << "Possible candidate of photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " associated to a muon with pT= " << RECOMU_PT[iL_loose_mu[l]]<< endl;
	    min_deltaR = deltaR;
	    l_min_deltaR = l;
	    tag_min_deltaR = 0;
	  }
	  
	}//end loop on muons  
	
	for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
//	  if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue;  //loose ID + SIP cut
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[l]],2) );
	  if(debug) cout << "DeltaR= " << deltaR << " " <<  deltaR/pow(RECOPFPHOT_PT[iLp[i]],2) << endl;
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    if(debug) cout << "Possible candidate of photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " associated to an electron with pT= " << RECOELE_PT[iL_loose_e[l]]<< endl;
	    min_deltaR = deltaR;
	    l_min_deltaR = l;
	    tag_min_deltaR = 1;
	  }
	  
	}//end loop on electrons  

	
	if( min_deltaR < 0.5 ){
	  if (tag_min_deltaR==0) iLp_l[ i ] = iL_loose_mu[l_min_deltaR];
	  if (tag_min_deltaR==1) iLp_l[ i ] = iL_loose_e[l_min_deltaR];
	  iLp_tagEM[ i ] = tag_min_deltaR;
	  RECOPFPHOT_DR[iLp[i]]=min_deltaR; 	
	}

      }
      
      if( debug ) cout << "Indeces of loose leptons associated to photons: "
		       << "\niLp_l[0]: " << iLp_l[0]
		       << "\niLp_l[1]: " << iLp_l[1]
		       << "\niLp_l[2]: " << iLp_l[2]
		       << "\niLp_l[3]: " << iLp_l[3]
		       << "\niLp_l[4]: " << iLp_l[4]
		       << "\niLp_l[5]: " << iLp_l[5]
		       << "\niLp_l[6]: " << iLp_l[6]
		       << "\niLp_l[7]: " << iLp_l[7]
		       << endl ;
      
      if( debug ) cout << "Tag of leptons associated to photons: (0: mu , 1:ele)"
 		  << "\niLp_tagEM[0]: " << iLp_tagEM[0]
  		  << "\niLp_tagEM[1]: " << iLp_tagEM[1]
 		  << "\niLp_tagEM[2]: " << iLp_tagEM[2]
		  << "\niLp_tagEM[3]: " << iLp_tagEM[3]
 		  << "\niLp_tagEM[4]: " << iLp_tagEM[4]
  		  << "\niLp_tagEM[5]: " << iLp_tagEM[5]
 		  << "\niLp_tagEM[6]: " << iLp_tagEM[6]
		  << "\niLp_tagEM[7]: " << iLp_tagEM[7]
		  << endl ;


      // Multiple photons associated to the same lepton: the lowest-ΔR(γ,l)/ETγ2 has to be selected.
      double min_deltaR_ET2=1000;
      int p_min_deltaR_ET2=-1;

      for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
//	if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue; //loose ID + SIP cut
	min_deltaR_ET2=1000;
	p_min_deltaR_ET2=-1;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	    if(debug) cout <<  "index muon" << iL_loose_mu[l] << endl;
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	    double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	    if (deltaR_ET2<min_deltaR_ET2) {
	      min_deltaR_ET2=deltaR_ET2;
	      RECOPFPHOT_DR[iLp[p]]=deltaR;
	      p_min_deltaR_ET2=p;
	    }
	  }
	}
	
	if (p_min_deltaR_ET2!=-1){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	      if (p!=p_min_deltaR_ET2){
		iLp_l[ p ] = -1;
		iLp_tagEM[ p ] = -1;
	      }
	    }
	  }
	}
	
      }

   
      
      //
      min_deltaR_ET2=1000;
      p_min_deltaR_ET2=-1;
      
      for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
//	if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue; //loose ID + SIP cut
	min_deltaR_ET2=1000;
	p_min_deltaR_ET2=-1;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	    if(debug) cout <<  "index electron" << iL_loose_e[l] << endl;
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOELE_ETA[iL_loose_e[l]],2));
	    double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	    if(debug) cout << " deltaR_ET2= " << deltaR_ET2 <<endl;
	    if (deltaR_ET2<min_deltaR_ET2){
	      min_deltaR_ET2=deltaR_ET2;
	      RECOPFPHOT_DR[iLp[p]]=deltaR;
	      p_min_deltaR_ET2=p;
	      if(debug) cout << " p_min_deltaR_ET2= " << p_min_deltaR_ET2 <<endl;
	    }
	  }	  
	}
	
	if (p_min_deltaR_ET2!=-1){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	      if (p!=p_min_deltaR_ET2){
		iLp_l[ p ] = -1;
		iLp_tagEM[ p ] = -1;
	      }
	    }
	  }	  
	}
	
      }	
     
      
      if( debug ) cout << "Indeces of loose leptons associated to the photon with lowest DeltaR/ET2: "
		       << "\niLp_l[0]: " << iLp_l[0]
		       << "\niLp_l[1]: " << iLp_l[1]
		       << "\niLp_l[2]: " << iLp_l[2]
		       << "\niLp_l[3]: " << iLp_l[3]
		       << "\niLp_l[4]: " << iLp_l[4]
		       << "\niLp_l[5]: " << iLp_l[5]
		       << "\niLp_l[6]: " << iLp_l[6]
		       << "\niLp_l[7]: " << iLp_l[7]
		       << endl ;
      
      if( debug ) cout << "Tag of leptons associated to the photon with lowest DetaR/ET2: (0: mu , 1:ele)"
		       << "\niLp_tagEM[0]: " << iLp_tagEM[0]
		       << "\niLp_tagEM[1]: " << iLp_tagEM[1]
		       << "\niLp_tagEM[2]: " << iLp_tagEM[2]
		       << "\niLp_tagEM[3]: " << iLp_tagEM[3]
		       << "\niLp_tagEM[4]: " << iLp_tagEM[4]
		       << "\niLp_tagEM[5]: " << iLp_tagEM[5]
		       << "\niLp_tagEM[6]: " << iLp_tagEM[6]
		       << "\niLp_tagEM[7]: " << iLp_tagEM[7]
		       << endl ;

      
      for(int i=0.;i<Nphotons;i++) {
	if (iLp_l[i]!=-1 && iLp_tagEM[i]==0&&debug) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[i]] << endl;
	if (iLp_l[i]!=-1 && iLp_tagEM[i]==1&&debug) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to an electron with pT= " << RECOELE_PT[iLp_l[i]] << endl;
      };

       // Exclude that photon from the isolation cone all leptons in the event passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons
      
      if(debug) cout << "Rho for electron pileup isolation correction is= " << RHO_ele << endl;
      double EffectiveArea=-9999.;

	    
      for(int i=0.;i<Nphotons;i++) {
	if (iLp_l[i]==-1) continue;
	
	for(int e = 0; e < N_loose_e; ++e){
	  //if(!( iLp_l[i] == iL_loose_e[e] && iLp_tagEM[i] == 1 ) ) continue;
	  if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;
	  //double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[e]],2) );
	  if(debug) cout << "deltaR for photon subtraction= " << deltaR << endl;
	  if( deltaR<=0.3 && (RECOELE_scl_Eta[iL_loose_e[e]]< 1.479 || deltaR>0.08) ){ // 0.3 in 76x              
	    if( debug )cout << "Subtracting the photon isolation from the electron isolation value " << endl;
	    
	    EffectiveArea=EAele(iL_loose_e[e],tag_2011);
	    RECOELE_PFX_rho_new[iL_loose_e[e]]=
              (RECOELE_PFchHad[iL_loose_e[e]]+
               max(0.,RECOELE_PFneuHad[iL_loose_e[e]]+
                   (RECOELE_PFphoton[iL_loose_e[e]]-RECOPFPHOT_PT[iLp[i]] )-
                   max(RHO_ele,0.0)*(EffectiveArea)))/RECOELE_PT[iL_loose_e[e]];	    
	  }
	} // end loop on ele
	
	for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
	  //if(!( iLp_l[i] == iL_loose_mu[l] && iLp_tagEM[i] == 0 ) ) continue;
          if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;
          double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	  if(debug) cout << "deltaR for photon subtraction= " << deltaR << endl;
	  if( deltaR<=0.3 && deltaR>0.01){ // 0.3 is the isolation cone for muons in 76x
	    if( debug )cout << "Subtracting the photon isolation from the muon isolation value " << endl;
	    
	    RECOMU_PFX_dB_new[iL_loose_mu[l]]=
              (RECOMU_PFchHad[iL_loose_mu[l]]+
               max(0.,RECOMU_PFneuHad[iL_loose_mu[l]]+
                   (RECOMU_PFphoton[iL_loose_mu[l]]-RECOPFPHOT_PT[iLp[i]] )-
                   0.5*RECOMU_PFPUchAllPart[iL_loose_mu[l]]))/RECOMU_PT[iL_loose_mu[l]];
	    
	  }
	} // end loop on mu
	
	
	
      }	
      
      
      // *** end FSR
      
      
      // **** Step 3:
      struct candidateZ {
	float massvalue;
	int ilept1;
	float pt1;
	float isol1;
	bool ilept1_FSR;
	float eta1;
	float phi1;
	int charge1;
	int charge2;
	int ilept2;
	float pt2;
	float isol2;
	bool ilept2_FSR;
	float eta2;
	float phi2;
	float pxZ;
	float pyZ;
	float pzZ;
	float EZ;
	bool withFSR;
	float ptFSR;
	int tag;
      };

      vector<candidateZ> Zcandvector;
      Zcandvector.clear();
      vector<candidateZ> Zcandisolvector;
      Zcandisolvector.clear();

      // a) pair #1: mass closest to Z1
      // b) mLL in ] 40,120 [
      if( debug ) cout  << "\nStep 3: Number of good leptons: " << N_good+Ne_good << endl;

      if(Ne_loose_4+N_loose_mu < 4) continue;
      ++N_5 ;
 
      if( N_good + Ne_good < 4) continue ; 	
      ++N_2 ;  // fill counter
      N_2_w=N_2_w+newweight;

      int Zxx_tag = 0;    // 1: Zmumu  ,  2: Zee

      int i1 = -1; //index of the first lepton (from Z1)
      int j1 = -1; //index of the second lepton (from Z1)
      int pi1 = -1; 
      int pj1 = -1;
      
      bool has_FSR_Z1 = 0;
      TLorentzVector Lepton1,Lepton2,DiLepton,LeptonCorrection;

      for(int i = 0; i < N_good; ++i){
        for(int j = i + 1; j < N_good; ++j){
//	  if (fabs(RECOMU_SIP[iL[i]])>=4.) continue; // SIP cut
//	  if (fabs(RECOMU_SIP[iL[j]])>=4.) continue;
	  if (fabs(RECOMU_PFX_dB_new[iL[i]])>=0.20) continue; // Isolation
	  if (fabs(RECOMU_PFX_dB_new[iL[j]])>=0.20) continue;
	  
	  if(RECOMU_CHARGE[ iL[j] ] == RECOMU_CHARGE[ iL[i] ]) continue; // opposite charge
//same sign (qier)(don't forget to change it back)
//          if(RECOMU_CHARGE[ iL[j] ] != RECOMU_CHARGE[ iL[i] ]) continue;

	  if(jentry%5000 == 0) cout << "\n Pairing muons with pT= " << RECOMU_PT[ iL[i] ] << " and " <<  RECOMU_PT[ iL[j] ] << endl;
		  
	  // evaluate the mass &
	  double pxZ, pyZ, pzZ;
	  double EZ;
	  double massZ;
	  double massZ_noFSR = 0;
	  
	  int tempphotid=-1;
	  int templepid=-1;

	  float pTphot=-999.;
	  Lepton1.SetPtEtaPhiM(RECOMU_PT[iL[i]], RECOMU_ETA[iL[i]], RECOMU_PHI[iL[i]], 0.105);
	  Lepton2.SetPtEtaPhiM(RECOMU_PT[iL[j]], RECOMU_ETA[iL[j]], RECOMU_PHI[iL[j]], 0.105);
	  DiLepton=Lepton1+Lepton2;	  
	  massZ = DiLepton.M();	  
	  massZ_noFSR = massZ;
	  if (debug) cout << "Mass Z= " << massZ << endl;
	  pxZ=DiLepton.Px();
	  pyZ=DiLepton.Py();
	  pzZ=DiLepton.Pz();
	  EZ=DiLepton.E();	  
	 
	  Zxx_tag=1;	 

	  // ** Association of FSR to Z
	  if( debug ) cout  << "Step Z+FSR  " << endl;
	  
	  bool has_FSR_Z = 0;
	  int N_FSR_Z = 0;
	  double max_pt_FSR_Z = -1.;
	  int pi = -1; 
	  int pj = -1;
	  
	  
	  for( int p = 0; p < Nphotons; ++p ){
	    
	    if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )  {  // exists a photon associated to a lepton mu
	      if(debug) cout << "Attaching a photon to muon and then to the Z" << endl;
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton1=Lepton1+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();	      	    

	      //cout << mllp << " " << Zmass << " " << massZ << endl;
	      pi = p; 
	      has_FSR_Z = 1;
	      ++N_FSR_Z;
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      if(debug) cout << "Mass Z with FSR= "<< massZ << endl;

	    }
	    if( iLp_l[ p ] == iL[j] && iLp_tagEM[ p ] == 0 )  { 
	      if(debug) cout << "Attaching a photon to muon and then to the Z" << endl;
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton2=Lepton2+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();

	      //cout << mllp << " " << Zmass << " " << massZ << endl;
	      pj = p;
	      has_FSR_Z = 1;
	      ++N_FSR_Z; 
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      if(debug) cout << "Mass Z with FSR= "<< massZ << endl;

	    }
	  } // end loop on FSR photons

	 
	  
	  
	  if( debug && has_FSR_Z) {
	    cout  << " Z has FSR! " << endl;
	    cout  << "  N_FSR_Z " << N_FSR_Z << endl;
	    cout  << "  max_pt of photon FSR_Z " << max_pt_FSR_Z << endl;
	    if( pi > -1 ) cout  << "  pi " << pi << " --> index photon: " << iLp[pi] << " associated lepton: " << iLp_l[pi] << " (= "<< iL[i]<<" ? )  tag: " << iLp_tagEM[pi] << endl;
	    if( pj > -1 ) cout  << "  pj " << pj << " --> index photon: " << iLp[pj] << " associated lepton: " << iLp_l[pj] << " (= "<< iL[j]<<" ? )  tag: " << iLp_tagEM[pj] << endl;
	  }
	  else {
	    if(debug) cout << "No FSR photon attached" << endl;
	  }
	  
	  
	  if( has_FSR_Z ){ // if Z has FSR
	    
	    ++N_3_FSR; // fill the counter
	    N_3_FSR_w=N_3_FSR_w+newweight;
	    	    
	    // do not recompute isolation here
	    if( debug ) cout  << "Z Isolation (not corrected for photon): "
                              << "\n RECOMU_PFX_dB[ iL[i] ] " << RECOMU_PFX_dB[ iL[i] ]
                              << "\n RECOMU_PFX_dB[ iL[j] ] " << RECOMU_PFX_dB[ iL[j] ]
                              << endl;
	    	   
	    if( pi != -1 ){
	      //double deltaR_i = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pi]] , RECOMU_PHI[iL[i]] ),2) + pow(RECOPFPHOT_ETA[iLp[pi]] - RECOMU_ETA[iL[i]],2) );
              //if( deltaR_i < 0.4 && deltaR_i > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pi]];	      
	    }
	    else if( pj != -1 ){	      
              //double deltaR_j = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pj]] , RECOMU_PHI[iL[j]] ),2) + pow(RECOPFPHOT_ETA[iLp[pj]] - RECOMU_ETA[iL[j]],2) );
              //if( deltaR_j < 0.4 && deltaR_j > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pj]];	      
	    }

 
	  } // end if has FSR
	  else{
	    if( debug ) cout  << "Z Isolation: "  
			      << "\n RECOMU_PFX_dB_new[ iL[i] ] " << RECOMU_PFX_dB_new[ iL[i] ]
			      << "\n RECOMU_PFX_dB_new[ iL[j] ] " << RECOMU_PFX_dB_new[ iL[j] ]
			      << endl;	    
	  }
	  // ** end association of FSR to Z
	  
	  //if( massZ == 0 || i1 == -1 || j1 == -1) continue;
	  
	  if(jentry%5000 == 0) cout << "2mu: " << Zxx_tag << endl; 
	  
	  if(jentry%5000 == 0) cout << "Filling a struct for Z" << endl; 
	  candidateZ *Z = new candidateZ;
	  Z->massvalue=massZ;
	  Z->ilept1=iL[i];
	  Z->ilept2=iL[j];
	  Z->pt1=RECOMU_PT[iL[i]];
	  Z->pt2=RECOMU_PT[iL[j]];
	  Z->eta1=RECOMU_ETA[iL[i]];
	  Z->eta2=RECOMU_ETA[iL[j]];
	  Z->phi1=RECOMU_PHI[iL[i]];
	  Z->phi2=RECOMU_PHI[iL[j]];
	  Z->charge1=RECOMU_CHARGE[iL[i]];
	  Z->charge2=RECOMU_CHARGE[iL[j]];
	  Z->isol1=RECOMU_PFX_dB_new[ iL[i] ];
	  Z->isol2=RECOMU_PFX_dB_new[ iL[j] ];
	  if( pi != -1 ) Z->ilept1_FSR=true;
	  if( pj != -1 ) Z->ilept2_FSR=true;
	  Z->pxZ=pxZ;
	  Z->pyZ=pyZ;
	  Z->pzZ=pzZ;
	  Z->EZ=EZ;
	  if( has_FSR_Z ) {
	    Z->withFSR=1;
	    Z->ptFSR=pTphot;	    
	  }	      
	  else {
	    Z->withFSR=0;
	    Z->ptFSR=0.;
	  }
	  Z->tag=Zxx_tag;
	  	 	  
	  Zcandvector.push_back(*Z);	  
	  
	}
      } // end loop on pairs

      // 2mu2e
      for(int i = 0; i < Ne_good; ++i){
        for(int j = i + 1; j < Ne_good; ++j){
//	  if (fabs(RECOELE_SIP[iLe[i]])>=4.) continue; // SIP cut
//	  if (fabs(RECOELE_SIP[iLe[j]])>=4.) continue;
	  if (fabs(RECOELE_PFX_rho_new[iLe[i]])>=0.20) continue; // Isolation cut
	  if (fabs(RECOELE_PFX_rho_new[iLe[j]])>=0.20) continue;
	  
	  if(RECOELE_CHARGE[ iLe[j] ] == RECOELE_CHARGE[ iLe[i] ]) continue; // opposite charge

	  if(debug) cout << "\n Pairing electrons with pT= " << RECOELE_PT[ iLe[i] ] << " and " <<  RECOELE_PT[ iLe[j] ] << endl;
	  
	  // evaluate the mass &
	  double pxZ, pyZ, pzZ;
	  double EZ;
	  double massZ;
	  double massZ_noFSR = 0;
	  
	 
	  int tempphotid=-1;
	  int templepid=-1;
	  
	  float pTphot=-999.;
	  Lepton1.SetPtEtaPhiM(RECOELE_PT[iLe[i]], RECOELE_ETA[iLe[i]], RECOELE_PHI[iLe[i]], 0.000511);
	  Lepton2.SetPtEtaPhiM(RECOELE_PT[iLe[j]], RECOELE_ETA[iLe[j]], RECOELE_PHI[iLe[j]], 0.000511);
	  DiLepton=Lepton1+Lepton2;	  
	  massZ = DiLepton.M();	  
	  massZ_noFSR = massZ;
	  if (debug) cout << "Mass Z= " << massZ << endl;
	  pxZ=DiLepton.Px();
	  pyZ=DiLepton.Py();
	  pzZ=DiLepton.Pz();
	  EZ=DiLepton.E();
	  
	  Zxx_tag=2;
	
	  // ** Association of FSR to Z
	  if( debug ) cout  << "Step Z+FSR  " << endl;
	  
	  bool has_FSR_Z = 0;
	  int N_FSR_Z = 0;
	  double max_pt_FSR_Z = -1.;
	  int pi = -1; 
	  int pj = -1;
	  
	  
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {  // exit a photon associated to a lepton electron
	      
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton1=Lepton1+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();	      	      	   
	      
	      has_FSR_Z = 1; 
	      pi = p; 
	      ++N_FSR_Z;
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      if(debug) cout << "Mass Z with FSR= "<< massZ << endl;
		      
	    }
	    
	    if( iLp_l[ p ] == iLe[j] && iLp_tagEM[ p ] == 1 )  { 
	      
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton2=Lepton2+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();
	      	    
	      pj = p;
	      has_FSR_Z = 1;
	      ++N_FSR_Z; 
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      if(debug) cout << "Mass Z with FSR= "<< massZ << endl;

	    }
	  } // end loop on FSR photons
	  
	  
	  
	  //if( has_FSR_Z ) debug = 1;
	  
	  if( debug && has_FSR_Z) {
	    cout  << " Z has FSR! " << endl;
	    cout  << "  N_FSR_Z " << N_FSR_Z << endl;
	    cout  << "  max_pt of photon FSR_Z " << max_pt_FSR_Z << endl;
	    if( pi > -1 ) cout  << "  pi " << pi << " --> index photon: " << iLp[pi] << " associated lepton: " << iLp_l[pi] << " (= "<< iLe[i]<<" ? )  tag: " << iLp_tagEM[pi] << endl;
	    if( pj > -1 ) cout  << "  pj " << pj << " --> index photon: " << iLp[pj] << " associated lepton: " << iLp_l[pj] << " (= "<< iLe[j]<<" ? )  tag: " << iLp_tagEM[pj] << endl;
	  }
	  else {
	    if(debug)  cout << "No FSR photon attached" << endl;
	  }
	  
	  
	  if( has_FSR_Z ){ // if Z has FSR
	    
	    ++N_3_FSR; // fill the counter
	    N_3_FSR_w=N_3_FSR_w+newweight;	    
	    
	    // do not recompute isolation here
	    if( debug ) cout  << "Z Isolation (not corrected for photon): "
                              << "\n RECOELE_PFX_rho[ iLe[i] ] " << RECOELE_PFX_rho[ iLe[i] ]
                              << "\n RECOELE_PFX_rho[ iLe[j] ] " << RECOELE_PFX_rho[ iLe[j] ]
                              << endl;
	    
	    if( pi != -1 ){
	      //double deltaR_i = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pi]] , RECOELE_PHI[iLe[i]] ),2) + pow(RECOPFPHOT_ETA[iLp[pi]] - RECOELE_ETA[iLe[i]],2) );
              //if( deltaR_i < 0.4 && deltaR_i > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pi]];	      
	    }
	    else if( pj != -1 ){	      
              //double deltaR_j = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pj]] , RECOELE_PHI[iLe[j]] ),2) + pow(RECOPFPHOT_ETA[iLp[pj]] - RECOELE_ETA[iLe[j]],2) );
              //if( deltaR_j < 0.4 && deltaR_j > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pj]];	      
	    }
	    
	    
	  } // end if has FSR
	  else{
	    
	    if( debug ) cout  << "Z Isolation: "  
			      << "\n RECOELE_PFX_rho_new[ iLe[i] ] " << RECOELE_PFX_rho_new[ iLe[i] ]
			      << "\n RECOELE_PFX_rho_new[ iLe[j] ] " << RECOELE_PFX_rho_new[ iLe[j] ]
			      << endl;	    
	  }
	  // ** end association of FSR to Z
	  
	  //if( massZ == 0 || i1 == -1 || j1 == -1) continue;
	  
	  if(jentry%5000 == 0) cout << "2e2mu: " << Zxx_tag << endl;
	   
	  if(jentry%5000 == 0) cout << "Filling a struct for Z" << endl; 
	  candidateZ *Z = new candidateZ;
	  Z->massvalue=massZ;
	  Z->ilept1=iLe[i];
	  Z->ilept2=iLe[j];
	  Z->pt1=RECOELE_PT[iLe[i]];
	  Z->pt2=RECOELE_PT[iLe[j]];
	  Z->eta1=RECOELE_ETA[iLe[i]];
	  Z->eta2=RECOELE_ETA[iLe[j]];
	  Z->phi1=RECOELE_PHI[iLe[i]];
	  Z->phi2=RECOELE_PHI[iLe[j]];
	  Z->charge1=RECOELE_CHARGE[iLe[i]];
	  Z->charge2=RECOELE_CHARGE[iLe[j]];
	  Z->isol1=RECOELE_PFX_rho_new[ iLe[i] ];
	  Z->isol2=RECOELE_PFX_rho_new[ iLe[j] ];
	  if( pi != -1 ) Z->ilept1_FSR=true;
	  if( pj != -1 ) Z->ilept2_FSR=true;
	  Z->pxZ=pxZ;
	  Z->pyZ=pyZ;
	  Z->pzZ=pzZ;
	  Z->EZ=EZ;
	  if( has_FSR_Z ) {
	    Z->withFSR=1;
	    Z->ptFSR=pTphot;
	  }	      
	  else {
	    Z->withFSR=0;
	    Z->ptFSR=0.;
	  }
	  	 
	  Z->tag=Zxx_tag;
	  Zcandvector.push_back(*Z);	  
	  
	}
      } // end loop on couples
      
      
      if (Zcandvector.size()<2) {
	if(debug) cout << "Less than two Z pairs with isolated leptons...exiting" << endl;
	continue; 
      }
      

      ++N_3a ;  // fill counter
      N_3a_w=N_3a_w+newweight;
      

      // Mass cut on Z
      vector<candidateZ> Zcandisolmassvector;
      Zcandisolmassvector.clear();

      for (int index=0; index<Zcandvector.size();index++){
//	if (!(Zcandvector.at(index).massvalue > 60 && Zcandvector.at(index).massvalue < 120)) continue;
	if(debug) cout << "Z passing the 60 < mll < 120 cut with mass= " << Zcandvector.at(index).massvalue<< endl;
	Zcandisolmassvector.push_back(Zcandvector.at(index));
      };
      
      if (Zcandisolmassvector.size()<2) {
	if(debug) cout << "No Z passing the mass cut"<< endl;
	continue;
      }

      if(jentry%5000 == 0) cout << "Number of Z passing the isolation and the 60 < mll < 120 cut is= " << Zcandisolmassvector.size() << endl;

      ++N_3b ;  // fill counter
      N_3b_w=N_3b_w+newweight;


      // **** Step 4:
       // a) 4 leptons
      // b) pair #2
      // c) highest pt
      // d) mZ2 in ] 4,120 [

      int issamesign = 0;

      //if( debug ) cout  << "\nStep 4: Number of good leptons: " << N_good << endl;

      int N_Z2_pairs = 0;

      int i2 = -1; //index of the first lepton (from Z1)
      int j2 = -1; //index of the second lepton (from Z1)
      int pi2 = -1; 
      int pj2 = -1; 
      
      bool has_FSR_Z2 = 0;
      
      
      // PT,20/10 for any di-lepton
/*
      vector<candidateZ> firstpTcleanedgoodZ;    
      vector<float> leptonspTcleaned;
      
      for (int l=0;l<Zcandisolmassvector.size();l++){
        leptonspTcleaned.clear();
        leptonspTcleaned.push_back((Zcandisolmassvector.at(l)).pt1);
        leptonspTcleaned.push_back((Zcandisolmassvector.at(l)).pt2);
        std::sort(leptonspTcleaned.rbegin(),leptonspTcleaned.rend());
      }
*/
      vector<candidateZ> pTcleanedgoodZ;
      pTcleanedgoodZ=Zcandisolmassvector;
            

      if (pTcleanedgoodZ.size()<2) {
        if(debug) cout << "No Z (PT required) passing the mass cut"<< endl;
        continue;
      }     
 
      // Z1 selection
      double pxZ1 = 0;  //Z1 kinematics
      double pyZ1 = 0;
      double pzZ1 = 0;
      double ptZ1 = 0;
      double EZ1 = 0;
      double Y_Z1 = -9;
      double massZ1 = 0;
      double massZ1_noFSR = 0;
      double sum_ptZ1 = 0.;
      int indexlep1Z1 = -1;
      int indexlep2Z1 = -1;
      int indexZ1= -1;
      int Z1tag=-999;

      // Z2 selection
      double pxZ2 = 0;  //Z2 kinematics
      double pyZ2 = 0;
      double pzZ2 = 0;
      double ptZ2 = 0;
      double EZ2 = 0;
      double Y_Z2 = -9;
      double massZ2 = 0;
      double massZ2_noFSR = 0;
      double sum_ptZ2 = 0.;
      int indexlep1Z2 = -1;
      int indexlep2Z2 = -1;
      int indexZ2= -1;
      int Z2tag=-999;
      int Z1charge=0;
      int Z2charge=0;
      
      // Choice of Z1 and Z2 as the closest to the 2Z mass using minimum Ki
      for (int i=0;i<pTcleanedgoodZ.size();++i){
        for(int j=i+1; j<pTcleanedgoodZ.size();++j){
         //skip Z that contains same lepton
         if((pTcleanedgoodZ.at(i).tag==pTcleanedgoodZ.at(j).tag)&&
            (pTcleanedgoodZ.at(i).ilept1==pTcleanedgoodZ.at(j).ilept1
             ||pTcleanedgoodZ.at(i).ilept1==pTcleanedgoodZ.at(j).ilept2
             ||pTcleanedgoodZ.at(i).ilept2==pTcleanedgoodZ.at(j).ilept1
             ||pTcleanedgoodZ.at(i).ilept2==pTcleanedgoodZ.at(j).ilept2)) continue;
	
        // ki square minimum
	if( pow((pTcleanedgoodZ.at(i).massvalue - Zmass),2)+pow((pTcleanedgoodZ.at(j).massvalue - Zmass),2) < pow((massZ1 - Zmass),2)+pow((massZ2 - Zmass),2)){
	  int a,b;
          double zpti=sqrt(pTcleanedgoodZ.at(i).pxZ*pTcleanedgoodZ.at(i).pxZ+pTcleanedgoodZ.at(i).pyZ*pTcleanedgoodZ.at(i).pyZ);
          double zptj=sqrt(pTcleanedgoodZ.at(j).pxZ*pTcleanedgoodZ.at(j).pxZ+pTcleanedgoodZ.at(j).pyZ*pTcleanedgoodZ.at(j).pyZ);
          if(zpti > zptj){ a=i; b=j;}
          else a=j,b=i;

	  massZ1 = pTcleanedgoodZ.at(a).massvalue;
          massZ2 = pTcleanedgoodZ.at(b).massvalue;
	  indexZ1=a;
          indexZ2=b;
	  
	  pxZ1 = pTcleanedgoodZ.at(a).pxZ;
	  pyZ1 = pTcleanedgoodZ.at(a).pyZ;
	  pzZ1 = pTcleanedgoodZ.at(a).pzZ;
	  EZ1  = pTcleanedgoodZ.at(a).EZ;

          pxZ2 = pTcleanedgoodZ.at(b).pxZ;
          pyZ2 = pTcleanedgoodZ.at(b).pyZ;
          pzZ2 = pTcleanedgoodZ.at(b).pzZ;
          EZ2  = pTcleanedgoodZ.at(b).EZ;

	  
	  ptZ1 = sqrt( pxZ1*pxZ1 + pyZ1*pyZ1 );
	  sum_ptZ1 = pTcleanedgoodZ.at(a).pt1+pTcleanedgoodZ.at(a).pt2;

          ptZ2 = sqrt( pxZ2*pxZ2 + pyZ2*pyZ2 );
          sum_ptZ2 = pTcleanedgoodZ.at(b).pt1+pTcleanedgoodZ.at(b).pt2;
	  
	  Y_Z1 = 0.5 * log ( (EZ1 + pzZ1)/(EZ1 - pzZ1) );
	  indexlep1Z1=pTcleanedgoodZ.at(a).ilept1;
	  indexlep2Z1=pTcleanedgoodZ.at(a).ilept2;
	  Z1tag=pTcleanedgoodZ.at(a).tag;
          if(pTcleanedgoodZ.at(a).charge1!=pTcleanedgoodZ.at(a).charge2) Z1charge=-1;
          else Z1charge=1;

          Y_Z2 = 0.5 * log ( (EZ2 + pzZ2)/(EZ2 - pzZ2) );
          indexlep1Z2=pTcleanedgoodZ.at(b).ilept1;
          indexlep2Z2=pTcleanedgoodZ.at(b).ilept2;
          Z2tag=pTcleanedgoodZ.at(b).tag;
          if(pTcleanedgoodZ.at(b).charge1!=pTcleanedgoodZ.at(b).charge2) Z2charge=-1;
          else Z2charge=1;
	}
      }
      }
     
      
 
/*      if (massZ1 < 40.&& massZ2 < 40. ) {
	if(debug) cout << "The mass of Z1 is < 40 GeV...exiting" << endl;
	continue;
      } 
     
*/
      if(!((massZ1>60&&massZ1<120)||(massZ2>60&&massZ2<120))) continue;

      if( debug ) cout  << "\n Final Z1 properties: "
			<< "\n pxZ1 " << pxZ1
			<< "\n pyZ1 " << pyZ1
			<< "\n pzZ1 " << pzZ1
			<< "\n ptZ1 " << ptZ1
			<< "\n EZ1 "  << EZ1
			<< "\n Y_Z1 " << Y_Z1
			<< "\n massZ1 " << massZ1
			<< "\n indexlep1 " << indexlep1Z1
			<< "\n indexlep2 " << indexlep2Z1
			<< "\n indexZ1 " << indexZ1 
			<< "\n Z1 tag (1 for 2mu and 2 for 2e) " << Z1tag
		    
			<< endl;

      if( debug ) cout  << "\n Final Z2 properties: "
                        << "\n pxZ2 " << pxZ2
                        << "\n pyZ2 " << pyZ2
                        << "\n pzZ2 " << pzZ2
                        << "\n ptZ2 " << ptZ2
                        << "\n EZ2 "  << EZ2
                        << "\n Y_Z2 " << Y_Z2
                        << "\n massZ2 " << massZ2
                        << "\n indexlep1 " << indexlep1Z2
                        << "\n indexlep2 " << indexlep2Z2
                        << "\n indexZ2 " << indexZ2
                        << "\n Z2 tag (1 for 2mu and 2 for 2e) " << Z2tag

                        << endl;

      
      
      ++N_4b ;  // fill counter
      N_4b_w=N_4b_w+newweight;

      
      // **** Step 5:

      
       // Execute Efficiency Reweighting
 
      int z1lept[2]={indexlep1Z1,indexlep2Z1};
      int z2lept[2]={indexlep1Z2,indexlep2Z2};

      Double_t eff_weight = 1.;
      
      if (Z1tag==1){ 	
        for(int i = 0; i < 2; ++i){
	  Double_t Pt = RECOMU_PT[ z1lept[i] ]; 
	  Double_t Eta = RECOMU_ETA[ z1lept[i] ]; 
	  
	  if( MC_type == "Spring16" && DATA_type == "NO"){
            int biny1 = mu_scale_factors_id_p1->GetYaxis()->FindBin(Pt);
            int binx1 = mu_scale_factors_id_p1->GetXaxis()->FindBin(abs(Eta));
  
            double sf_id=mu_scale_factors_id_p1->GetBinContent(binx1,biny1);
            if (sf_id>0.) eff_weight*=sf_id;
  
            int biny2 = mu_scale_factors_iso_p1->GetYaxis()->FindBin(Pt);
            int binx2 = mu_scale_factors_iso_p1->GetXaxis()->FindBin(abs(Eta));
  
            double sf_iso = mu_scale_factors_iso_p1->GetBinContent(binx2,biny2);
            if (sf_iso>0.)  eff_weight*=sf_iso;
  
//            double tk_sf = mu_scale_factors_tk->Eval(Eta);
//            if(mu_scale_factors_tk->Eval(Eta)>0) eff_weight*=tk_sf;
  
//            cout << " id weight = " << sf_id << "\n iso weight = " << sf_iso << "\n tk weight = " << tk_sf << endl;
          }
	}
      }
      else if (Z1tag==2){ 
	for(int i = 0; i < 2; ++i){
	  Double_t Pt = RECOELE_PT[ z1lept[i] ]; 
	  Double_t Eta = RECOELE_ETA[ z1lept[i] ]; 
	  
	  if( MC_type == "Spring16" && DATA_type == "NO"){
            if(debug) cout << "Pt= " << Pt << " Eta= " << Eta << endl;
            int biny4 = ele_scale_factors_reco->GetYaxis()->FindBin(Pt);
            int binx4 = ele_scale_factors_reco->GetXaxis()->FindBin(Eta);
            if (ele_scale_factors_reco->GetBinContent(binx4,biny4)>0.) eff_weight*=ele_scale_factors_reco->GetBinContent(binx4,biny4);
            if(debug) cout << "ele reco sf = " << ele_scale_factors_reco->GetBinContent(binx4,biny4) << endl;

            int biny5 = ele_scale_factors_wp90->GetYaxis()->FindBin(Pt);
            int binx5 = ele_scale_factors_wp90->GetXaxis()->FindBin(Eta);
            if (ele_scale_factors_wp90->GetBinContent(binx5,biny5)>0.) eff_weight*=ele_scale_factors_wp90->GetBinContent(binx5,biny5);
            if(debug) cout << "ele wp90 sf = " << ele_scale_factors_wp90->GetBinContent(binx5,biny5) << endl;


          }
        }
	    
      }

      if (Z2tag==1){
        for(int i = 0; i < 2; ++i){
          Double_t Pt = RECOMU_PT[ z2lept[i] ];
          Double_t Eta = RECOMU_ETA[ z2lept[i] ];
          if( MC_type == "Spring16" && DATA_type == "NO"){
          int biny1 = mu_scale_factors_id_p1->GetYaxis()->FindBin(Pt);
          int binx1 = mu_scale_factors_id_p1->GetXaxis()->FindBin(abs(Eta));

          double sf_id=mu_scale_factors_id_p1->GetBinContent(binx1,biny1);
          if (sf_id>0.) eff_weight*=sf_id;

          int biny2 = mu_scale_factors_iso_p1->GetYaxis()->FindBin(Pt);
          int binx2 = mu_scale_factors_iso_p1->GetXaxis()->FindBin(abs(Eta));

          double sf_iso = mu_scale_factors_iso_p1->GetBinContent(binx2,biny2);
          if (sf_iso>0.)  eff_weight*=sf_iso;

//          double tk_sf = mu_scale_factors_tk->Eval(Eta);
//          if(mu_scale_factors_tk->Eval(Eta)>0) eff_weight*=tk_sf;

//          cout << " id weight = " << sf_id << "\n iso weight = " << sf_iso << "\n tk weight = " << tk_sf << endl;

         }
        }
      }
      else if (Z2tag==2){
        for(int i = 0; i < 2; ++i){
          Double_t Pt = RECOELE_PT[ z2lept[i] ];
          Double_t Eta = RECOELE_ETA[ z2lept[i] ];

          if( MC_type == "Spring16" && DATA_type == "NO"){
            if(debug) cout << "Pt= " << Pt << " Eta= " << Eta << endl;
            int biny4 = ele_scale_factors_reco->GetYaxis()->FindBin(Pt);
            int binx4 = ele_scale_factors_reco->GetXaxis()->FindBin(Eta);
            if (ele_scale_factors_reco->GetBinContent(binx4,biny4)>0.) eff_weight*=ele_scale_factors_reco->GetBinContent(binx4,biny4); 
            if(debug) cout << "ele reco sf = " << ele_scale_factors_reco->GetBinContent(binx4,biny4) << endl; 
            
            int biny5 = ele_scale_factors_wp90->GetYaxis()->FindBin(Pt);
            int binx5 = ele_scale_factors_wp90->GetXaxis()->FindBin(Eta);
            if (ele_scale_factors_wp90->GetBinContent(binx5,biny5)>0.) eff_weight*=ele_scale_factors_wp90->GetBinContent(binx5,biny5);
            if(debug) cout << "ele wp90 sf = " << ele_scale_factors_wp90->GetBinContent(binx5,biny5) << endl; 

          }
        }
      }

//trigger SF 
      //finding leading two muons
       double leadpt=-1;
       double subpt=-1;
       double leadeta,subeta;
       int indexleptonfinal[4]={indexlep1Z1,indexlep2Z1,indexlep1Z2,indexlep2Z2};

      int indextwomu[2]={-1,-1};
      int indextwoele[2]={-1,-1};
      if(Z1tag==1&&Z2tag==1){
        for(int i=0; i<4; i++){
          if(RECOMU_PT[indexleptonfinal[i]]>leadpt) {leadpt=RECOMU_PT[indexleptonfinal[i]]; leadeta=RECOMU_ETA[indexleptonfinal[i]];indextwomu[0]=indexleptonfinal[i];}
        }
        for(int i=0; i<4; i++){
          if(RECOMU_PT[indexleptonfinal[i]]==leadpt) continue;
          if(RECOMU_PT[indexleptonfinal[i]]>subpt) {subpt=RECOMU_PT[indexleptonfinal[i]];subeta= RECOMU_ETA[indexleptonfinal[i]];indextwomu[1]=indexleptonfinal[i];}
        }
      }
      else if(Z1tag==1&&Z2tag!=1){
        indextwomu[0]=indexleptonfinal[0];
        indextwomu[1]=indexleptonfinal[1];
        indextwoele[0]=indexleptonfinal[2];
        indextwoele[1]=indexleptonfinal[3];
      }
      else if(Z1tag!=1&&Z2tag==1){
        indextwomu[0]=indexleptonfinal[2];
        indextwomu[1]=indexleptonfinal[3];
        indextwoele[0]=indexleptonfinal[0];
        indextwoele[1]=indexleptonfinal[1];
      }
     else if(Z1tag==2&&Z2tag==2){
        for(int i=0; i<4; i++){
          if(RECOELE_PT[indexleptonfinal[i]]>leadpt) {leadpt=RECOELE_PT[indexleptonfinal[i]]; leadeta=RECOELE_ETA[indexleptonfinal[i]];indextwoele[0]=indexleptonfinal[i];}
        }
        for(int i=0; i<4; i++){
          if(RECOELE_PT[indexleptonfinal[i]]==leadpt) continue;
          if(RECOELE_PT[indexleptonfinal[i]]>subpt) {subpt=RECOELE_PT[indexleptonfinal[i]];subeta= RECOELE_ETA[indexleptonfinal[i]];indextwoele[1]=indexleptonfinal[i];}
        }
      }


     bool nomatch=false;
     bool noleg1=true;
//dm_trig matching (qier)

    if(dm_trig){ 
    for(int i=0; i<2; i++){
       if( MC_type == "Spring16" && DATA_type == "NO"){
          Double_t Pt = RECOMU_PT[indextwomu[i]];
          Double_t Eta = RECOMU_ETA[indextwomu[i]];
          if(RECOMU_dm_MuHLTMatch[indextwomu[i]]== 2){
            int biny3 = mu_scale_factors_hlt_p1->GetYaxis()->FindBin(Pt);
            int binx3 = mu_scale_factors_hlt_p1->GetXaxis()->FindBin(abs(Eta));

//            int biny32 = mu_scale_factors_hlt_p2->GetYaxis()->FindBin(Pt);
//            int binx32 = mu_scale_factors_hlt_p2->GetXaxis()->FindBin(abs(Eta));
//            double sf_hlt = mu_scale_factors_hlt_p1->GetBinContent(binx3,biny3)*19.666/35.812+mu_scale_factors_hlt_p2->GetBinContent(binx32,biny32)*16.146/35.812;
            double sf_hlt = mu_scale_factors_hlt_p1->GetBinContent(binx3,biny3);
            if (sf_hlt>0.) eff_weight*=sf_hlt;
            if(debug) cout << "l1trigger matching mu17 leg weight = " << sf_hlt << endl;
           }
         else{
            int biny32 = mu_scale_factors_hlt_p2->GetYaxis()->FindBin(Pt);
            int binx32 = mu_scale_factors_hlt_p2->GetXaxis()->FindBin(abs(Eta));

            double sf_hlt = mu_scale_factors_hlt_p2->GetBinContent(binx32,biny32);
            if (sf_hlt>0.) eff_weight*=sf_hlt;
            if(debug) cout << "l1trigger matching mu8 leg weight = " << sf_hlt << endl;
          }
       }
     if(RECOMU_dm_MuHLTMatch[indextwomu[i]]<0) nomatch=true;
     if(RECOMU_dm_MuHLTMatch[indextwomu[i]]==2) noleg1=false;
     }
   }


//de_trig matching(qier)

   if(de_trig){
    for(int i=0; i<2; i++){
       if( MC_type == "Spring16" && DATA_type == "NO"){
          Double_t Pt = RECOELE_PT[indextwoele[i]];
          Double_t Eta = RECOELE_ETA[indextwoele[i]];
          if(RECOELE_de_EleHLTMatch[indextwoele[i]]== 2){
            int biny4 = ele_scale_factors_leg1->GetYaxis()->FindBin(Pt);
            int binx4 = ele_scale_factors_leg1->GetXaxis()->FindBin(Eta);

            double sf_hlt = ele_scale_factors_leg1->GetBinContent(binx4,biny4);
            if (sf_hlt>0.) eff_weight*=sf_hlt;
            if(debug) cout << "hlt matching ele leg1 weight = " << sf_hlt << endl;
           }
         else{
            int biny42 = ele_scale_factors_leg2->GetYaxis()->FindBin(Pt);
            int binx42 = ele_scale_factors_leg2->GetXaxis()->FindBin(Eta);

            double sf_hlt = ele_scale_factors_leg2->GetBinContent(binx42,biny42);
            if (sf_hlt>0.) eff_weight*=sf_hlt;
            if(debug) cout << "hlt matching ele leg2 weight = " << sf_hlt << endl;
          }
       }
     if(RECOELE_de_EleHLTMatch[indextwoele[i]]<0) nomatch=true; 
     if(RECOELE_de_EleHLTMatch[indextwoele[i]]==2) noleg1=false;
     }
  }

     if(nomatch||noleg1) continue;
     if(debug) cout << "eff_weight" << eff_weight << endl;
      
      // // Changing the weight for pileup and efficiency
      if (eff_weight>0.) newweight=newweight*eff_weight;
      
      if(jentry%5000 == 0) cout << "Starting weight + pileup + efficiency= " << newweight << endl;
      if(debug) cout << "Efficiency Weight for the 4l: " << eff_weight << " Final weight= " << newweight << endl;

      // sort index by pt (kinematics not corrected for FSR)
      int ipt[4] ;
      double tmp_pt[4];
      int tmp_type[4];
      int lep_type[4];

      //cout << "PTs= " << RECOMU_PT[indexleptonfinal[0]] << " " << RECOMU_PT[indexleptonfinal[1]] << " " <<  RECOMU_PT[indexleptonfinal[2]] << " " << RECOMU_PT[indexleptonfinal[3]]<< endl;

      for(int i = 0; i < 2; ++i){
        if (Z1tag==1) {tmp_pt[i] =  RECOMU_PT[indexleptonfinal[i]]; tmp_type[i]=1;}
        if (Z1tag==2) {tmp_pt[i] =  RECOELE_PT[indexleptonfinal[i]]; tmp_type[i]=2;}
        if(debug) cout << tmp_pt[i] << endl;
      }
      for(int i = 2; i < 4; ++i){
        if (Z2tag==1) {tmp_pt[i] =  RECOMU_PT[indexleptonfinal[i]]; tmp_type[i]=1;}
        if (Z2tag==2) {tmp_pt[i] =  RECOELE_PT[indexleptonfinal[i]]; tmp_type[i]=2;}
        if(debug) cout << tmp_pt[i] << endl;
      }

      
      float sortedpT[4];

      for(int i = 0; i < 4; ++i){		
        double tmp_max_pt = 0;
      	int jj = i;
        int type = 0;
        for(int j = 0; j < 4; ++j){
      		if( tmp_pt[j] > tmp_max_pt ){
      			tmp_max_pt = tmp_pt[j];
      			jj  = j;
                        type = tmp_type[j];
      		}
      	}
	sortedpT[i]=tmp_max_pt;
      	ipt[i] = indexleptonfinal[jj];
        lep_type[i]=type;
      	tmp_pt[jj] = 0;	
      }
      //end sorting
     
      //4 lepton pt selection
      if(debug) cout << "lepton Pt= "<< sortedpT[0] <<" "<< sortedpT[1] <<" " << sortedpT[2] <<" " << sortedpT[3] << endl;
      if(!(sortedpT[0]>30 && sortedpT[1]>20 && sortedpT[2] >10 && sortedpT[3]>10)) continue;
//      if(!(sortedpT[0]>20 && sortedpT[1]>10 && sortedpT[2] >5 && sortedpT[3]>5)) continue;
      //4mu selection qier
//      if(Z1tag!=1||Z2tag!=1) continue;
     if(!((Z1tag==1&&Z2tag==2)||(Z1tag==2&&Z2tag==1))) continue;
/*      int os=0;
      int ss=0;
      if(Z1charge==1) ss++;
      else os++;
      if(Z2charge==1) ss++;
      else os++;

      if(ss!=1) continue;
 */
      TLorentzVector Z1P4,Z2P4;
      Z1P4.SetPxPyPzE(pxZ1,pyZ1,pzZ1,EZ1);
      Z2P4.SetPxPyPzE(pxZ2,pyZ2,pzZ2,EZ2);

      hMZ1_5->Fill( massZ1,newweight );
      hPtZ1_5->Fill( ptZ1,newweight );
      hYZ1_5->Fill( Y_Z1,newweight );

      hMZ1_5_pdf_up->Fill( massZ1,newweight*(1+PDF_weighting_un));
      hPtZ1_5_pdf_up->Fill( ptZ1,newweight*(1+PDF_weighting_un));
      hYZ1_5_pdf_up->Fill( Y_Z1,newweight*(1+PDF_weighting_un));

      hMZ1_5_pdf_dow->Fill( massZ1,newweight*(1-PDF_weighting_un));
      hPtZ1_5_pdf_dow->Fill( ptZ1,newweight*(1-PDF_weighting_un));
      hYZ1_5_pdf_dow->Fill( Y_Z1,newweight*(1-PDF_weighting_un));

      hMZ1_5_pu_up->Fill( massZ1,newweight*pu_up);
      hPtZ1_5_pu_up->Fill( ptZ1,newweight*pu_up);
      hYZ1_5_pu_up->Fill( Y_Z1,newweight*pu_up);

      hMZ1_5_pu_dow->Fill( massZ1,newweight*pu_dow);
      hPtZ1_5_pu_dow->Fill( ptZ1,newweight*pu_dow);
      hYZ1_5_pu_dow->Fill( Y_Z1,newweight*pu_dow);

      hMZ1_5_mc_1->Fill( massZ1,newweight*mc_weight_un[1] );
      hPtZ1_5_mc_1->Fill( ptZ1,newweight*mc_weight_un[1]);
      hYZ1_5_mc_1->Fill( Y_Z1,newweight*mc_weight_un[1]);

      hMZ1_5_mc_2->Fill( massZ1,newweight*mc_weight_un[2] );
      hPtZ1_5_mc_2->Fill( ptZ1,newweight*mc_weight_un[2]);
      hYZ1_5_mc_2->Fill( Y_Z1,newweight*mc_weight_un[2]);

      hMZ1_5_mc_3->Fill( massZ1,newweight*mc_weight_un[3] );
      hPtZ1_5_mc_3->Fill( ptZ1,newweight*mc_weight_un[3]);
      hYZ1_5_mc_3->Fill( Y_Z1,newweight*mc_weight_un[3]);

      hMZ1_5_mc_4->Fill( massZ1,newweight*mc_weight_un[4] );
      hPtZ1_5_mc_4->Fill( ptZ1,newweight*mc_weight_un[4]);
      hYZ1_5_mc_4->Fill( Y_Z1,newweight*mc_weight_un[4]);

      hMZ1_5_mc_5->Fill( massZ1,newweight*mc_weight_un[5] );
      hPtZ1_5_mc_5->Fill( ptZ1,newweight*mc_weight_un[5]);
      hYZ1_5_mc_5->Fill( Y_Z1,newweight*mc_weight_un[5]);

      hMZ1_5_mc_6->Fill( massZ1,newweight*mc_weight_un[6] );
      hPtZ1_5_mc_6->Fill( ptZ1,newweight*mc_weight_un[6]);
      hYZ1_5_mc_6->Fill( Y_Z1,newweight*mc_weight_un[6]);

      hMZ1_5_mc_7->Fill( massZ1,newweight*mc_weight_un[7] );
      hPtZ1_5_mc_7->Fill( ptZ1,newweight*mc_weight_un[7]);
      hYZ1_5_mc_7->Fill( Y_Z1,newweight*mc_weight_un[7]);

      hMZ1_5_mc_8->Fill( massZ1,newweight*mc_weight_un[8] );
      hPtZ1_5_mc_8->Fill( ptZ1,newweight*mc_weight_un[8]);
      hYZ1_5_mc_8->Fill( Y_Z1,newweight*mc_weight_un[8]);

      hMZ2_5->Fill( massZ2,newweight );
      hPtZ2_5->Fill( ptZ2,newweight );
      hYZ2_5->Fill( Y_Z2,newweight );

      hMZ2_5_pdf_up->Fill( massZ2,newweight*(1+PDF_weighting_un));
      hPtZ2_5_pdf_up->Fill( ptZ2,newweight*(1+PDF_weighting_un));
      hYZ2_5_pdf_up->Fill( Y_Z2,newweight*(1+PDF_weighting_un));

      hMZ2_5_pdf_dow->Fill( massZ2,newweight*(1-PDF_weighting_un));
      hPtZ2_5_pdf_dow->Fill( ptZ2,newweight*(1-PDF_weighting_un));
      hYZ2_5_pdf_dow->Fill( Y_Z2,newweight*(1-PDF_weighting_un));

      hMZ2_5_pu_up->Fill( massZ2,newweight*pu_up);
      hPtZ2_5_pu_up->Fill( ptZ2,newweight*pu_up);
      hYZ2_5_pu_up->Fill( Y_Z2,newweight*pu_up);

      hMZ2_5_pu_dow->Fill( massZ2,newweight*pu_dow);
      hPtZ2_5_pu_dow->Fill( ptZ2,newweight*pu_dow);
      hYZ2_5_pu_dow->Fill( Y_Z2,newweight*pu_dow);

      hMZ2_5_mc_1->Fill( massZ2,newweight*mc_weight_un[1] );
      hPtZ2_5_mc_1->Fill( ptZ2,newweight*mc_weight_un[1]);
      hYZ2_5_mc_1->Fill( Y_Z2,newweight*mc_weight_un[1]);

      hMZ2_5_mc_2->Fill( massZ2,newweight*mc_weight_un[2] );
      hPtZ2_5_mc_2->Fill( ptZ2,newweight*mc_weight_un[2]);
      hYZ2_5_mc_2->Fill( Y_Z2,newweight*mc_weight_un[2]);

      hMZ2_5_mc_3->Fill( massZ2,newweight*mc_weight_un[3] );
      hPtZ2_5_mc_3->Fill( ptZ2,newweight*mc_weight_un[3]);
      hYZ2_5_mc_3->Fill( Y_Z2,newweight*mc_weight_un[3]);

      hMZ2_5_mc_4->Fill( massZ2,newweight*mc_weight_un[4] );
      hPtZ2_5_mc_4->Fill( ptZ2,newweight*mc_weight_un[4]);
      hYZ2_5_mc_4->Fill( Y_Z2,newweight*mc_weight_un[4]);

      hMZ2_5_mc_5->Fill( massZ2,newweight*mc_weight_un[5] );
      hPtZ2_5_mc_5->Fill( ptZ2,newweight*mc_weight_un[5]);
      hYZ2_5_mc_5->Fill( Y_Z2,newweight*mc_weight_un[5]);

      hMZ2_5_mc_6->Fill( massZ2,newweight*mc_weight_un[6] );
      hPtZ2_5_mc_6->Fill( ptZ2,newweight*mc_weight_un[6]);
      hYZ2_5_mc_6->Fill( Y_Z2,newweight*mc_weight_un[6]);

      hMZ2_5_mc_7->Fill( massZ2,newweight*mc_weight_un[7] );
      hPtZ2_5_mc_7->Fill( ptZ2,newweight*mc_weight_un[7]);
      hYZ2_5_mc_7->Fill( Y_Z2,newweight*mc_weight_un[7]);

      hMZ2_5_mc_8->Fill( massZ2,newweight*mc_weight_un[8] );
      hPtZ2_5_mc_8->Fill( ptZ2,newweight*mc_weight_un[8]);
      hYZ2_5_mc_8->Fill( Y_Z2,newweight*mc_weight_un[8]);

      if(debug) cout <<  "mc_weight_un=" << mc_weight_un[3] << endl;  
    // Format lepton syncronization                                                                                                                                      
      // {run}:{lumi}:{event}:{pdgId}:{pT:.2f}:{eta:.2f}:{phi:.2f}{SIP:.2f}:{PFChargedHadIso:.2f}:{PFNeutralHadIso:.2f}:{PFPhotonIso:.2f}:{PUCorr:.2f}:{combRelIsoPF:.3f}:{eleBDT:.3f}:{photpT:.2f}:{photDR:.2f}:{photRelIso:.2f}          
      
      Char_t leptformat[20000];
/*
      for(int i = 0; i < 4; ++i){      
	bool ismuon=false,iselectron=false;
	float dummy=0.;
	int flagFSR_tag=-999;
	int pfsr=-999;
	
	for (int j=0;j<RECO_NMU;j++) {	    
	  if (RECOMU_PT[j]==sortedpT[i]) {
	    ismuon=true;
	    break;
	  }
	  else {
	    for (int j=0;j<RECO_NELE;j++) {	    
	      if (RECOELE_PT[j]==sortedpT[i]) {
		iselectron=true;
		break;
	      }
	    }
	  }
	}
	
	for (int j=0;j<RECO_NELE;j++) {	    
	  if (RECOELE_PT[j]==sortedpT[i]) iselectron=true;
	}
	
	
	cout << "isMuon= " << ismuon << " and isElectron= " << iselectron << endl;

	if (ismuon){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == ipt[i] && iLp_tagEM[ p ] == 0 )  {
	      cout << "Muon with pT= " << RECOMU_PT[ipt[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	      flagFSR_tag=0;
	      pfsr=p;
	      break;
	    }
	  }
	}
	    
	else if (iselectron){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == ipt[i] && iLp_tagEM[ p ] == 1 )  {
	      cout << "Electron with pT= " << RECOELE_PT[ipt[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	      flagFSR_tag=1;
	      pfsr=p;
	      break;
	    }
	  }
	}
	

	if (ismuon && flagFSR_tag==0){
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f:%.2f:%.2f:%.2f",
		   Run,LumiSection,Event,
		   int(-13*RECOMU_CHARGE[ipt[i]]),
		   RECOMU_PT[ipt[i]],RECOMU_ETA[ipt[i]],RECOMU_PHI[ipt[i]],RECOMU_SIP[ipt[i]],
		   RECOMU_PFchHad[ipt[i]],RECOMU_PFneuHad[ipt[i]],RECOMU_PFphoton[ipt[i]],RECOMU_PFPUchAllPart[ipt[i]],RECOMU_PFX_dB[ipt[i]],dummy,
		   RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_DR[iLp[pfsr]],RECOPFPHOT_PFX_rho[iLp[pfsr]]
		   );
	}
	else if (iselectron && flagFSR_tag==1){
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f:%.2f:%.2f:%.2f",
		   Run,LumiSection,Event,
		   int(-11*RECOELE_CHARGE[ipt[i]]),
		   RECOELE_PT[ipt[i]],RECOELE_ETA[ipt[i]],RECOELE_PHI[ipt[i]],RECOELE_SIP[ipt[i]],
		   RECOELE_PFchHad[ipt[i]],RECOELE_PFneuHad[ipt[i]],RECOELE_PFphoton[ipt[i]],RHO_ele,RECOELE_PFX_rho[ipt[i]],RECOELE_mvaNonTrigV0[ipt[i]],
		   RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_DR[iLp[pfsr]],RECOPFPHOT_PFX_rho[iLp[pfsr]]
		   );
	}	
	else if (ismuon){	  	  	      
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f",
		   Run,LumiSection,Event,
		   int(-13*RECOMU_CHARGE[ipt[i]]),
		   RECOMU_PT[ipt[i]],RECOMU_ETA[ipt[i]],RECOMU_PHI[ipt[i]],RECOMU_SIP[ipt[i]],
		   RECOMU_PFchHad[ipt[i]],RECOMU_PFneuHad[ipt[i]],RECOMU_PFphoton[ipt[i]],RECOMU_PFPUchAllPart[ipt[i]],RECOMU_PFX_dB[ipt[i]],dummy
		   );
	}
	else if (iselectron){
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f",
		   Run,LumiSection,Event,
		   int(-11*RECOELE_CHARGE[ipt[i]]),
		   RECOELE_PT[ipt[i]],RECOELE_ETA[ipt[i]],RECOELE_PHI[ipt[i]],RECOELE_SIP[ipt[i]],
		   RECOELE_PFchHad[ipt[i]],RECOELE_PFneuHad[ipt[i]],RECOELE_PFphoton[ipt[i]],RHO_ele,RECOELE_PFX_rho[ipt[i]],RECOELE_mvaNonTrigV0[ipt[i]]
		   );
	}	  
	  
	output_txt  << leptformat << endl;

       
      }//end fill leptons
*/      
      
      
      //hIso_5->Fill( Iso_max,newweight ) ;
      //hSip_5->Fill( Sip_max,newweight ) ;
      //hIp_5->Fill( Ip_max,newweight ) ;

      // N.B. Do NOT Update the Isolation values and correct the 4 momenta of leptons for FSR
      for(int i = 0; i < N_good; ++i){
	int flagFSR=0;
	int pfsr=-999;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if (iLp[p]==-1) continue;
	  if (iLp_l[p]==-1) continue;
	  
	  if(debug) cout << "Index of lepton with photon ISR= " << iLp_l[ p ] << " and final lepton index= " << iL[i] << endl;
	  if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )  {
	    if(debug) cout << "Muon with pT= " << RECOMU_PT[iL[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	    
	    flagFSR=1;
	    pfsr=p;
	  }
	}
	
	
	if (flagFSR==1){
	  if(debug) cout << "Before correcting for FSR; muon pT= " << RECOMU_PT[iL[i]] << " Eta= " << RECOMU_ETA[iL[i]] << " Phi= " << RECOMU_PHI[iL[i]] << endl;
	  TLorentzVector Lept,LeptCorrection;
	  Lept.SetPtEtaPhiM(RECOMU_PT[iL[i]], RECOMU_ETA[iL[i]], RECOMU_PHI[iL[i]], 0.105);
	  LeptCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_ETA[iLp[pfsr]],RECOPFPHOT_PHI[iLp[pfsr]],0);
	  Lept+=LeptCorrection;
	  RECOMU_PT[iL[i]]=Lept.Pt();
	  RECOMU_ETA[iL[i]]=Lept.Eta();
	  RECOMU_PHI[iL[i]]=Lept.Phi();
	  if(debug) cout << "After correcting for FSR; muon pT= " << RECOMU_PT[iL[i]] << " Eta= " << RECOMU_ETA[iL[i]] << " Phi= " << RECOMU_PHI[iL[i]] << endl;
	}
      }
      
      // N.B. DO NOT Update the Isolation values and correct the 4 momenta of leptons for FSR
      for(int i = 0; i < Ne_good; ++i){
	int flagFSR=0;
	int pfsr=-999;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if (iLp[p]==-1) continue;
	  if (iLp_l[p]==-1) continue;
	  
	  if(debug) cout << "Index of lepton with photon ISR= " << iLp_l[ p ] << " and final lepton index= " << iLe[i] << endl;
	  if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {
	    if(debug) cout << "Electron with pT= " << RECOELE_PT[iLe[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	    // RECOELE_PFX_rho_new[iLe[i]]=
	    //   (RECOELE_PFchHad[iLe[i]]+
	    //    max(0.,RECOELE_PFneuHad[iLe[i]]+
	    // 	   (RECOELE_PFphoton[iLe[i]]-RECOPFPHOT_PT[iLp[p]] )-
	    // 	   max(RHO_ele,0.0)*(EffectiveArea)))/RECOELE_PT[iLe[i]];	    
	    flagFSR=1;
	    pfsr=p;
	  }
	}
	
	if (flagFSR==1){
	  if(debug) cout << "Before correcting for FSR; electron pT= " << RECOELE_PT[iLe[i]] << " Eta= " << RECOELE_ETA[iLe[i]] << " Phi= " << RECOELE_PHI[iLe[i]] << endl;
	  TLorentzVector Lept,LeptCorrection;
	  Lept.SetPtEtaPhiM(RECOELE_PT[iLe[i]], RECOELE_ETA[iLe[i]], RECOELE_PHI[iLe[i]], 0.105);
	  LeptCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_ETA[iLp[pfsr]],RECOPFPHOT_PHI[iLp[pfsr]],0);
	  Lept+=LeptCorrection;
	  RECOELE_PT[iLe[i]]=Lept.Pt();
	  RECOELE_ETA[iLe[i]]=Lept.Eta();
	  RECOELE_PHI[iLe[i]]=Lept.Phi();
	  if(debug) cout << "After correcting for FSR; muon pT= " << RECOELE_PT[iLe[i]] << " Eta= " << RECOELE_ETA[iLe[i]] << " Phi= " << RECOELE_PHI[iLe[i]] << endl;
	}
      }
      
      if(debug) cout << "Kinematics of leptons corrected for FSR photons (if existing)" << endl;
      
      
      // // **** Step 6:
      //  // QCD suppression: mll>4 GeV cut on all OS-SF pairs (4/4)           
     //if( min_mass_2L <= 4 ) continue ;
     
     ++N_6 ;  // fill counter
     N_6_w=N_6_w+newweight;

     // **** Step 7:
     // mass4l > 70 


   //  cout << "lep_type=" << lep_type[2] << endl; 
     // Leptons PT, ETA, Phi, Isol corrected for FSR
     for(int i = 0; i < 4; ++i){
       
       if( i == 0 && lep_type[i]==1){
	 hPtLep1_7->Fill( RECOMU_PT[ ipt[i] ],newweight ) ;
	 hYLep1_7->Fill(  RECOMU_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep1_7->Fill( RECOMU_PFX_dB_new[ ipt[i] ],newweight ) ;
	 hSipLep1_7->Fill( RECOMU_SIP[ ipt[i] ],newweight ) ;
	 hIpLep1_7->Fill( RECOMU_IP[ ipt[i] ],newweight ) ;
	 hIpErLep1_7->Fill( RECOMU_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 1 && lep_type[i]==1){
	 hPtLep2_7->Fill( RECOMU_PT[ ipt[i] ],newweight ) ;
	 hYLep2_7->Fill(  RECOMU_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep2_7->Fill( RECOMU_PFX_dB_new[ ipt[i] ],newweight ) ;
	 hSipLep2_7->Fill( RECOMU_SIP[ ipt[i] ],newweight ) ;
	 hIpLep2_7->Fill( RECOMU_IP[ ipt[i] ],newweight ) ;
	 hIpErLep2_7->Fill( RECOMU_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 2 && lep_type[i]==1){
         hPtLep3_7->Fill( RECOMU_PT[ ipt[i] ],newweight ) ;
         hYLep3_7->Fill(  RECOMU_ETA[ ipt[i] ],newweight ) ;
         hIsoLep3_7->Fill( RECOMU_PFX_dB_new[ ipt[i] ],newweight ) ;
         hSipLep3_7->Fill( RECOMU_SIP[ ipt[i] ],newweight ) ;
         hIpLep3_7->Fill( RECOMU_IP[ ipt[i] ],newweight ) ;
         hIpErLep3_7->Fill( RECOMU_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 3 && lep_type[i]==1){
         hPtLep4_7->Fill( RECOMU_PT[ ipt[i] ],newweight ) ;
         hYLep4_7->Fill(  RECOMU_ETA[ ipt[i] ],newweight ) ;
         hIsoLep4_7->Fill( RECOMU_PFX_dB_new[ ipt[i] ],newweight ) ;
         hSipLep4_7->Fill( RECOMU_SIP[ ipt[i] ],newweight ) ;
         hIpLep4_7->Fill( RECOMU_IP[ ipt[i] ],newweight ) ;
         hIpErLep4_7->Fill( RECOMU_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 0 && lep_type[i]==2){
         hPtLep1_7->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
         hYLep1_7->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
         hIsoLep1_7->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
         hSipLep1_7->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
         hIpLep1_7->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
         hIpErLep1_7->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 1 && lep_type[i]==2){
         hPtLep2_7->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
         hYLep2_7->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
         hIsoLep2_7->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
         hSipLep2_7->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
         hIpLep2_7->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
         hIpErLep2_7->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 2 && lep_type[i]==2){
         hPtLep3_7->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
         hYLep3_7->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
         hIsoLep3_7->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
         hSipLep3_7->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
         hIpLep3_7->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
         hIpErLep3_7->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 3 && lep_type[i]==2){
         hPtLep4_7->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
         hYLep4_7->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
         hIsoLep4_7->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
         hSipLep4_7->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
         hIpLep4_7->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
         hIpErLep4_7->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       

     }//end fill leptons
     
     //hIso_7->Fill( Iso_max,newweight ) ;
     //hSip_7->Fill( Sip_max,newweight ) ;
     //hIp_7->Fill( Ip_max,newweight ) ;
     
     
     
     
     
     ++N_8_PFMET;
     N_8_PFMET_w=N_8_PFMET_w+newweight;
     
     hPFMET_8->Fill(RECO_CORMETMUONS,newweight);
     hPFMET_8_pdf_up->Fill( RECO_CORMETMUONS,newweight*(1+PDF_weighting_un));
     hPFMET_8_pdf_dow->Fill( RECO_CORMETMUONS,newweight*(1-PDF_weighting_un));
     hPFMET_8_pu_up->Fill( RECO_CORMETMUONS,newweight*pu_up);
     hPFMET_8_pu_dow->Fill( RECO_CORMETMUONS,newweight*pu_dow);

     hPFMET_8_mc_1->Fill(RECO_CORMETMUONS,newweight*mc_weight_un[1]);
     hPFMET_8_mc_2->Fill(RECO_CORMETMUONS,newweight*mc_weight_un[2]);
     hPFMET_8_mc_3->Fill(RECO_CORMETMUONS,newweight*mc_weight_un[3]);
     hPFMET_8_mc_4->Fill(RECO_CORMETMUONS,newweight*mc_weight_un[4]);
     hPFMET_8_mc_5->Fill(RECO_CORMETMUONS,newweight*mc_weight_un[5]);
     hPFMET_8_mc_6->Fill(RECO_CORMETMUONS,newweight*mc_weight_un[6]);
     hPFMET_8_mc_7->Fill(RECO_CORMETMUONS,newweight*mc_weight_un[7]);
     hPFMET_8_mc_8->Fill(RECO_CORMETMUONS,newweight*mc_weight_un[8]);
 
     
     
     //Basic cuts to jets AND delta R section
     int njets_pass=0;
     int nbtag_pass=0;
     int nhb_pass=0,nlb_pass=0;
     TLorentzVector JET1,JET2,BOT1,BOT2;
     TLorentzVector JET1_up,JET2_up,BOT1_up,BOT2_up,JET1_jer_up,JET2_jer_up,BOT1_jer_up,BOT2_jer_up;
     TLorentzVector JET1_dow,JET2_dow,BOT1_dow,BOT2_dow,JET1_jer_dow,JET2_jer_dow,BOT1_jer_dow,BOT2_jer_dow;
     bool b1=false;
     bool b2=false;
     int jet1=-999,jet2=-999,bot1=-999,bot2=-999;      
     int jetfail[100];
     float GOOD_JET_PT_MAX = 0.;
     float JET_PHI_PT_MAX = 0.;
     float max_dphi_jet_met = 0.;
     float min_dphi_jet_met = 999.;
     vector <int>  v_good_jets_index;
     float nbjetweight=1;
     float scaleFactor1=1, scaleFactor2=1,scaleFactor1_up=1,scaleFactor1_dow=1,scaleFactor2_up=1,scaleFactor2_dow=1,scale_h_up=0,scale_h_dow=0,scale_l_up=0,scale_l_dow=0,scale_a_up=0,scale_a_dow=0;


     for(int i=0;i<100;i++) jetfail[i]=0;
     
     for(int i=0;i<RECO_PFJET_N;i++){

          if(RECO_PFJET_PT[i]<-100) continue;

       if(debug) cout<<i<<" Jet with pt= "<<RECO_PFJET_PT[i]<<" ETA "<<RECO_PFJET_ETA[i]<<" PUID "<<RECO_PFJET_PUID[i] << " PUID_MVA "<< RECO_PFJET_PUID_MVA[i]<<endl;

       double scaleFactor_1=1;
       double scaleFactor_1_up=1;
       double scaleFactor_1_dow=1;

       double Pt=RECO_PFJET_PT[i];
       double Eta=RECO_PFJET_ETA[i];
      if(RECO_PFJET_PT[i]>25. && fabs(RECO_PFJET_ETA[i])<2.4){
       hPtJet_7->Fill(RECO_PFJET_PT[i],newweight);
       hYJet_7->Fill(RECO_PFJET_ETA[i],newweight);
              if( MC_type == "Spring16" && DATA_type == "NO"){
                 if(RECOBOT_MatchingMCTruth[i+1]==2) scaleFactor_1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, Eta, Pt);
                 else if(RECOBOT_MatchingMCTruth[i+1]==1) scaleFactor_1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, Eta, Pt);
                 else if(RECOBOT_MatchingMCTruth[i+1]==5) scaleFactor_1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                 else scaleFactor_1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
              }
          if(cSV_BTagJet_DISCR[i]> 0.5426) {
          hPtBot_7->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1);}
       }
  
       float dphi_jet_met=0.;
       dphi_jet_met=RECO_PFJET_PHI[i]-RECO_PFMET_PHI;

       // JET smearing

        bool match_jet = false ;
        double jer_scale_up=1;
        double jer_scale_dow=1;
        double jer_scale=1;
        int aa=0;
        for(int p=0; p<9; p++){
          if(abs(RECO_PFJET_ETA[i])>=eta_jer[p]&&abs(RECO_PFJET_ETA[i])<eta_jer[p+1]) aa=p;
        }

        for(int k=0;k<100;k++){

             if(MC_GENJET_PT[k]<-1) continue;

             double deltaR = sqrt( pow( DELTAPHI( RECO_PFJET_PHI[i] , MC_GENJET_PHI[k] ),2) + pow(RECO_PFJET_ETA[i] - MC_GENJET_ETA[k],2) );
         if(deltaR<0.2) {
             match_jet=true;
             jer_scale_up=1+un_jer[aa]*abs(RECO_PFJET_PT[i]-MC_GENJET_PT[k])/RECO_PFJET_PT[i];
             jer_scale_dow=1-un_jer[aa]*abs(RECO_PFJET_PT[i]-MC_GENJET_PT[k])/RECO_PFJET_PT[i];
           }
       }

       if(!match_jet){

           double scale = std::sqrt(std::max(scale_jer[aa]*scale_jer[aa] - 1.0 , 0.0));

           TRandom3 rand1(0);
           jer_scale_up=1+un_jer[aa]*scale;
           jer_scale_dow=1-un_jer[aa]*scale;

       }

       bool goodjet = RECO_PFJET_NHF[i] < 0.99 &&
                      RECO_PFJET_NEF[i] < 0.99 &&
                      RECO_PFJET_CHF[i] < 0.99 &&
                      RECO_PFJET_CEF[i] < 0.99 &&
                      RECO_PFJET_nconstituents[i] > 1 &&
                      RECO_PFJET_NCH > 0;
       
 
       if(RECO_PFJET_PT[i]>25. && fabs(RECO_PFJET_ETA[i])<2.4 /*&& goodjet==1*/){
       
      	 //Check that jet has deltaR>0.4 away from any tight lepton corrected for FSR
	 for(int mu = 0; mu < N_good; ++mu){
	//   if (fabs(RECOMU_SIP[iL[mu]])>=4.) continue;
      	   if (RECOMU_PFX_dB_new[iL[mu]]>=0.20) continue;
	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOMU_PHI[iL[mu]]),2) + pow(RECO_PFJET_ETA[i] - RECOMU_ETA[iL[mu]],2));
	   if(debug) cout << "1st lepton muon: " << " pT=" << RECOMU_PT[iL[mu]] <<" deltaR "<< deltaR <<endl;	   
	   if (deltaR<0.4){
	     jetfail[i]=1;
     	     if(debug) cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
      	 for(int ele = 0; ele < Ne_good; ++ele){
//    	   if (fabs(RECOELE_SIP[iLe[ele]])>=4.) continue;
	   if (RECOELE_PFX_rho_new[iLe[ele]]>=0.20) continue;
      	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOELE_PHI[iLe[ele]]),2) + pow(RECO_PFJET_ETA[i] - RECOELE_ETA[iLe[ele]],2));
     	   if(debug) cout << "1st lepton electron: " << " pT=" << RECOELE_PT[iLe[ele]] <<" deltaR "<< deltaR <<endl;
	   if (deltaR<0.4){
     	     jetfail[i]=1;
     	     if(debug) cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
	 // cleaning w.r.t FSR photons attached to leptons
	 for(int j=0.;j<Nphotons;j++) {
           if (iLp_l[j]!=-1 && (iLp_tagEM[j]==0 || iLp_tagEM[j]==1) ) {
	     if (iLp_tagEM[j]==0&&debug) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[j]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[j]] << endl;
	     if (iLp_tagEM[j]==1&&debug) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[j]] << " attached to a electron with pT= " << RECOELE_PT[iLp_l[j]] << endl;
	     double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOPFPHOT_PHI[iLp[j]]),2) + pow(RECO_PFJET_ETA[i] - RECOPFPHOT_ETA[iLp[j]],2));
	     if (deltaR<0.4){
	       jetfail[i]=1;
	       if(debug) cout << " jetfail " << jetfail[i] <<endl;
	       break;
	     }
	   }
         }
	 // 


	 if (jetfail[i]==0){
	   if(debug) cout<< " PASS jet " <<i<<" PT= "<<RECO_PFJET_PT[i]<<" ETA= "<<RECO_PFJET_ETA[i]<<" PUID= "<<RECO_PFJET_PUID[i]<<endl;
	   njets_pass++;

           hPtJet_8->Fill(RECO_PFJET_PT[i],newweight);
           hYJet_8->Fill(RECO_PFJET_ETA[i],newweight);
           hPtJet_8_up->Fill(RECO_PFJET_PT_UP[i],newweight);
           hPtJet_8_dow->Fill(RECO_PFJET_PT_DOW[i],newweight);

           hPtJet_8_pdf_up->Fill(RECO_PFJET_PT[i],newweight*(1+PDF_weighting_un));
           hPtJet_8_pdf_dow->Fill(RECO_PFJET_PT[i],newweight*(1-PDF_weighting_un));
           hPtJet_8_pu_up->Fill(RECO_PFJET_PT[i],newweight*pu_up);
           hPtJet_8_pu_dow->Fill(RECO_PFJET_PT[i],newweight*pu_dow);

           hYJet_8_pdf_up->Fill(RECO_PFJET_ETA[i],newweight*(1+PDF_weighting_un));
           hYJet_8_pdf_dow->Fill(RECO_PFJET_ETA[i],newweight*(1-PDF_weighting_un));

           hPtJet_8_jer_up->Fill(RECO_PFJET_PT[i]*jer_scale_up,newweight);
           hPtJet_8_jer_dow->Fill(RECO_PFJET_PT[i]*jer_scale_dow,newweight);

           hPtJet_8_mc_1->Fill(RECO_PFJET_PT[i],newweight*mc_weight_un[1]);
           hYJet_8_mc_1->Fill(RECO_PFJET_ETA[i],newweight*mc_weight_un[1]);

           hPtJet_8_mc_2->Fill(RECO_PFJET_PT[i],newweight*mc_weight_un[2]);
           hYJet_8_mc_2->Fill(RECO_PFJET_ETA[i],newweight*mc_weight_un[2]);

           hPtJet_8_mc_3->Fill(RECO_PFJET_PT[i],newweight*mc_weight_un[3]);
           hYJet_8_mc_3->Fill(RECO_PFJET_ETA[i],newweight*mc_weight_un[3]);

           hPtJet_8_mc_4->Fill(RECO_PFJET_PT[i],newweight*mc_weight_un[4]);
           hYJet_8_mc_4->Fill(RECO_PFJET_ETA[i],newweight*mc_weight_un[4]);

           hPtJet_8_mc_5->Fill(RECO_PFJET_PT[i],newweight*mc_weight_un[5]);
           hYJet_8_mc_5->Fill(RECO_PFJET_ETA[i],newweight*mc_weight_un[5]);

           hPtJet_8_mc_6->Fill(RECO_PFJET_PT[i],newweight*mc_weight_un[6]);
           hYJet_8_mc_6->Fill(RECO_PFJET_ETA[i],newweight*mc_weight_un[6]);

           hPtJet_8_mc_7->Fill(RECO_PFJET_PT[i],newweight*mc_weight_un[7]);
           hYJet_8_mc_7->Fill(RECO_PFJET_ETA[i],newweight*mc_weight_un[7]);

           hPtJet_8_mc_8->Fill(RECO_PFJET_PT[i],newweight*mc_weight_un[8]);
           hYJet_8_mc_8->Fill(RECO_PFJET_ETA[i],newweight*mc_weight_un[8]);

           //b-tagging

         if( MC_type == "Spring16" && DATA_type == "NO"){
            if(RECOBOT_MatchingMCTruth[i+1]==1){
             hPtJet_8_b->Fill(RECO_PFJET_PT[i],newweight);
            }
            else if(RECOBOT_MatchingMCTruth[i+1]==2){
             hPtJet_8_c->Fill(RECO_PFJET_PT[i],newweight);
            }
            else if(RECOBOT_MatchingMCTruth[i+1]==3){
             hPtJet_8_l->Fill(RECO_PFJET_PT[i],newweight);
            }
            else if(RECOBOT_MatchingMCTruth[i+1]==5){
             hPtJet_8_o->Fill(RECO_PFJET_PT[i],newweight);
            }
           }

           double scaleFactor=1;
           double scaleFactor_up=1;
           double scaleFactor_dow=1;


           if(cSV_BTagJet_DISCR[i]> 0.5426){
              nbtag_pass++;
              if( MC_type == "Spring16" && DATA_type == "NO"){
                 if(RECOBOT_MatchingMCTruth[i+1]==2){
                    scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, Eta, Pt);
                    scaleFactor_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_C, Eta, Pt);
                    scaleFactor_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_C, Eta, Pt);
                    hPtBot_8_c->Fill(RECO_PFJET_PT[i],newweight*scaleFactor);
                    hPtBot_8_c_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_up);
                    hPtBot_8_c_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_dow);
                    scale_h_up+=(scaleFactor_up/scaleFactor-1)*(scaleFactor_up/scaleFactor-1);
                    scale_h_dow+=(1-scaleFactor_dow/scaleFactor)*(1-scaleFactor_dow/scaleFactor);
                    nhb_pass++;
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==1){
                    scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, Eta, Pt);
                    scaleFactor_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_B, Eta, Pt);
                    scaleFactor_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_B, Eta, Pt);
                    hPtBot_8_b->Fill(RECO_PFJET_PT[i],newweight*scaleFactor);
                    hPtBot_8_b_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_up);
                    hPtBot_8_b_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_dow);
                    scale_h_up+=(scaleFactor_up/scaleFactor-1)*(scaleFactor_up/scaleFactor-1);
                    scale_h_dow+=(1-scaleFactor_dow/scaleFactor)*(1-scaleFactor_dow/scaleFactor); 
                   nhb_pass++;
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==3){
                    double sf_l=1;
//                    int binxb = b_eff_zz_p5->GetXaxis()->FindBin(Pt);
//                    if( datasetName.Contains("ZZTo4L")|| datasetName.Contains("GluGluTo")) sf_l=b_eff_zz_p5->GetBinContent(binxb);
                    if(sf_l<=0) sf_l=1;
                      
                    scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor*=sf_l;
                    scaleFactor_up = sf_l*reader.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_dow = sf_l*reader.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    hPtBot_8_l->Fill(RECO_PFJET_PT[i],newweight*scaleFactor);
                    hPtBot_8_l_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_up);
                    hPtBot_8_l_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_dow);
                    scale_l_up+=(scaleFactor_up/scaleFactor-1)*(scaleFactor_up/scaleFactor-1);
                    scale_l_dow+=(1-scaleFactor_dow/scaleFactor)*(1-scaleFactor_dow/scaleFactor);
                    nlb_pass++;
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==5){
                    scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    hPtBot_8_o->Fill(RECO_PFJET_PT[i],newweight*scaleFactor);
                    hPtBot_8_o_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_up);
                    hPtBot_8_o_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_dow);
                    scale_l_up+=(scaleFactor_up/scaleFactor-1)*(scaleFactor_up/scaleFactor-1);
                    scale_l_dow+=(1-scaleFactor_dow/scaleFactor)*(1-scaleFactor_dow/scaleFactor);
                    nlb_pass++;
                 }
                    scale_a_up+=(scaleFactor_up/scaleFactor-1)*(scaleFactor_up/scaleFactor-1);
                    scale_a_dow+=(1-scaleFactor_dow/scaleFactor)*(1-scaleFactor_dow/scaleFactor);

           }
           hPtBot_8->Fill(RECO_PFJET_PT[i],newweight*scaleFactor);
           hPtBot_8_up->Fill(RECO_PFJET_PT_UP[i],newweight*scaleFactor);
           hPtBot_8_dow->Fill(RECO_PFJET_PT_DOW[i],newweight*scaleFactor);
           hPtBot_8_pdf_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*(1+PDF_weighting_un));
           hPtBot_8_pdf_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*(1-PDF_weighting_un));
           hPtBot_8_pu_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*pu_up);
           hPtBot_8_pu_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*pu_dow);
           hPtBot_8_jer_up->Fill(RECO_PFJET_PT[i]*jer_scale_up,newweight);
           hPtBot_8_jer_dow->Fill(RECO_PFJET_PT[i]*jer_scale_dow,newweight);

           hPtBot_8_mc_1->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*mc_weight_un[1]);
           hPtBot_8_mc_2->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*mc_weight_un[2]);
           hPtBot_8_mc_3->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*mc_weight_un[3]);
           hPtBot_8_mc_4->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*mc_weight_un[4]);
           hPtBot_8_mc_5->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*mc_weight_un[5]);
           hPtBot_8_mc_6->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*mc_weight_un[6]);
           hPtBot_8_mc_7->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*mc_weight_un[7]);
           hPtBot_8_mc_8->Fill(RECO_PFJET_PT[i],newweight*scaleFactor*mc_weight_un[8]);

           hYBot_8->Fill(RECO_PFJET_ETA[i],newweight*scaleFactor);


           }
          else{
             if( MC_type == "Spring16" && DATA_type == "NO"){
                 double b_sf=1;
                 double b_eff=0;
                 int binxb = b_eff_p1->GetXaxis()->FindBin(Pt);
                 if(RECOBOT_MatchingMCTruth[i+1]==2){
                    b_sf = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, Eta, Pt);
                    b_eff=b_eff_p2->GetBinContent(binxb);
                    if(b_eff>=1||b_eff<=0) b_eff=0;
                    scaleFactor=(1-b_eff*b_sf)/(1-b_eff);
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==1){
                    b_sf = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, Eta, Pt);
                    b_eff=b_eff_p1->GetBinContent(binxb);
                    if(b_eff>=1||b_eff<=0) b_eff=0;
                    scaleFactor=(1-b_eff*b_sf)/(1-b_eff);
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==3){
                    b_sf = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    b_eff=b_eff_p3->GetBinContent(binxb);
                    if(b_eff>=1||b_eff<=0) b_eff=0;
                    scaleFactor=(1-b_eff*b_sf)/(1-b_eff);
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==5){
                    b_sf = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    b_eff=b_eff_p4->GetBinContent(binxb);
                    if(b_eff>=1||b_eff<=0) b_eff=0;
                    scaleFactor=(1-b_eff*b_sf)/(1-b_eff);
                 }
              if(scaleFactor>50) scaleFactor=1;
              if(scaleFactor<0) scaleFactor=0;
             }
          }
         nbjetweight*=scaleFactor;

	   if (njets_pass==1){
	     jet1=i;
	     JET1.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

             JET1_jer_up.SetPtEtaPhiE(RECO_PFJET_PT[i]*jer_scale_up,RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             JET1_jer_dow.SetPtEtaPhiE(RECO_PFJET_PT[i]*jer_scale_dow,RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

             JET1_up.SetPtEtaPhiE(RECO_PFJET_PT_UP[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             JET1_dow.SetPtEtaPhiE(RECO_PFJET_PT_DOW[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

	     GOOD_JET_PT_MAX=RECO_PFJET_PT[i];
	     JET_PHI_PT_MAX=RECO_PFJET_PHI[i];
             bdiscr_5_lead->Fill(cSV_BTagJet_DISCR[i],newweight);
	     if(debug) cout<<"Among the jets that pass the jet with the highet pt is the jet of index "<< i <<". It has pt "<<GOOD_JET_PT_MAX<<". The corresponding value of phi is " << JET_PHI_PT_MAX <<endl;
  
	   }
	   if (njets_pass==2){
	     jet2=i;
	     JET2.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

             JET2_jer_up.SetPtEtaPhiE(RECO_PFJET_PT[i]*jer_scale_up,RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             JET2_jer_dow.SetPtEtaPhiE(RECO_PFJET_PT[i]*jer_scale_dow,RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

             JET2_up.SetPtEtaPhiE(RECO_PFJET_PT_UP[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             JET2_dow.SetPtEtaPhiE(RECO_PFJET_PT_DOW[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

             bdiscr_5_sub->Fill(cSV_BTagJet_DISCR[i],newweight);

	   }
           //b-tagging

           if (nbtag_pass==1&&cSV_BTagJet_DISCR[i]> 0.5426){
             bot1=i;
             if(RECOBOT_MatchingMCTruth[i+1]==1||RECOBOT_MatchingMCTruth[i+1]==2) b1=true;

                 if(RECOBOT_MatchingMCTruth[i+1]==2){
                    scaleFactor1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, Eta, Pt);
                    scaleFactor1_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_C, Eta, Pt);
                    scaleFactor1_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_C, Eta, Pt);
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==1){
                    scaleFactor1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, Eta, Pt);
                    scaleFactor1_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_B, Eta, Pt);
                    scaleFactor1_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_B, Eta, Pt);
                 }
                 else{
                    scaleFactor1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor1_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor1_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                 }

             BOT1.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             BOT1_up.SetPtEtaPhiE(RECO_PFJET_PT_UP[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             BOT1_dow.SetPtEtaPhiE(RECO_PFJET_PT_DOW[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));          
             BOT1_jer_up.SetPtEtaPhiE(RECO_PFJET_PT[i]*jer_scale_up,RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             BOT1_jer_dow.SetPtEtaPhiE(RECO_PFJET_PT[i]*jer_scale_dow,RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
           }
           if (nbtag_pass==2&&cSV_BTagJet_DISCR[i]> 0.5426){
             bot2=i;
             if(RECOBOT_MatchingMCTruth[i+1]==1||RECOBOT_MatchingMCTruth[i+1]==2) b2=true;
                 if(RECOBOT_MatchingMCTruth[i+1]==2){
                    scaleFactor2 = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, Eta, Pt);
                    scaleFactor2_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_C, Eta, Pt);
                    scaleFactor2_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_C, Eta, Pt);
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==1){
                    scaleFactor2 = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, Eta, Pt);
                    scaleFactor2_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_B, Eta, Pt);
                    scaleFactor2_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_B, Eta, Pt);
                 }
                 else{
                    scaleFactor2 = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor2_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor2_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                 }

             BOT2.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             BOT2_up.SetPtEtaPhiE(RECO_PFJET_PT_UP[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             BOT2_dow.SetPtEtaPhiE(RECO_PFJET_PT_DOW[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

             BOT2_jer_up.SetPtEtaPhiE(RECO_PFJET_PT[i]*jer_scale_up,RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             BOT2_jer_dow.SetPtEtaPhiE(RECO_PFJET_PT[i]*jer_scale_dow,RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

           }
	 }

       }
       else{
      	 jetfail[i]=1;
       }
       //cout<<" JETFAIL "<<jetfail[i]<<endl;
     }


     hNjets_8->Fill(njets_pass,newweight);
     hNbjets_8_btag_up->Fill(nbtag_pass,newweight*nbjetweight*(1+sqrt(scale_a_up)));
     hNbjets_8_btag_dow->Fill(nbtag_pass,newweight*nbjetweight*(1-sqrt(scale_a_dow)));
     hNjets_8_pdf_up->Fill(njets_pass,newweight*(1+PDF_weighting_un));
     hNjets_8_pdf_dow->Fill(njets_pass,newweight*(1-PDF_weighting_un));
     hNjets_8_pu_up->Fill(njets_pass,newweight*pu_up);
     hNjets_8_pu_dow->Fill(njets_pass,newweight*pu_dow);
     hNjets_8_mc_1->Fill(njets_pass,newweight*mc_weight_un[1]);
     hNjets_8_mc_2->Fill(njets_pass,newweight*mc_weight_un[2]);
     hNjets_8_mc_3->Fill(njets_pass,newweight*mc_weight_un[3]);
     hNjets_8_mc_4->Fill(njets_pass,newweight*mc_weight_un[4]);
     hNjets_8_mc_5->Fill(njets_pass,newweight*mc_weight_un[5]);
     hNjets_8_mc_6->Fill(njets_pass,newweight*mc_weight_un[6]);
     hNjets_8_mc_7->Fill(njets_pass,newweight*mc_weight_un[7]);
     hNjets_8_mc_8->Fill(njets_pass,newweight*mc_weight_un[8]);

     hNbjets_8->Fill(nbtag_pass,newweight*nbjetweight);
     hNbjets_8_pdf_up->Fill(nbtag_pass,newweight*nbjetweight*(1+PDF_weighting_un));
     hNbjets_8_pdf_dow->Fill(nbtag_pass,newweight*nbjetweight*(1-PDF_weighting_un));
     hNbjets_8_pu_up->Fill(nbtag_pass,newweight*nbjetweight*pu_up);
     hNbjets_8_pu_dow->Fill(nbtag_pass,newweight*nbjetweight*pu_dow);
     hNbjets_8_mc_1->Fill(nbtag_pass,newweight*nbjetweight*mc_weight_un[1]);
     hNbjets_8_mc_2->Fill(nbtag_pass,newweight*nbjetweight*mc_weight_un[2]);
     hNbjets_8_mc_3->Fill(nbtag_pass,newweight*nbjetweight*mc_weight_un[3]);
     hNbjets_8_mc_4->Fill(nbtag_pass,newweight*nbjetweight*mc_weight_un[4]);
     hNbjets_8_mc_5->Fill(nbtag_pass,newweight*nbjetweight*mc_weight_un[5]);
     hNbjets_8_mc_6->Fill(nbtag_pass,newweight*nbjetweight*mc_weight_un[6]);
     hNbjets_8_mc_7->Fill(nbtag_pass,newweight*nbjetweight*mc_weight_un[7]);
     hNbjets_8_mc_8->Fill(nbtag_pass,newweight*nbjetweight*mc_weight_un[8]);

     if(njets_pass<2) continue;
     double Mjj = (JET1+JET2).M();

     Mjj_6->Fill(Mjj,newweight);
     Mjj_6_up->Fill((JET1_up+JET2_up).M(),newweight);
     Mjj_6_dow->Fill((JET1_dow+JET2_dow).M(),newweight);

     Mjj_6_pdf_up->Fill(Mjj,newweight*(1+PDF_weighting_un));
     Mjj_6_pdf_dow->Fill(Mjj,newweight*(1-PDF_weighting_un));
     Mjj_6_pu_up->Fill(Mjj,newweight*pu_up);
     Mjj_6_pu_dow->Fill(Mjj,newweight*pu_dow);

     Mjj_6_jer_up->Fill((JET1_jer_up+JET2_jer_up).M(),newweight);
     Mjj_6_jer_dow->Fill((JET1_jer_dow+JET2_jer_dow).M(),newweight);

     Mjj_6_mc_1->Fill(Mjj,newweight*mc_weight_un[1]);
     Mjj_6_mc_2->Fill(Mjj,newweight*mc_weight_un[2]);
     Mjj_6_mc_3->Fill(Mjj,newweight*mc_weight_un[3]);
     Mjj_6_mc_4->Fill(Mjj,newweight*mc_weight_un[4]);
     Mjj_6_mc_5->Fill(Mjj,newweight*mc_weight_un[5]);
     Mjj_6_mc_6->Fill(Mjj,newweight*mc_weight_un[6]);
     Mjj_6_mc_7->Fill(Mjj,newweight*mc_weight_un[7]);
     Mjj_6_mc_8->Fill(Mjj,newweight*mc_weight_un[8]);

     if(nbtag_pass==0){
     hMZ1_61->Fill(massZ1,newweight);
     hPtZ1_61->Fill( ptZ1,newweight );
     hYZ1_61->Fill( Y_Z1,newweight );

     hMZ2_61->Fill(massZ2,newweight);
     hPtZ2_61->Fill( ptZ2,newweight );
     hYZ2_61->Fill( Y_Z2,newweight );

     Mjj_7->Fill(Mjj,newweight);
     Mjj_7_up->Fill((JET1_up+JET2_up).M(),newweight);
     Mjj_7_dow->Fill((JET1_dow+JET2_dow).M(),newweight);

     Mjj_7_pdf_up->Fill(Mjj,newweight*(1+PDF_weighting_un));
     Mjj_7_pdf_dow->Fill(Mjj,newweight*(1-PDF_weighting_un));
     Mjj_7_pu_up->Fill(Mjj,newweight*pu_up);
     Mjj_7_pu_dow->Fill(Mjj,newweight*pu_dow);

     Mjj_7_jer_up->Fill((JET1_jer_up+JET2_jer_up).M(),newweight);
     Mjj_7_jer_dow->Fill((JET1_jer_dow+JET2_jer_dow).M(),newweight);

     Mjj_7_mc_1->Fill(Mjj,newweight*mc_weight_un[1]);
     Mjj_7_mc_2->Fill(Mjj,newweight*mc_weight_un[2]);
     Mjj_7_mc_3->Fill(Mjj,newweight*mc_weight_un[3]);
     Mjj_7_mc_4->Fill(Mjj,newweight*mc_weight_un[4]);
     Mjj_7_mc_5->Fill(Mjj,newweight*mc_weight_un[5]);
     Mjj_7_mc_6->Fill(Mjj,newweight*mc_weight_un[6]);
     Mjj_7_mc_7->Fill(Mjj,newweight*mc_weight_un[7]);
     Mjj_7_mc_8->Fill(Mjj,newweight*mc_weight_un[8]);
     }

     if(nbtag_pass==1){
      double Mbj,Mbj_jer_up,Mbj_jer_dow,Mbj_up,Mbj_dow;
      if(jet1==bot1){
         Mbj=(BOT1+JET2).M();
         Mbj_jer_up=(BOT1_jer_up+JET2_jer_up).M();
         Mbj_jer_dow=(BOT1_jer_dow+JET2_jer_dow).M();
         Mbj_up=(BOT1_up+JET2_up).M();
         Mbj_dow=(BOT1_dow+JET2_dow).M();
      }
      else{ 
         Mbj=(JET1+BOT1).M();
         Mbj_jer_up=(BOT1_jer_up+JET1_jer_up).M();
         Mbj_jer_dow=(BOT1_jer_dow+JET1_jer_dow).M();
         Mbj_up=(BOT1_up+JET1_up).M();
         Mbj_dow=(BOT1_dow+JET1_dow).M();
      }
      Mbj_6->Fill(Mbj,newweight);

      hMZ1_62->Fill(massZ1,newweight);
      hPtZ1_62->Fill( ptZ1,newweight );
      hYZ1_62->Fill( Y_Z1,newweight );

      hMZ2_62->Fill(massZ2,newweight);
      hPtZ2_62->Fill( ptZ2,newweight );
      hYZ2_62->Fill( Y_Z2,newweight );

      Mbj_6_jer_up->Fill(Mbj_jer_up,newweight*nbjetweight);
      Mbj_6_jer_dow->Fill(Mbj_jer_dow,newweight*nbjetweight);

      Mbj_6_up->Fill(Mbj_up,newweight*nbjetweight);
      Mbj_6_dow->Fill(Mbj_dow,newweight*nbjetweight);

      Mbj_6_pdf_up->Fill(Mbj,newweight*nbjetweight*(1+PDF_weighting_un));
      Mbj_6_pdf_dow->Fill(Mbj,newweight*nbjetweight*(1-PDF_weighting_un));

      Mbj_6_pu_up->Fill(Mbj,newweight*nbjetweight*pu_up);
      Mbj_6_pu_dow->Fill(Mbj,newweight*nbjetweight*pu_dow);

      Mbj_6_mc_1->Fill(Mbj,newweight*nbjetweight*mc_weight_un[1]);
      Mbj_6_mc_2->Fill(Mbj,newweight*nbjetweight*mc_weight_un[2]);
      Mbj_6_mc_3->Fill(Mbj,newweight*nbjetweight*mc_weight_un[3]);
      Mbj_6_mc_4->Fill(Mbj,newweight*nbjetweight*mc_weight_un[4]);
      Mbj_6_mc_5->Fill(Mbj,newweight*nbjetweight*mc_weight_un[5]);
      Mbj_6_mc_6->Fill(Mbj,newweight*nbjetweight*mc_weight_un[6]);
      Mbj_6_mc_7->Fill(Mbj,newweight*nbjetweight*mc_weight_un[7]);
      Mbj_6_mc_8->Fill(Mbj,newweight*nbjetweight*mc_weight_un[8]);
      }

     hMZ1_6->Fill(massZ1,newweight);
     hPtZ1_6->Fill( ptZ1,newweight );
     hYZ1_6->Fill( Y_Z1,newweight );

      hMZ1_6_pdf_up->Fill( massZ1,newweight*(1+PDF_weighting_un));
      hPtZ1_6_pdf_up->Fill( ptZ1,newweight*(1+PDF_weighting_un));
      hYZ1_6_pdf_up->Fill( Y_Z1,newweight*(1+PDF_weighting_un));

      hMZ1_6_pdf_dow->Fill( massZ1,newweight*(1-PDF_weighting_un));
      hPtZ1_6_pdf_dow->Fill( ptZ1,newweight*(1-PDF_weighting_un));
      hYZ1_6_pdf_dow->Fill( Y_Z1,newweight*(1-PDF_weighting_un));

      hMZ1_6_pu_up->Fill( massZ1,newweight*pu_up);
      hPtZ1_6_pu_up->Fill( ptZ1,newweight*pu_up);
      hYZ1_6_pu_up->Fill( Y_Z1,newweight*pu_up);

      hMZ1_6_pu_dow->Fill( massZ1,newweight*pu_dow);
      hPtZ1_6_pu_dow->Fill( ptZ1,newweight*pu_dow);
      hYZ1_6_pu_dow->Fill( Y_Z1,newweight*pu_dow);

      hMZ1_6_mc_1->Fill( massZ1,newweight*mc_weight_un[1] );
      hPtZ1_6_mc_1->Fill( ptZ1,newweight*mc_weight_un[1]);
      hYZ1_6_mc_1->Fill( Y_Z1,newweight*mc_weight_un[1]);

      hMZ1_6_mc_2->Fill( massZ1,newweight*mc_weight_un[2] );
      hPtZ1_6_mc_2->Fill( ptZ1,newweight*mc_weight_un[2]);
      hYZ1_6_mc_2->Fill( Y_Z1,newweight*mc_weight_un[2]);

      hMZ1_6_mc_3->Fill( massZ1,newweight*mc_weight_un[3] );
      hPtZ1_6_mc_3->Fill( ptZ1,newweight*mc_weight_un[3]);
      hYZ1_6_mc_3->Fill( Y_Z1,newweight*mc_weight_un[3]);

      hMZ1_6_mc_4->Fill( massZ1,newweight*mc_weight_un[4] );
      hPtZ1_6_mc_4->Fill( ptZ1,newweight*mc_weight_un[4]);
      hYZ1_6_mc_4->Fill( Y_Z1,newweight*mc_weight_un[4]);

      hMZ1_6_mc_5->Fill( massZ1,newweight*mc_weight_un[5] );
      hPtZ1_6_mc_5->Fill( ptZ1,newweight*mc_weight_un[5]);
      hYZ1_6_mc_5->Fill( Y_Z1,newweight*mc_weight_un[5]);

      hMZ1_6_mc_6->Fill( massZ1,newweight*mc_weight_un[6] );
      hPtZ1_6_mc_6->Fill( ptZ1,newweight*mc_weight_un[6]);
      hYZ1_6_mc_6->Fill( Y_Z1,newweight*mc_weight_un[6]);

      hMZ1_6_mc_7->Fill( massZ1,newweight*mc_weight_un[7] );
      hPtZ1_6_mc_7->Fill( ptZ1,newweight*mc_weight_un[7]);
      hYZ1_6_mc_7->Fill( Y_Z1,newweight*mc_weight_un[7]);

      hMZ1_6_mc_8->Fill( massZ1,newweight*mc_weight_un[8] );
      hPtZ1_6_mc_8->Fill( ptZ1,newweight*mc_weight_un[8]);
      hYZ1_6_mc_8->Fill( Y_Z1,newweight*mc_weight_un[8]);

      hMZ2_6->Fill(massZ2,newweight);
      hPtZ2_6->Fill( ptZ2,newweight );
      hYZ2_6->Fill( Y_Z2,newweight );

      hMZ2_6_pdf_up->Fill( massZ2,newweight*(1+PDF_weighting_un));
      hPtZ2_6_pdf_up->Fill( ptZ2,newweight*(1+PDF_weighting_un));
      hYZ2_6_pdf_up->Fill( Y_Z2,newweight*(1+PDF_weighting_un));

      hMZ2_6_pdf_dow->Fill( massZ2,newweight*(1-PDF_weighting_un));
      hPtZ2_6_pdf_dow->Fill( ptZ2,newweight*(1-PDF_weighting_un));
      hYZ2_6_pdf_dow->Fill( Y_Z2,newweight*(1-PDF_weighting_un));

      hMZ2_6_pu_up->Fill( massZ2,newweight*pu_up);
      hPtZ2_6_pu_up->Fill( ptZ2,newweight*pu_up);
      hYZ2_6_pu_up->Fill( Y_Z2,newweight*pu_up);

      hMZ2_6_pu_dow->Fill( massZ2,newweight*pu_dow);
      hPtZ2_6_pu_dow->Fill( ptZ2,newweight*pu_dow);
      hYZ2_6_pu_dow->Fill( Y_Z2,newweight*pu_dow);
   
      hMZ2_6_mc_1->Fill( massZ2,newweight*mc_weight_un[1] );
      hPtZ2_6_mc_1->Fill( ptZ2,newweight*mc_weight_un[1]);
      hYZ2_6_mc_1->Fill( Y_Z2,newweight*mc_weight_un[1]);

      hMZ2_6_mc_2->Fill( massZ2,newweight*mc_weight_un[2] );
      hPtZ2_6_mc_2->Fill( ptZ2,newweight*mc_weight_un[2]);
      hYZ2_6_mc_2->Fill( Y_Z2,newweight*mc_weight_un[2]);

      hMZ2_6_mc_3->Fill( massZ2,newweight*mc_weight_un[3] );
      hPtZ2_6_mc_3->Fill( ptZ2,newweight*mc_weight_un[3]);
      hYZ2_6_mc_3->Fill( Y_Z2,newweight*mc_weight_un[3]);

      hMZ2_6_mc_4->Fill( massZ2,newweight*mc_weight_un[4] );
      hPtZ2_6_mc_4->Fill( ptZ2,newweight*mc_weight_un[4]);
      hYZ2_6_mc_4->Fill( Y_Z2,newweight*mc_weight_un[4]);

      hMZ2_6_mc_5->Fill( massZ2,newweight*mc_weight_un[5] );
      hPtZ2_6_mc_5->Fill( ptZ2,newweight*mc_weight_un[5]);
      hYZ2_6_mc_5->Fill( Y_Z2,newweight*mc_weight_un[5]);

      hMZ2_6_mc_6->Fill( massZ2,newweight*mc_weight_un[6] );
      hPtZ2_6_mc_6->Fill( ptZ2,newweight*mc_weight_un[6]);
      hYZ2_6_mc_6->Fill( Y_Z2,newweight*mc_weight_un[6]);

      hMZ2_6_mc_7->Fill( massZ2,newweight*mc_weight_un[7] );
      hPtZ2_6_mc_7->Fill( ptZ2,newweight*mc_weight_un[7]);
      hYZ2_6_mc_7->Fill( Y_Z2,newweight*mc_weight_un[7]);

      hMZ2_6_mc_8->Fill( massZ2,newweight*mc_weight_un[8] );
      hPtZ2_6_mc_8->Fill( ptZ2,newweight*mc_weight_un[8]);
      hYZ2_6_mc_8->Fill( Y_Z2,newweight*mc_weight_un[8]);

     if(nbtag_pass==0) N_bjets_w++;
     if(nbtag_pass==1) N_bjets_cut_w++;
     if(nbtag_pass>=2) N_njets_cut_w++;

// qier (pick one/two bot)
     if(nbtag_pass<2) continue;

     if(debug) {cout <<"GOOD JET PT MAX is: "<< GOOD_JET_PT_MAX << endl;
     cout << "The max value of Deltaphi(jet,MET) among all jets in the event is "<<max_dphi_jet_met<< endl;
     cout << "The min value of Deltaphi(jet,MET) among all jets in the event is "<<min_dphi_jet_met<< endl;
     }
     double Mbb = (BOT1+BOT2).M();
     double PTbb = (BOT1+BOT2).Pt();

     double wu1=scaleFactor1_up/scaleFactor1-1;
     double wu2=scaleFactor2_up/scaleFactor2-1;
     double wd1=1-scaleFactor1_dow/scaleFactor1;
     double wd2=1-scaleFactor2_dow/scaleFactor2;
     double w_up=1+sqrt(wu1*wu1+wu2*wu2);
     double w_dow=1-sqrt(wd1*wd1+wd2*wd2);

     double btag_weight=1;

     if(b1&&b2) {Mbb_6_hh->Fill(Mbb,newweight*nbjetweight);
                 Mbb_6_hh_up->Fill(Mbb,newweight*nbjetweight*w_up);
                 Mbb_6_hh_dow->Fill(Mbb,newweight*nbjetweight*w_dow);
                 ptbb_6_hh->Fill(PTbb,newweight*nbjetweight);
                 ptbb_6_hh_up->Fill(PTbb,newweight*nbjetweight*w_up);
                 ptbb_6_hh_dow->Fill(PTbb,newweight*nbjetweight*w_dow);
                }
     else if((!b1)&&(!b2)) {Mbb_6_ll->Fill(Mbb,newweight*nbjetweight);
                            Mbb_6_ll_up->Fill(Mbb,newweight*nbjetweight*w_up);
                            Mbb_6_ll_dow->Fill(Mbb,newweight*nbjetweight*w_dow);
                            ptbb_6_ll->Fill(PTbb,newweight*nbjetweight);
                            ptbb_6_ll_up->Fill(PTbb,newweight*nbjetweight*w_up);
                            ptbb_6_ll_dow->Fill(PTbb,newweight*nbjetweight*w_dow);
                           }
     else {Mbb_6_hl->Fill(Mbb,newweight*nbjetweight);
           Mbb_6_hl_up->Fill(Mbb,newweight*nbjetweight*w_up);
           Mbb_6_hl_dow->Fill(Mbb,newweight*nbjetweight*w_dow);
           ptbb_6_hl->Fill(PTbb,newweight*nbjetweight);
           ptbb_6_hl_up->Fill(PTbb,newweight*nbjetweight*w_up);
           ptbb_6_hl_dow->Fill(PTbb,newweight*nbjetweight*w_dow);
          }

     Mbb_6->Fill(Mbb,newweight);

     Mbb_6_jer_up->Fill((BOT1_jer_up+BOT2_jer_up).M(),newweight*nbjetweight);
     Mbb_6_jer_dow->Fill((BOT1_jer_dow+BOT2_jer_dow).M(),newweight*nbjetweight);

     Mbb_6_up->Fill((BOT1_up+BOT2_up).M(),newweight*nbjetweight);
     Mbb_6_dow->Fill((BOT1_dow+BOT2_dow).M(),newweight*nbjetweight);

     Mbb_6_pdf_up->Fill(Mbb,newweight*nbjetweight*(1+PDF_weighting_un));
     Mbb_6_pdf_dow->Fill(Mbb,newweight*nbjetweight*(1-PDF_weighting_un));

     Mbb_6_pu_up->Fill(Mbb,newweight*nbjetweight*pu_up);
     Mbb_6_pu_dow->Fill(Mbb,newweight*nbjetweight*pu_dow);

     Mbb_6_mc_1->Fill(Mbb,newweight*nbjetweight*mc_weight_un[1]);
     Mbb_6_mc_2->Fill(Mbb,newweight*nbjetweight*mc_weight_un[2]);
     Mbb_6_mc_3->Fill(Mbb,newweight*nbjetweight*mc_weight_un[3]);
     Mbb_6_mc_4->Fill(Mbb,newweight*nbjetweight*mc_weight_un[4]);
     Mbb_6_mc_5->Fill(Mbb,newweight*nbjetweight*mc_weight_un[5]);
     Mbb_6_mc_6->Fill(Mbb,newweight*nbjetweight*mc_weight_un[6]);
     Mbb_6_mc_7->Fill(Mbb,newweight*nbjetweight*mc_weight_un[7]);
     Mbb_6_mc_8->Fill(Mbb,newweight*nbjetweight*mc_weight_un[8]);

     ptbb_6->Fill(PTbb,newweight*nbjetweight);
     ptbb_6_pdf_up->Fill(PTbb,newweight*nbjetweight*(1+PDF_weighting_un));
     ptbb_6_pdf_dow->Fill(PTbb,newweight*nbjetweight*(1-PDF_weighting_un));
     ptbb_6_pu_up->Fill(PTbb,newweight*nbjetweight*pu_up);
     ptbb_6_pu_dow->Fill(PTbb,newweight*nbjetweight*pu_dow);

     ptbb_6_jer_up->Fill((BOT1_jer_up+BOT2_jer_up).Pt(),newweight*nbjetweight);
     ptbb_6_jer_dow->Fill((BOT1_jer_dow+BOT2_jer_dow).Pt(),newweight*nbjetweight);

     ptbb_6_up->Fill((BOT1_up+BOT2_up).Pt(),newweight*nbjetweight);
     ptbb_6_dow->Fill((BOT1_dow+BOT2_dow).Pt(),newweight*nbjetweight);

     ptbb_6_mc_1->Fill(PTbb,newweight*nbjetweight*mc_weight_un[1]);
     ptbb_6_mc_2->Fill(PTbb,newweight*nbjetweight*mc_weight_un[2]);
     ptbb_6_mc_3->Fill(PTbb,newweight*nbjetweight*mc_weight_un[3]);
     ptbb_6_mc_4->Fill(PTbb,newweight*nbjetweight*mc_weight_un[4]);
     ptbb_6_mc_5->Fill(PTbb,newweight*nbjetweight*mc_weight_un[5]);
     ptbb_6_mc_6->Fill(PTbb,newweight*nbjetweight*mc_weight_un[6]);
     ptbb_6_mc_7->Fill(PTbb,newweight*nbjetweight*mc_weight_un[7]);
     ptbb_6_mc_8->Fill(PTbb,newweight*nbjetweight*mc_weight_un[8]);

     newweight=newweight*nbjetweight;

     hMZ1_7->Fill(massZ1,newweight);
     hPtZ1_7->Fill( ptZ1,newweight );
     hYZ1_7->Fill( Y_Z1,newweight );

      hMZ1_7_pdf_up->Fill( massZ1,newweight*(1+PDF_weighting_un));
      hPtZ1_7_pdf_up->Fill( ptZ1,newweight*(1+PDF_weighting_un));
      hYZ1_7_pdf_up->Fill( Y_Z1,newweight*(1+PDF_weighting_un));

      hMZ1_7_pdf_dow->Fill( massZ1,newweight*(1-PDF_weighting_un));
      hPtZ1_7_pdf_dow->Fill( ptZ1,newweight*(1-PDF_weighting_un));
      hYZ1_7_pdf_dow->Fill( Y_Z1,newweight*(1-PDF_weighting_un));

      hMZ1_7_pu_up->Fill( massZ1,newweight*pu_up);
      hPtZ1_7_pu_up->Fill( ptZ1,newweight*pu_up);
      hYZ1_7_pu_up->Fill( Y_Z1,newweight*pu_up);

      hMZ1_7_pu_dow->Fill( massZ1,newweight*pu_dow);
      hPtZ1_7_pu_dow->Fill( ptZ1,newweight*pu_dow);
      hYZ1_7_pu_dow->Fill( Y_Z1,newweight*pu_dow);

      hMZ1_7_mc_1->Fill( massZ1,newweight*mc_weight_un[1] );
      hPtZ1_7_mc_1->Fill( ptZ1,newweight*mc_weight_un[1]);
      hYZ1_7_mc_1->Fill( Y_Z1,newweight*mc_weight_un[1]);

      hMZ1_7_mc_2->Fill( massZ1,newweight*mc_weight_un[2] );
      hPtZ1_7_mc_2->Fill( ptZ1,newweight*mc_weight_un[2]);
      hYZ1_7_mc_2->Fill( Y_Z1,newweight*mc_weight_un[2]);

      hMZ1_7_mc_3->Fill( massZ1,newweight*mc_weight_un[3] );
      hPtZ1_7_mc_3->Fill( ptZ1,newweight*mc_weight_un[3]);
      hYZ1_7_mc_3->Fill( Y_Z1,newweight*mc_weight_un[3]);

      hMZ1_7_mc_4->Fill( massZ1,newweight*mc_weight_un[4] );
      hPtZ1_7_mc_4->Fill( ptZ1,newweight*mc_weight_un[4]);
      hYZ1_7_mc_4->Fill( Y_Z1,newweight*mc_weight_un[4]);

      hMZ1_7_mc_5->Fill( massZ1,newweight*mc_weight_un[5] );
      hPtZ1_7_mc_5->Fill( ptZ1,newweight*mc_weight_un[5]);
      hYZ1_7_mc_5->Fill( Y_Z1,newweight*mc_weight_un[5]);

      hMZ1_7_mc_6->Fill( massZ1,newweight*mc_weight_un[6] );
      hPtZ1_7_mc_6->Fill( ptZ1,newweight*mc_weight_un[6]);
      hYZ1_7_mc_6->Fill( Y_Z1,newweight*mc_weight_un[6]);

      hMZ1_7_mc_7->Fill( massZ1,newweight*mc_weight_un[7] );
      hPtZ1_7_mc_7->Fill( ptZ1,newweight*mc_weight_un[7]);
      hYZ1_7_mc_7->Fill( Y_Z1,newweight*mc_weight_un[7]);

      hMZ1_7_mc_8->Fill( massZ1,newweight*mc_weight_un[8] );
      hPtZ1_7_mc_8->Fill( ptZ1,newweight*mc_weight_un[8]);
      hYZ1_7_mc_8->Fill( Y_Z1,newweight*mc_weight_un[8]);

     hMZ2_7->Fill(massZ2,newweight);
     hPtZ2_7->Fill( ptZ2,newweight );
     hYZ2_7->Fill( Y_Z2,newweight );

      hMZ2_7_pdf_up->Fill( massZ2,newweight*(1+PDF_weighting_un));
      hPtZ2_7_pdf_up->Fill( ptZ2,newweight*(1+PDF_weighting_un));
      hYZ2_7_pdf_up->Fill( Y_Z2,newweight*(1+PDF_weighting_un));

      hMZ2_7_pdf_dow->Fill( massZ2,newweight*(1-PDF_weighting_un));
      hPtZ2_7_pdf_dow->Fill( ptZ2,newweight*(1-PDF_weighting_un));
      hYZ2_7_pdf_dow->Fill( Y_Z2,newweight*(1-PDF_weighting_un));

      hMZ2_7_pu_up->Fill( massZ2,newweight*pu_up);
      hPtZ2_7_pu_up->Fill( ptZ2,newweight*pu_up);
      hYZ2_7_pu_up->Fill( Y_Z2,newweight*pu_up);

      hMZ2_7_pu_dow->Fill( massZ2,newweight*pu_dow);
      hPtZ2_7_pu_dow->Fill( ptZ2,newweight*pu_dow);
      hYZ2_7_pu_dow->Fill( Y_Z2,newweight*pu_dow);

      hMZ2_7_mc_1->Fill( massZ2,newweight*mc_weight_un[1] );
      hPtZ2_7_mc_1->Fill( ptZ2,newweight*mc_weight_un[1]);
      hYZ2_7_mc_1->Fill( Y_Z2,newweight*mc_weight_un[1]);

      hMZ2_7_mc_2->Fill( massZ2,newweight*mc_weight_un[2] );
      hPtZ2_7_mc_2->Fill( ptZ2,newweight*mc_weight_un[2]);
      hYZ2_7_mc_2->Fill( Y_Z2,newweight*mc_weight_un[2]);

      hMZ2_7_mc_3->Fill( massZ2,newweight*mc_weight_un[3] );
      hPtZ2_7_mc_3->Fill( ptZ2,newweight*mc_weight_un[3]);
      hYZ2_7_mc_3->Fill( Y_Z2,newweight*mc_weight_un[3]);

      hMZ2_7_mc_4->Fill( massZ2,newweight*mc_weight_un[4] );
      hPtZ2_7_mc_4->Fill( ptZ2,newweight*mc_weight_un[4]);
      hYZ2_7_mc_4->Fill( Y_Z2,newweight*mc_weight_un[4]);

      hMZ2_7_mc_5->Fill( massZ2,newweight*mc_weight_un[5] );
      hPtZ2_7_mc_5->Fill( ptZ2,newweight*mc_weight_un[5]);
      hYZ2_7_mc_5->Fill( Y_Z2,newweight*mc_weight_un[5]);

      hMZ2_7_mc_6->Fill( massZ2,newweight*mc_weight_un[6] );
      hPtZ2_7_mc_6->Fill( ptZ2,newweight*mc_weight_un[6]);
      hYZ2_7_mc_6->Fill( Y_Z2,newweight*mc_weight_un[6]);

      hMZ2_7_mc_7->Fill( massZ2,newweight*mc_weight_un[7] );
      hPtZ2_7_mc_7->Fill( ptZ2,newweight*mc_weight_un[7]);
      hYZ2_7_mc_7->Fill( Y_Z2,newweight*mc_weight_un[7]);

      hMZ2_7_mc_8->Fill( massZ2,newweight*mc_weight_un[8] );
      hPtZ2_7_mc_8->Fill( ptZ2,newweight*mc_weight_un[8]);
      hYZ2_7_mc_8->Fill( Y_Z2,newweight*mc_weight_un[8]);

     double M2l2b = (Z2P4+BOT1+BOT2).M();
     double M4l2b = (Z1P4+Z2P4+BOT1+BOT2).M();

     if(b1&&b2) {Mllllbb_hh->Fill(M4l2b,newweight);
                 Mllllbb_hh_up->Fill(M4l2b,newweight*w_up);
                 Mllllbb_hh_dow->Fill(M4l2b,newweight*w_dow);
                 Mllbb_hh->Fill(M2l2b,newweight);
                 Mllbb_hh_up->Fill(M2l2b,newweight*w_up);
                 Mllbb_hh_dow->Fill(M2l2b,newweight*w_dow);
                 mass2D_hh->Fill(M2l2b,M4l2b,newweight);
                 mass2D_hh_up->Fill(M2l2b,M4l2b,newweight*w_up);
                 mass2D_hh_dow->Fill(M2l2b,M4l2b,newweight*w_dow);
                }
     else if((!b1)&&(!b2)) {Mllllbb_ll->Fill(M4l2b,newweight);
                            Mllllbb_ll_up->Fill(M4l2b,newweight*w_up);
                            Mllllbb_ll_dow->Fill(M4l2b,newweight*w_dow);
                            Mllbb_ll->Fill(M2l2b,newweight);
                            Mllbb_ll_up->Fill(M2l2b,newweight*w_up);
                            Mllbb_ll_dow->Fill(M2l2b,newweight*w_dow);
                            mass2D_ll->Fill(M2l2b,M4l2b,newweight);
                            mass2D_ll_up->Fill(M2l2b,M4l2b,newweight*w_up);
                            mass2D_ll_dow->Fill(M2l2b,M4l2b,newweight*w_dow);
                           }
     else {Mllllbb_hl->Fill(M4l2b,newweight);
           Mllllbb_hl_up->Fill(M4l2b,newweight*w_up);
           Mllllbb_hl_dow->Fill(M4l2b,newweight*w_dow);
           Mllbb_hl->Fill(M2l2b,newweight);
           Mllbb_hl_up->Fill(M2l2b,newweight*w_up);
           Mllbb_hl_dow->Fill(M2l2b,newweight*w_dow);
           mass2D_hl->Fill(M2l2b,M4l2b,newweight);
           mass2D_hl_up->Fill(M2l2b,M4l2b,newweight*w_up);
           mass2D_hl_dow->Fill(M2l2b,M4l2b,newweight*w_dow);
          }

     Mllllbb->Fill(M4l2b,newweight);

     Mllllbb_jer_up->Fill((Z1P4+Z2P4+BOT1_jer_up+BOT2_jer_up).M(),newweight);
     Mllllbb_jer_dow->Fill((Z1P4+Z2P4+BOT1_jer_dow+BOT2_jer_dow).M(),newweight);

     Mllllbb_up->Fill((Z1P4+Z2P4+BOT1_up+BOT2_up).M(),newweight);
     Mllllbb_dow->Fill((Z1P4+Z2P4+BOT1_dow+BOT2_dow).M(),newweight);

     Mllllbb_pdf_up->Fill(M4l2b,newweight*(1+PDF_weighting_un));
     Mllllbb_pdf_dow->Fill(M4l2b,newweight*(1-PDF_weighting_un));

     Mllllbb_pu_up->Fill(M4l2b,newweight*pu_up);
     Mllllbb_pu_dow->Fill(M4l2b,newweight*pu_dow);

     Mllllbb_mc_1->Fill(M4l2b,newweight*mc_weight_un[1]);
     Mllllbb_mc_2->Fill(M4l2b,newweight*mc_weight_un[2]);
     Mllllbb_mc_3->Fill(M4l2b,newweight*mc_weight_un[3]);
     Mllllbb_mc_4->Fill(M4l2b,newweight*mc_weight_un[4]);
     Mllllbb_mc_5->Fill(M4l2b,newweight*mc_weight_un[5]);
     Mllllbb_mc_6->Fill(M4l2b,newweight*mc_weight_un[6]);
     Mllllbb_mc_7->Fill(M4l2b,newweight*mc_weight_un[7]);
     Mllllbb_mc_8->Fill(M4l2b,newweight*mc_weight_un[8]);

     Mllbb->Fill(M2l2b,newweight);

     Mllbb_jer_up->Fill((Z2P4+BOT1_jer_up+BOT2_jer_up).M(),newweight);
     Mllbb_jer_dow->Fill((Z2P4+BOT1_jer_dow+BOT2_jer_dow).M(),newweight);

     Mllbb_up->Fill((Z2P4+BOT1_up+BOT2_up).M(),newweight);
     Mllbb_dow->Fill((Z2P4+BOT1_dow+BOT2_dow).M(),newweight);

     Mllbb_pdf_up->Fill(M2l2b,newweight*(1+PDF_weighting_un));
     Mllbb_pdf_dow->Fill(M2l2b,newweight*(1-PDF_weighting_un));

     Mllbb_pu_up->Fill(M2l2b,newweight*pu_up);
     Mllbb_pu_dow->Fill(M2l2b,newweight*pu_dow);

     Mllbb_mc_1->Fill(M2l2b,newweight*mc_weight_un[1]);
     Mllbb_mc_2->Fill(M2l2b,newweight*mc_weight_un[2]);
     Mllbb_mc_3->Fill(M2l2b,newweight*mc_weight_un[3]);
     Mllbb_mc_4->Fill(M2l2b,newweight*mc_weight_un[4]);
     Mllbb_mc_5->Fill(M2l2b,newweight*mc_weight_un[5]);
     Mllbb_mc_6->Fill(M2l2b,newweight*mc_weight_un[6]);
     Mllbb_mc_7->Fill(M2l2b,newweight*mc_weight_un[7]);
     Mllbb_mc_8->Fill(M2l2b,newweight*mc_weight_un[8]);

     mass2D->Fill(M2l2b,newweight);

     mass2D_jer_up->Fill((Z1P4+BOT1_jer_up+BOT2_jer_up).M(),(Z1P4+Z2P4+BOT1_jer_up+BOT2_jer_up).M(),newweight);
     mass2D_jer_dow->Fill((Z1P4+BOT1_jer_dow+BOT2_jer_dow).M(),(Z1P4+Z2P4+BOT1_jer_dow+BOT2_jer_dow).M(),newweight);

     mass2D_up->Fill((Z1P4+BOT1_up+BOT2_up).M(),(Z1P4+Z2P4+BOT1_up+BOT2_up).M(),newweight);
     mass2D_dow->Fill((Z1P4+BOT1_dow+BOT2_dow).M(),(Z1P4+Z2P4+BOT1_dow+BOT2_dow).M(),newweight);

     mass2D_pdf_up->Fill(M2l2b,M4l2b,newweight*(1+PDF_weighting_un));
     mass2D_pdf_dow->Fill(M2l2b,M4l2b,newweight*(1-PDF_weighting_un));

     mass2D_pu_up->Fill(M2l2b,M4l2b,newweight*pu_up);
     mass2D_pu_dow->Fill(M2l2b,M4l2b,newweight*pu_dow);

     mass2D_mc_1->Fill(M2l2b,M4l2b,newweight*mc_weight_un[1]);
     mass2D_mc_2->Fill(M2l2b,M4l2b,newweight*mc_weight_un[2]);
     mass2D_mc_3->Fill(M2l2b,M4l2b,newweight*mc_weight_un[3]);
     mass2D_mc_4->Fill(M2l2b,M4l2b,newweight*mc_weight_un[4]);
     mass2D_mc_5->Fill(M2l2b,M4l2b,newweight*mc_weight_un[5]);
     mass2D_mc_6->Fill(M2l2b,M4l2b,newweight*mc_weight_un[6]);
     mass2D_mc_7->Fill(M2l2b,M4l2b,newweight*mc_weight_un[7]);
     mass2D_mc_8->Fill(M2l2b,M4l2b,newweight*mc_weight_un[8]);

     f_Z1mass=massZ1;
     f_Z2mass=massZ2;
     f_Nbjets=nbtag_pass;
     f_pfmet=RECO_CORMETMUONS; 
     f_mbb=Mbb;
     f_m2l2b=M2l2b;
     f_m4l2b=M4l2b;
     //Number of jets and mJJ,delta eta cuts // categories

     

     //  exactly 4 leptons + at least 2 jets with Djet>0.5 + at most 1 b-tag jet in the event  - category 2


     // filling branches in the reduced tree
     f_weight = newweight;
 
     f_int_weight = -1;
     
     f_pu_weight = pu_weight;
     f_eff_weight = eff_weight;

     f_run = Run;
     f_event = Event;
     f_lumi = LumiSection;

     if(debug) cout << "test final f_weight=" << newweight << endl;
     if(debug) cout << "eff_weight=" << eff_weight << endl;
     //cout << "era BF" << endl;
 
     if (Z1tag==1){
       f_lept1_pt  = RECOMU_PT[indexleptonfinal[0]] ;
       f_lept1_eta = RECOMU_ETA[indexleptonfinal[0]] ;
       f_lept1_phi = RECOMU_PHI[indexleptonfinal[0]];
       f_lept1_charge = RECOMU_CHARGE[indexleptonfinal[0]];
       f_lept1_pfx = RECOMU_PFX_dB_new[indexleptonfinal[0]];
       f_lept1_sip = RECOMU_SIP[indexleptonfinal[0]];
       //    f_lept1_mvaid = RECOMU_mvaNonTrigV0[indexleptonfinal[0]];
       f_lept2_pt  = RECOMU_PT[indexleptonfinal[1]] ;
       f_lept2_eta = RECOMU_ETA[indexleptonfinal[1]] ;
       f_lept2_phi = RECOMU_PHI[indexleptonfinal[1]];
       f_lept2_charge = RECOMU_CHARGE[indexleptonfinal[1]];
       f_lept2_pfx = RECOMU_PFX_dB_new[indexleptonfinal[1]];
       f_lept2_sip = RECOMU_SIP[indexleptonfinal[1]];
       //    f_lept2_mvaid = RECOMU_mvaNonTrigV0[indexleptonfinal[1]];
     }
     else if (Z1tag==2) {
       f_lept1_pt = RECOELE_PT[indexleptonfinal[0]] ;
       f_lept1_eta = RECOELE_ETA[indexleptonfinal[0]] ;
       f_lept1_phi = RECOELE_PHI[indexleptonfinal[0]];
       f_lept1_charge = RECOELE_CHARGE[indexleptonfinal[0]];
       f_lept1_pfx = RECOELE_PFX_rho_new[indexleptonfinal[0]];
       f_lept1_sip = RECOELE_SIP[indexleptonfinal[0]];
       f_lept1_mvaid = RECOELE_mvaNonTrigV0[indexleptonfinal[0]];
       f_lept2_pt = RECOELE_PT[indexleptonfinal[1]] ;
       f_lept2_eta = RECOELE_ETA[indexleptonfinal[1]] ;
       f_lept2_phi = RECOELE_PHI[indexleptonfinal[1]];
       f_lept2_charge = RECOELE_CHARGE[indexleptonfinal[1]];
       f_lept2_pfx = RECOELE_PFX_rho_new[indexleptonfinal[1]];
       f_lept2_sip = RECOELE_SIP[indexleptonfinal[1]];
       f_lept2_mvaid = RECOELE_mvaNonTrigV0[indexleptonfinal[1]];

     }
     
     //f_iso_max = Iso_max;
     //f_sip_max = Sip_max;
     f_angle_costhetastar = -99.;
     f_angle_costheta1 = -99.;
     f_angle_costheta2 = -99.;
     f_angle_phi = -99.;
     f_angle_phistar1 = -99.;
     f_njets_pass = njets_pass;
     f_deltajj = -999.;
     f_massjj = -999.;
     f_D_jet = -999.;
     
     
     ++N_8; 
     // fill final tree
//     finaltree->Fill();
     newtree->Fill();
     if(jentry%5000 == 0) cout << "filling tree" << endl;
     

   } // end loop on entries

   // write on output txt file:

   

   cout << "N_0 "  << N_0  << " \n" 
	      << "N_01 " << N_01 << " \n"	
	      << "N_02 " << N_02 << " \n"	
	      << "N_1 "  << N_1  << " \n"	
	      << "N_2 "  << N_2  << " \n"	
	      << "N_3a " << N_3a << " \n"	
	      << "N_3_FSR " << N_3_FSR << " \n"	
	      << "N_3b " << N_3b << " \n"	
	      << "N_4a " << N_4a << " \n"	
	      << "N_4b " << N_4b << " \n"	
	      << "N_4c " << N_4c << " \n"	
	      << "N_4d " << N_4d << " \n"	
	      << "N_5 "  << N_5  << " \n"	
	      << "N_6 "  << N_6  << " \n"	
	      << "N_7 "  << N_7  << " \n"	
	      << "bjet=2 "  << N_njets_cut_w  << " \n"
	      << "bjet=0 "<< N_bjets_w<< " \n"
      	      << "bjet=1 "  << N_bjets_cut_w  << " \n";

   
   nEvent_4l->GetXaxis()->SetBinLabel(1,"Init.");
   nEvent_4l->GetXaxis()->SetBinLabel(2,"MCTruth: 4mu");
   nEvent_4l->GetXaxis()->SetBinLabel(3,"MCTruth: Acc");
   nEvent_4l->GetXaxis()->SetBinLabel(4,"Init");
   nEvent_4l->GetXaxis()->SetBinLabel(5,"At least 2 tight lept.");
   nEvent_4l->GetXaxis()->SetBinLabel(6,"Z+#gamma");
   nEvent_4l->GetXaxis()->SetBinLabel(7,"Z candidate");
   nEvent_4l->GetXaxis()->SetBinLabel(8,"m_{Z}");
   nEvent_4l->GetXaxis()->SetBinLabel(9,"Z_{1} and Z_{2}");
   nEvent_4l->GetXaxis()->SetBinLabel(10,"Z_{2} pairs > 1");
   nEvent_4l->GetXaxis()->SetBinLabel(11,"m4l > 70");
   nEvent_4l->GetXaxis()->SetBinLabel(12,"MELA KD > 0.1 and mH>100");
   nEvent_4l->GetXaxis()->SetBinLabel(13,"MET>100");
   nEvent_4l->GetXaxis()->SetBinLabel(14,"one Z+#gamma");
   nEvent_4l->GetXaxis()->SetBinLabel(15,"two Z+#gamma");

   nEvent_4l_w->GetXaxis()->SetBinLabel(1,"Init.");
   nEvent_4l_w->GetXaxis()->SetBinLabel(2,"MCTruth: 4mu");
   nEvent_4l_w->GetXaxis()->SetBinLabel(3,"MCTruth: Acc");
   nEvent_4l_w->GetXaxis()->SetBinLabel(4,"Init");
   nEvent_4l_w->GetXaxis()->SetBinLabel(5,"At least 2 tight lept.");
   nEvent_4l_w->GetXaxis()->SetBinLabel(6,"Z+#gamma");
   nEvent_4l_w->GetXaxis()->SetBinLabel(7,"Z candidate");
   nEvent_4l_w->GetXaxis()->SetBinLabel(8,"m_{Z}");
   nEvent_4l_w->GetXaxis()->SetBinLabel(9,"Z_{1} and Z_{2}");
   nEvent_4l_w->GetXaxis()->SetBinLabel(10,"Z_{2} pairs > 1");
   nEvent_4l_w->GetXaxis()->SetBinLabel(11,"m4l > 70");
   nEvent_4l_w->GetXaxis()->SetBinLabel(12,"MELA KD > 0.1 and mH>100");
   nEvent_4l_w->GetXaxis()->SetBinLabel(13,"MET>100");
   nEvent_4l_w->GetXaxis()->SetBinLabel(14,"one Z+#gamma");
   nEvent_4l_w->GetXaxis()->SetBinLabel(15,"two Z+#gamma");

   nEvent_4l->SetBinContent(1,N_0);
   nEvent_4l->SetBinContent(2,N_01);
   nEvent_4l->SetBinContent(3,N_02);
   nEvent_4l->SetBinContent(4,N_1);
   nEvent_4l->SetBinContent(5,N_2);
   nEvent_4l->SetBinContent(6,N_3_FSR);
   nEvent_4l->SetBinContent(7,N_3a);
   nEvent_4l->SetBinContent(8,N_3b);
   nEvent_4l->SetBinContent(9,N_4b);
   nEvent_4l->SetBinContent(10,N_4c);
   nEvent_4l->SetBinContent(11,N_7);
   nEvent_4l->SetBinContent(12,N_9);
   nEvent_4l->SetBinContent(13,N_9_PFMET);
   nEvent_4l->SetBinContent(14,N_9_1FSR);
   nEvent_4l->SetBinContent(15,N_9_2FSR);nEvent_4l_w->SetBinContent(1,N_0_w);
   nEvent_4l_w->SetBinContent(2,N_01_w);
   nEvent_4l_w->SetBinContent(3,N_02_w);
   nEvent_4l_w->SetBinContent(4,N_1_w);
   nEvent_4l_w->SetBinContent(5,N_2_w);
   nEvent_4l_w->SetBinContent(6,N_3_FSR_w);
   nEvent_4l_w->SetBinContent(7,N_3a_w);
   nEvent_4l_w->SetBinContent(8,N_3b_w);
   nEvent_4l_w->SetBinContent(9,N_4b_w);
   nEvent_4l_w->SetBinContent(10,N_4c_w);
   nEvent_4l_w->SetBinContent(11,N_7_w);
   nEvent_4l_w->SetBinContent(12,N_9_w);
   nEvent_4l_w->SetBinContent(13,N_9_PFMET_w);
   nEvent_4l_w->SetBinContent(14,N_9_1FSR_w);
   nEvent_4l_w->SetBinContent(15,N_9_2FSR_w);

   // for(int i=0; i<61; i++){
   //   Char_t PFMET_cut_char[500];
   //   sprintf (PFMET_cut_char,"%s",cut_PFMET.at(i));
   //   nEvent_CUT->GetXaxis()->SetBinLabel(i+1,PFMET_cut_char);
   //   //cout<<PFMET_cut_char<<endl;
   // }

   // for(int i=0; i<61; i++){
   //   Char_t PFMET_cut_char[500];
   //   sprintf (PFMET_cut_char,"%s",cut_PFMET.at(i));
   //   nEvent_CUT_w->GetXaxis()->SetBinLabel(i+1,PFMET_cut_char);
   //   //cout<<PFMET_cut_char<<endl;
   // }


   // write on output root file:
   _filePU->Close();
   theFile->cd();
   //z1tree->Write();
   newtree->Write(0,TObject::kOverwrite);
   theFile->Write();
   theFile->Close();
  // finaltree->Write();
} // end main


double HZZ4LeptonsAnalysis::EAele(int index,bool use2011EA){
  
  double EffectiveArea=0.;
  if (use2011EA){
    if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.18;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.20;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.15;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.19;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.21;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.22;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.4 )                                           EffectiveArea = 0.29;
  }
  //else { // 7_4_x use eta
  // if (fabs(RECOELE_ETA[index]) >= 0.0   && fabs(RECOELE_ETA[index]) < 0.8 )   EffectiveArea = 0.1830;
  // if (fabs(RECOELE_ETA[index]) >= 0.8   && fabs(RECOELE_ETA[index]) < 1.3 )   EffectiveArea = 0.1734;
  // if (fabs(RECOELE_ETA[index]) >= 1.3   && fabs(RECOELE_ETA[index]) < 2.0 )   EffectiveArea = 0.1077;
  // if (fabs(RECOELE_ETA[index]) >= 2.0   && fabs(RECOELE_ETA[index]) < 2.2 )   EffectiveArea = 0.1565;
  // if (fabs(RECOELE_ETA[index]) >= 2.2 )                                       EffectiveArea = 0.2680;
    //}                                                                                                                                                                             
  else { // 7_6_X use eta supercluster                                                                                                                                             
    if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
  }
  
  return EffectiveArea;
  
}

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

double invmass (float M1, float PT1, float ETA1, float PHI1, float M2, float PT2, float ETA2, float PHI2 ){ 
 float phi1=PHI1; 
 float eta1=ETA1; 
 float pt1=PT1; 
 float m1=M1; 

 float px1=pt1*cos(phi1); 
 float py1=pt1*sin(phi1); 
 float pz1=pt1/(2.*(exp(-1*eta1))/(1.0-exp(-2.*eta1))); 

 float phi2=PHI2; 
 float eta2=ETA2; 
 float pt2=PT2; 
 float m2=M2;

 float px2=pt2*cos(phi2); 
 float py2=pt2*sin(phi2); 
 float pz2=pt2/(2.*(exp(-1*eta2))/(1.0-exp(-2.*eta2))); 

 float e1sqr=pz1*pz1+pt1*pt1+m1*m1; 
 float e2sqr=pz2*pz2+pt2*pt2+m2*m2; 
 float e1e2=sqrt(e1sqr*e2sqr); 
 float p1dotp2=px1*px2+py1*py2+pz1*pz2; 

 float m=sqrt(m1*m1+m2*m2+2.*(e1e2-p1dotp2)); 
 //cout << "Invariant mass= " << m << endl;  
 return m;
} // float invmass closed

double HZZ4LeptonsAnalysis::masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr){ 

  //	if(Lep.size()!=4 or pterr.size()!=4) {std::cout<<" Lepsize="<<Lep.size()<<", "<<pterr.size()<<std::endl;}
  int debug_ = 0;
  TLorentzVector compositeParticle ;
  for(unsigned int i=0; i<Lep.size(); i++){
    compositeParticle+=Lep[i];
    if(debug_) std::cout<<" in mass error :  add lep  "<<i<<endl;
  }
  double mass  =  compositeParticle.M();
  
  if(debug_) std::cout<<" in mass error :  mass "<<mass<<endl;
  double masserr = 0;
  
  for(unsigned int i=0; i<Lep.size(); i++){
    if(debug_) std::cout<<" in mass error :  varying lep "<<i<<endl;
    TLorentzVector variedLep; // = Lep[i];
    
    if(debug_) std::cout<<" in mass error : pterr = "<<pterr[i]<<endl;
    variedLep.SetPtEtaPhiM(Lep[i].Pt()+ pterr[i], Lep[i].Eta(), Lep[i].Phi(), Lep[i].M());
    TLorentzVector compositeParticleVariation ;
    for(unsigned int j=0; j<Lep.size(); j++){
      if(i!=j)compositeParticleVariation+=Lep[j];
      else compositeParticleVariation+=variedLep;
    }
    
    masserr += (compositeParticleVariation.M()-mass)*(compositeParticleVariation.M()-mass);
    if(debug_) std::cout<<" in mass error :  intermediate masserr "<<masserr<<endl;
  }
  return sqrt(masserr);
}

float HZZ4LeptonsAnalysis::kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {        
        k+=1.515838921760*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
        k+=1.496256665410*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
        k+=1.495522061910*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
        k+=1.483273154250*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
        k+=1.465589701130*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
        k+=1.491500887510*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
        k+=1.441183580450*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
        k+=1.440830603990*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
        k+=1.414339019120*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
        k+=1.422534218560*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
        k+=1.401037066000*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
        k+=1.408539428810*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
        k+=1.381247744080*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
        k+=1.370553357430*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
        k+=1.347323316000*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
        k+=1.340113437450*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
        k+=1.312661036510*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
        k+=1.290055062010*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
        k+=1.255322614790*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
        k+=1.254455642450*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
        k+=1.224047664420*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
        k+=1.178816782670*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
        k+=1.162624827140*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
        k+=1.105401140940*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
        k+=1.074749265690*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
        k+=1.021864599380*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
        k+=0.946334793286*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
        k+=0.857458082628*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
        k+=0.716607670482*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
        k+=1.132841784840*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);
    }

    if (finalState==2) {
       k+=1.513834489150*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
       k+=1.541738780180*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
       k+=1.497829632510*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
       k+=1.534956782920*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
       k+=1.478217033060*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
       k+=1.504330859290*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
       k+=1.520626246850*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
       k+=1.507013090030*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
       k+=1.494243156250*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
       k+=1.450536096150*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
       k+=1.460812521660*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
       k+=1.471603622200*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
       k+=1.467700038200*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
       k+=1.422408690640*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
       k+=1.397184022730*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
       k+=1.375593447520*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
       k+=1.391901318370*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
       k+=1.368564350560*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
       k+=1.317884804290*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
       k+=1.314019950800*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
       k+=1.274641749910*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
       k+=1.242346606820*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
       k+=1.244727403840*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
       k+=1.146259351670*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
       k+=1.107804993520*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
       k+=1.042053646740*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
       k+=0.973608545141*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
       k+=0.872169942668*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
       k+=0.734505279177*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
       k+=1.163152837230*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);       
    }
    if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
    else return k;

}

float HZZ4LeptonsAnalysis::kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=1.23613311013*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.17550314639*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.17044565911*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.03141209689*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.05285574912*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.11287217794*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.13361441158*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10355603327*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10053981637*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10972676811*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.12069120525*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.11589101635*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.13906170314*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.14854594271*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.14616229031*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.14573157789*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.13829430515*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.15521193686*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.13679822698*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.13223956942*(abs(GENmassZZ)>475.0);
    }

    if (finalState==2) {
        k+=1.25094466582*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.22459455362*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.19287368979*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.04597506451*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.08323413771*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.09994968030*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.16698455800*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10399053155*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10592664340*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10690381480*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.11194928918*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.13522586553*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.11895090244*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.13898508615*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.15463977506*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.17341664594*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.20093349763*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.18915554919*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.18546007375*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.12864505708*(abs(GENmassZZ)>475.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}

float HZZ4LeptonsAnalysis::kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=0.64155491983*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
        k+=1.09985240531*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
        k+=1.29390628654*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
        k+=1.37859998571*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
        k+=1.42430263312*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
        k+=1.45038493266*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
        k+=1.47015377651*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
        k+=1.48828685748*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
        k+=1.50573440448*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
        k+=1.50211655928*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
        k+=1.50918720827*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
        k+=1.52463089491*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
        k+=1.52400838378*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
        k+=1.52418067701*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
        k+=1.55424382578*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
        k+=1.52544284222*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
        k+=1.57896384602*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
        k+=1.53034682567*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
        k+=1.56147329708*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
        k+=1.54468169268*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
        k+=1.57222952415*(abs(GENpTZZ)>100.0);
    }

    if (finalState==2) {
        k+=0.743602533303*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
        k+=1.14789453219*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
        k+=1.33815867892*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
        k+=1.41420044104*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
        k+=1.45511318916*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
        k+=1.47569225244*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
        k+=1.49053003693*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
        k+=1.50622827695*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
        k+=1.50328889799*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
        k+=1.52186945281*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
        k+=1.52043468754*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
        k+=1.53977869986*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
        k+=1.53491994434*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
        k+=1.51772882172*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
        k+=1.54494489131*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
        k+=1.57762411697*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
        k+=1.55078339014*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
        k+=1.57078191891*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
        k+=1.56162666568*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
        k+=1.54183774627*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
        k+=1.58485762205*(abs(GENpTZZ)>100.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}


