#define HZZ4LeptonsAnalysis_cxx
#include "HZZ4LeptonsAnalysis_llbb.h"
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
     cout << "era BF" << endl;
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

   if(era=="BH"){
     pufile="pileup.root";
     puhist="pileup_scale";
     ele_leg1="Leg1_GH_EGM2D.root";
     ele_leg2="Leg2_GH_EGM2D.root";
     id_sf="IDSF_GH.root";
     iso_sf="ISOSF_GH.root";
     mu17_leg="Mu17Leg_SF_GH.root";
     mu8_leg="Mu8Leg_SF_GH.root";
     mu_track="track_GH.root";
     btagcal="CSVv2_Moriond17_B_H.csv";
   }


   RoccoR  rc("/uscms/home/zwang4/nobackup/WORKSPCACE/ntuple/CMSSW_8_0_24/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros/roccor/rcdata.2016.v3"); //directory path as input for now; initialize only once, contains all variations

// setup calibration + reader
   BTagCalibration calib("CSVv2", "btag/"+btagcal);
   BTagCalibrationReader reader(BTagEntry::OP_LOOSE,"central",{"up","down"});      // other sys types
   reader.load(calib,BTagEntry::FLAV_B,"comb");               // measurement type
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

//  TFile *mu_scale_factors1_p1 = new TFile("SF_GH/IDSF_BCDEF.root");
//  TH2F *mu_scale_factors_id_p1 = (TH2F*)gDirectory->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile *mu_scale_factors2_p1 = new TFile("SF_GH/"+iso_sf);
  TH2F *mu_scale_factors_iso_p1 = (TH2F*)gDirectory->Get("LooseISO_LooseID_pt_eta/abseta_pt_ratio");

  TFile *mu_scale_factors3_p1 = new TFile("SF_GH/dm2/"+mu17_leg);
  TH2F *mu_scale_factors_hlt_p1 = (TH2F*)gDirectory->Get("abseta_pt_PLOT");

  TFile *mu_scale_factors4 = new TFile("SF_GH/"+mu_track); //just for GH
  TGraph *mu_scale_factors_tk = (TGraph*)gDirectory->Get("ratio_eff_eta3_dr030e030_corr");

  TFile *b_scale_factors3_p1 = new TFile("btag/loose_eff_all.root");
  TH1F *b_eff_p1 = (TH1F*)gDirectory->Get("hPtBot_8_b");
  TH1F *b_eff_p2 = (TH1F*)gDirectory->Get("hPtBot_8_c");
  TH1F *b_eff_p3 = (TH1F*)gDirectory->Get("hPtBot_8_l");
  TH1F *b_eff_p4 = (TH1F*)gDirectory->Get("hPtBot_8_o");
 
   // Book root file (for output):
   TFile * theFile = new TFile(output,"RECREATE");

   //TString Vars("Weight:Run:Event:LumiSection:massZ1:massZ2:mass4l:Iso_max:Sip_max:MELA:FSR");
   //TNtuple * thePlots=new TNtuple("Candidates","Candidates",Vars);

    // Clone tree for Z1
   //TTree *z1tree = fChain->CloneTree(0);
   

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
   double min_cut_PFMET = 0.; //initialize the min value of RECO_PFMET
   double max_cut_PFMET = 300.;//initialize the max value of RECO_PFMET
   vector <double> cut_PFMET;//declaration of a vector in which the cut values could be stored
   vector <double> counter_cut_PFMET;
   vector <double> counter_cut_PFMET_w;
   
   int cut_n = 61;//number of cuts on PFMET
   double step = (max_cut_PFMET - min_cut_PFMET) / (cut_n - 1);// width between two next cuts
   
   double cut_var_PFMET=0.;
   
   while(cut_var_PFMET <= max_cut_PFMET){
     cut_PFMET.push_back(cut_var_PFMET);//filling the vector with cut's value 
     cut_var_PFMET = cut_var_PFMET + step;//cuts
     
   }
   
   cout << "The size of the vector cut_PFMET is: "<<cut_PFMET.size()<<endl;
   
   //DEBUG
   // for(int i=0; i< cut_PFMET.size(); i++){
   //   cout <<"The value of the cut_PFMET vector at index "<< i <<" is "<<cut_PFMET.at(i)<<endl;
   // }
   //END DEBUG
   
   
   for(int i = 0; i < cut_PFMET.size(); i++){//initialize the counter
     counter_cut_PFMET.push_back(0);
   }
   
   for(int i=0; i< cut_PFMET.size(); i++){
     cout <<"The value of the counter_cut_PFMET vector BEFORE THE LOOP at index "<< i <<" is "<<counter_cut_PFMET.at(i)<<endl;
   }
   
   for(int i = 0; i < cut_PFMET.size(); i++){//initialize the counter
     counter_cut_PFMET_w.push_back(0);
   }
   
   //////////// END OF SENSITIVITY ON MET


   TH1F *hPUvertices             = new TH1F("hPUvertices", "hPUvertices",70,0.,70.);  
   TH1F *hPUvertices_ReWeighted  = new TH1F("hPUvertices_ReWeighted", "hPUvertices_ReWeighted",70,0.,70.);  

   TH1F * hPtJet_7 = new TH1F("hPtJet_7", "Pt of (no ID)jet after selection step 5", 300 ,  0 , 600 );
   hPtJet_7->SetXTitle("pt_jet  (GeV)");

   TH1F * hPtBot_7 = new TH1F("hPtBot_7", "Pt of bot after selection step 5", 300 ,  0 , 600 );
   hPtBot_7->SetXTitle("pt_bot  (GeV)");


   TH1F * hPtJet_7_b = new TH1F("hPtJet_7_b", "Pt of jet_b after selection step 5", 300 ,  0 , 600 );
   hPtJet_7_b->SetXTitle("pt_jet  (GeV)");

   TH1F * hPtBot_7_b = new TH1F("hPtBot_7_b", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_b->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_7_b_dow = new TH1F("hPtBot_7_b_dow", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_b_dow->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_7_b_up = new TH1F("hPtBot_7_b_up", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_b_up->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtJet_7_c = new TH1F("hPtJet_7_c", "Pt of jet_c after selection step 5", 300 ,  0 , 600 );
   hPtJet_7_c->SetXTitle("pt_jet  (GeV)");

   TH1F * hPtBot_7_c = new TH1F("hPtBot_7_c", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_c->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_7_c_up = new TH1F("hPtBot_7_c_up", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_c_up->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_7_c_dow = new TH1F("hPtBot_7_c_dow", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_c_dow->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtJet_7_l = new TH1F("hPtJet_7_l", "Pt of jet_l after selection step 5", 300 ,  0 , 600 );
   hPtJet_7_l->SetXTitle("pt_jet  (GeV)");

   TH1F * hPtBot_7_l = new TH1F("hPtBot_7_l", "Pt of bot_l after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_l->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_7_l_up = new TH1F("hPtBot_7_l_up", "Pt of bot_l after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_l_up->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_7_l_dow = new TH1F("hPtBot_7_l_dow", "Pt of bot_l after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_l_dow->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtJet_7_o = new TH1F("hPtJet_7_o", "Pt of jet_o after selection step 5", 300 ,  0 , 600 );
   hPtJet_7_o->SetXTitle("pt_jet  (GeV)");

   TH1F * hPtBot_7_o = new TH1F("hPtBot_7_o", "Pt of bot_o after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_o->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_7_o_up = new TH1F("hPtBot_7_o_up", "Pt of bot_o after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_o_up->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_7_o_dow = new TH1F("hPtBot_7_o_dow", "Pt of bot_o after selection step 5", 300 ,  0 , 600 );
   hPtBot_7_o_dow->SetXTitle("pt_bot  (GeV)");


   TH1F * hPtJet_8_b = new TH1F("hPtJet_8_b", "Pt of jet_b after selection step 5", 300 ,  0 , 600 );
   hPtJet_8_b->SetXTitle("pt_jet  (GeV)");

   TH1F * hYJet_8_b = new TH1F("hEtaJet_8_b", "Y of jet_b after selection step 5", 500 , -5. , 5. );
   hYJet_8_b->SetXTitle("Y of Jet8");

   TH1F * hPtBot_8_b = new TH1F("hPtBot_8_b", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_b->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_b_dow = new TH1F("hPtBot_8_b_dow", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_b_dow->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_b_up = new TH1F("hPtBot_8_b_up", "Pt of bot_b after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_b_up->SetXTitle("pt_bot  (GeV)");

   TH1F * hYBot_8_b = new TH1F("hEtaBot_8_b", "Y of bot_b after selection step 5", 500 , -5. , 5. );
   hYBot_8_b->SetXTitle("Y of Bot8");

   TH1F * hPtJet_8_c = new TH1F("hPtJet_8_c", "Pt of jet_c after selection step 5", 300 ,  0 , 600 );
   hPtJet_8_c->SetXTitle("pt_jet  (GeV)");

   TH1F * hYJet_8_c = new TH1F("hEtaJet_8_c", "Y of jet_c after selection step 5", 500 , -5. , 5. );
   hYJet_8_c->SetXTitle("Y of Jet8");

   TH1F * hPtBot_8_c = new TH1F("hPtBot_8_c", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_c->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_c_up = new TH1F("hPtBot_8_c_up", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_c_up->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_c_dow = new TH1F("hPtBot_8_c_dow", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_c_dow->SetXTitle("pt_bot  (GeV)");

   TH1F * hYBot_8_c = new TH1F("hEtaBot_8_c", "Y of bot_c after selection step 5", 500 , -5. , 5. );
   hYBot_8_c->SetXTitle("Y of Bot8");

   TH1F * hPtJet_8_l = new TH1F("hPtJet_8_l", "Pt of jet_l after selection step 5", 300 ,  0 , 600 );
   hPtJet_8_l->SetXTitle("pt_jet  (GeV)");

   TH1F * hYJet_8_l = new TH1F("hEtaJet_8_l", "Y of jet_l after selection step 5", 500 , -5. , 5. );
   hYJet_8_l->SetXTitle("Y of Jet8");

   TH1F * hPtBot_8_l = new TH1F("hPtBot_8_l", "Pt of bot_l after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_l->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_l_up = new TH1F("hPtBot_8_l_up", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_l_up->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_l_dow = new TH1F("hPtBot_8_l_dow", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_l_dow->SetXTitle("pt_bot  (GeV)");

   TH1F * hYBot_8_l = new TH1F("hEtaBot_8_l", "Y of bot_l after selection step 5", 500 , -5. , 5. );
   hYBot_8_l->SetXTitle("Y of Bot8");

   TH1F * hPtJet_8_o = new TH1F("hPtJet_8_o", "Pt of jet_o after selection step 5", 300 ,  0 , 600 );
   hPtJet_8_o->SetXTitle("pt_jet  (GeV)");

   TH1F * hYJet_8_o = new TH1F("hEtaJet_8_o", "Y of jet_o after selection step 5", 500 , -5. , 5. );
   hYJet_8_o->SetXTitle("Y of Jet8");

   TH1F * hPtBot_8_o = new TH1F("hPtBot_8_o", "Pt of bot_o after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_o->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_o_up = new TH1F("hPtBot_8_o_up", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_o_up->SetXTitle("pt_bot  (GeV)");

   TH1F * hPtBot_8_o_dow = new TH1F("hPtBot_8_o_dow", "Pt of bot_c after selection step 5", 300 ,  0 , 600 );
   hPtBot_8_o_dow->SetXTitle("pt_bot  (GeV)");

   TH1F * hYBot_8_o = new TH1F("hEtaBot_8_o", "Y of bot_o after selection step 5", 500 , -5. , 5. );
   hYBot_8_o->SetXTitle("Y of Bot8");

   TH1D * hNbjets_8 = new TH1D("hNbjets", "Number of b jets", 10, -0.5, 9.5);
   hNbjets_8->SetXTitle("# b-jets");

   TH1D * hNbjets_8_btag_up = new TH1D("hNbjets_btag_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_btag_dow = new TH1D("hNbjets_btag_dow", "Number of b jets", 10, -0.5, 9.5);

   TH1D * hNbjets_8_h = new TH1D("hNbjets_h", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_h_up = new TH1D("hNbjets_h_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_h_dow = new TH1D("hNbjets_h_dow", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_l = new TH1D("hNbjets_l", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_l_up = new TH1D("hNbjets_l_up", "Number of b jets", 10, -0.5, 9.5);
   TH1D * hNbjets_8_l_dow = new TH1D("hNbjets_l_dow", "Number of b jets", 10, -0.5, 9.5);

   TH1D * hNjets_9 = new TH1D("h_Nbjets_9", "Number of b jets botweights", 10, -0.5, 9.5);
   hNjets_9->SetXTitle("# jets");

   TH1D * hNjets_8 = new TH1D("hNjets_8", "Number of jets passing VBF", 10, -0.5, 9.5);
   hNjets_8->SetXTitle("# jets");
   
   TH1F * Mbb_6 = new TH1F("Mbb_6","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hh = new TH1F("Mbb_6_hh","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hl = new TH1F("Mbb_6_hl","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_ll = new TH1F("Mbb_6_ll","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hh_up = new TH1F("Mbb_6_hh_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hl_up = new TH1F("Mbb_6_hl_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_ll_up = new TH1F("Mbb_6_ll_up","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hh_dow = new TH1F("Mbb_6_hh_dow","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_hl_dow = new TH1F("Mbb_6_hl_dow","invariant mass of bottom pair after step 6",50,20,420);
   TH1F * Mbb_6_ll_dow = new TH1F("Mbb_6_ll_dow","invariant mass of bottom pair after step 6",50,20,420);

   TH1F * ptbb_6 = new TH1F("ptbb_6","pt of bottom pair after step 6",50,20,420);
   ptbb_6->SetXTitle("pt_{bb} (GeV)");
   TH1F * ptbb_6_hh = new TH1F("ptbb_6_hh","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hl = new TH1F("ptbb_6_hl","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_ll = new TH1F("ptbb_6_ll","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hh_up = new TH1F("ptbb_6_hh_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hl_up = new TH1F("ptbb_6_hl_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_ll_up = new TH1F("ptbb_6_ll_up","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hh_dow = new TH1F("ptbb_6_hh_dow","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_hl_dow = new TH1F("ptbb_6_hl_dow","pt of bottom pair after step 6",50,20,420);
   TH1F * ptbb_6_ll_dow = new TH1F("ptbb_6_ll_dow","pt of bottom pair after step 6",50,20,420); 
   // Add branches to output rootuple 
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

      if(jentry%1 == 5000) cout << "Analyzing entry: " << jentry << endl;
      

      if( RECO_NMU > 100 ) RECO_NMU = 100;
      if( RECO_NELE > 100 ) RECO_NELE = 100;
      if( RECO_NPFPHOT > 20 ) RECO_NPFPHOT = 20;
      
      bool debug=false;  //debug flag  -- default false
   
      bool zptweight=false;
      bool etaweight=false;
      bool botsf=false;
      if( datasetName.Contains("DYJetsToLL")){
        zptweight=true;
        cout << "DY correction" << endl;
      }
      newweight=weight;
      cout << "Starting weight= " << newweight << endl;
  
      // pileup reweighting 2012 and 2011
      if (DATA_type=="NO" && num_PU_vertices < 0) continue;                                                                                                                                              
      // pileup reweighting 2015
   
      double pu_weight=1.;
      if (MC_type == "Spring16"||MC_type == "Summer16"){
	Int_t binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
	if(debug) cout << " bin x= " << binx << " " << puweight->GetBinContent(binx) << endl;	
	pu_weight=double(puweight->GetBinContent(binx));
	
      }      
       
      
      //if (num_PU_vertices < 0) continue;

      // Changing the weight for pileup
      newweight=weight*pu_weight;
      cout << "Starting weight + pileup = " << newweight << endl;
           
      

      // Weight for MCNLO samples                                                                                      
      if( datasetName.Contains("amcatnlo")) {
        cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting << endl;
        newweight=weight*pu_weight*MC_weighting;
      }

      
      float pFill[11];for(int pf=0;pf<11;pf++)pFill[11]=-999.;

      // ** Step 0:
      // simply number of entries...
      if( debug ) cout << "\n** Step 0: \nAnalyzing entry: " << jentry << " Run: " << Run << " Event: " << Event << " LumiSection: " << LumiSection << endl ;
      ++N_0 ;  // fill counter
      N_0_w=N_0_w+newweight;
      
      // ** Step 0.1:
      // number of 4L (4mu)...
      //trigger requirements (qier)
      if(!de_trig) continue;
     
      bool is2e2mu=false;

      ++N_1 ;  // fill counter
      N_1_w=N_1_w+newweight;
      
      
      // Effective AREA
      bool tag_2011=false;
      
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
       	
 	if(/* ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0) )
	    && RECOMU_mubesttrkType[i]!=2
	    && RECOMU_PT[i] > 5. 
	    && fabs(RECOMU_ETA[i]) < 2.4 
	    && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.*/
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

       //test(qier) mark dR to closest jet
        double drmin = 10;

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
       	
 	if( RECOELE_PT[i] > 10. 
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

       //test(qier) mark dR to closest jet
        double drmin = 10;

	}
	
      }// end loop on electrons
      
      
      for(int e = 0; e < RECO_NELE; ++e)
      	for(int mu = 0; mu < RECO_NMU; ++mu){
	  
	  if(/*(RECOMU_isPFMu[mu] || (RECOMU_isTrackerHighPtMu[mu] && RECOMU_PT[mu] > 200.))
	      && (RECOMU_isGlobalMu[mu] || (RECOMU_isTrackerMu[mu] && RECOMU_numberOfMatches[mu]>0))
	      && RECOMU_mubesttrkType[mu]!=2
	      && RECOMU_PT[mu] > 5. 
	      && fabs(RECOMU_ETA[mu]) < 2.4 
	      && fabs(RECOMU_mubesttrkDxy[mu]) < .5 && fabs(RECOMU_mubesttrkDz[mu]) < 1. 
	      && fabs(RECOMU_SIP[mu])<4. // TightID + SIP cut*/
 //              RECOMU_isMedium[mu] &&//mediumID(qier) 
              ( RECOMU_isGlobalMu[mu] || RECOMU_isTrackerMu[mu] ) && RECOMU_isPFMu[mu] 
              && RECOMU_PT[mu] > 5. 
              && fabs(RECOMU_ETA[mu]) < 2.4 
              && fabs(RECOMU_mubesttrkDxy[mu]) < .5 && fabs(RECOMU_mubesttrkDz[mu]) < 1. //loose muon (qier)
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
			 << "\nfabs( RECOMU_SIP[i] ) " << fabs( RECOMU_SIP[i] )
			 << "\nfabs( RECOMU_mubesttrkDxy[i] ) " << fabs( RECOMU_mubesttrkDxy[i] )
			 << "\nfabs( RECOMU_mubesttrkDz[i] ) " << fabs( RECOMU_mubesttrkDz[i] )
			 << endl ;
	
       	// Tight muons
 	if( /*(RECOMU_isPFMu[i] || (RECOMU_isTrackerHighPtMu[i] && RECOMU_PT[i] > 200.))
	    && ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0))
	    && RECOMU_mubesttrkType[i]!=2	 
	    && RECOMU_PT[i] > 5. 
	    && fabs(RECOMU_ETA[i]) < 2.4 
	    && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.//tight muon */
 //             RECOMU_isMedium[i] &&//mediumID(qier)
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
          if(debug) cout << "calibration weight = " << mcSF << endl;
          RECOMU_PT[i]=RECOMU_PT[i]*mcSF;
          }
        if( MC_type == "NO" && DATA_type == "2016"){
          double dataSF = rc.kScaleDT(Q,Pt,Eta,Phi);
          RECOMU_PT[i]=RECOMU_PT[i]*dataSF;
          if(debug) cout << "calibration weight = " << dataSF << endl;
        }

	  iL[ N_good ] = i ;
	  ++N_good ;	  
	}
      } // end loop on muons
      
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
       	
 	if( RECOELE_PT[i] > 10. && fabs(RECOELE_ETA[i]) < 2.5 );
	  // && RECOELE_gsftrack_expected_inner_hits[i]<=1 ) /* ok */ ;
	else continue ;

        Ne_loose_4++;

	bool BDT_ok = 0; // Spring16 with CMSSW_8_0_x

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
  	    if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;  // loose ID + SIP cut	    
	    double deltaPhi = DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ) ;
	    double deltaEta = fabs( RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]] );
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	    cout << "debug: " << RECOELE_PT[iL_loose_e[e]] << " " << deltaPhi << " " << deltaEta << " " << deltaR << endl;
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
	  if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;  //loose ID + SIP cut
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	  cout << "DeltaR= " << deltaR << " " <<  deltaR/pow(RECOPFPHOT_PT[iLp[i]],2) << endl;
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    cout << "Possible candidate of photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " associated to a muon with pT= " << RECOMU_PT[iL_loose_mu[l]]<< endl;
	    min_deltaR = deltaR;
	    l_min_deltaR = l;
	    tag_min_deltaR = 0;
	  }
	  
	}//end loop on muons  
	
	for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
	  if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue;  //loose ID + SIP cut
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[l]],2) );
	  cout << "DeltaR= " << deltaR << " " <<  deltaR/pow(RECOPFPHOT_PT[iLp[i]],2) << endl;
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    cout << "Possible candidate of photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " associated to an electron with pT= " << RECOELE_PT[iL_loose_e[l]]<< endl;
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
	if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue; //loose ID + SIP cut
	min_deltaR_ET2=1000;
	p_min_deltaR_ET2=-1;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	    cout <<  "index muon" << iL_loose_mu[l] << endl;
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
	if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue; //loose ID + SIP cut
	min_deltaR_ET2=1000;
	p_min_deltaR_ET2=-1;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	    cout <<  "index electron" << iL_loose_e[l] << endl;
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOELE_ETA[iL_loose_e[l]],2));
	    double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	    cout << " deltaR_ET2= " << deltaR_ET2 <<endl;
	    if (deltaR_ET2<min_deltaR_ET2){
	      min_deltaR_ET2=deltaR_ET2;
	      RECOPFPHOT_DR[iLp[p]]=deltaR;
	      p_min_deltaR_ET2=p;
	      cout << " p_min_deltaR_ET2= " << p_min_deltaR_ET2 <<endl;
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
	if (iLp_l[i]!=-1 && iLp_tagEM[i]==0) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[i]] << endl;
	if (iLp_l[i]!=-1 && iLp_tagEM[i]==1) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to an electron with pT= " << RECOELE_PT[iLp_l[i]] << endl;
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
	  cout << "deltaR for photon subtraction= " << deltaR << endl;
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
	  cout << "deltaR for photon subtraction= " << deltaR << endl;
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
 
      if( N_good + Ne_good < 2 ) continue ; 	
      ++N_2 ;  // fill counter
      N_2_w=N_2_w+newweight;

      int Zxx_tag = 0;    // 1: Zmumu  ,  2: Zee

      int i1 = -1; //index of the first lepton (from Z1)
      int j1 = -1; //index of the second lepton (from Z1)
      int pi1 = -1; 
      int pj1 = -1;
      
      bool has_FSR_Z1 = 0;
      TLorentzVector Lepton1,Lepton2,DiLepton,LeptonCorrection;

      //pick 2 hight pt muon(qier)
//      if(N_good>2) N_good=2;
         // just pick highest pt muon to be pair (qier)
      for(int i = 0; i < N_good; ++i){ //1->N_good
          if(i!=0) continue;   //qier only pick highest pt muon
        for(int j = i + 1; j < N_good; ++j){
//	  if (fabs(RECOMU_SIP[iL[i]])>=4.) continue; // SIP cut (remove qier)
//	  if (fabs(RECOMU_SIP[iL[j]])>=4.) continue;
	  if (fabs(RECOMU_PFX_dB_new[iL[i]])>=0.20) continue; // Isolation
	  if (fabs(RECOMU_PFX_dB_new[iL[j]])>=0.20) continue;
     	  
	  if(RECOMU_CHARGE[ iL[j] ] == RECOMU_CHARGE[ iL[i] ]) continue; // opposite charge
//same sign (qier)(don't forget to change it back)
//          if(RECOMU_CHARGE[ iL[j] ] != RECOMU_CHARGE[ iL[i] ]) continue;

	  cout << "\n Pairing muons with pT= " << RECOMU_PT[ iL[i] ] << " and " <<  RECOMU_PT[ iL[j] ] << endl;
		  
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
	      cout << "Attaching a photon to muon and then to the Z" << endl;
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

	      cout << "Mass Z with FSR= "<< massZ << endl;

	    }
	    if( iLp_l[ p ] == iL[j] && iLp_tagEM[ p ] == 0 )  { 
	      cout << "Attaching a photon to muon and then to the Z" << endl;
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

	      cout << "Mass Z with FSR= "<< massZ << endl;

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
	    cout << "No FSR photon attached" << endl;
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

      for(int i = 0; i < Ne_good; ++i){
        if(i!=0) continue;  //qier always pick highest ele
        for(int j = i + 1; j < Ne_good; ++j){
//	  if (fabs(RECOELE_SIP[iLe[i]])>=4.) continue; // SIP cut
//	  if (fabs(RECOELE_SIP[iLe[j]])>=4.) continue;
	  if (fabs(RECOELE_PFX_rho_new[iLe[i]])>=0.20) continue; // Isolation cut
	  if (fabs(RECOELE_PFX_rho_new[iLe[j]])>=0.20) continue;
	  
	  if(RECOELE_CHARGE[ iLe[j] ] == RECOELE_CHARGE[ iLe[i] ]) continue; // opposite charge
          //don't forget to change it back (qier) same sign
//          if(RECOELE_CHARGE[ iLe[j] ] != RECOELE_CHARGE[ iLe[i] ]) continue;
	  cout << "\n Pairing electrons with pT= " << RECOELE_PT[ iLe[i] ] << " and " <<  RECOELE_PT[ iLe[j] ] << endl;
	  
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
	    if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {  // exist a photon associated to a lepton electron
	      
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

	      cout << "Mass Z with FSR= "<< massZ << endl;
		      
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

	      cout << "Mass Z with FSR= "<< massZ << endl;

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
	    cout << "No FSR photon attached" << endl;
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
      
      
      if (Zcandvector.size()<1) {
	cout << "Less than one Z pairs with isolated leptons...exiting" << endl;
	continue; 
      }
      

      ++N_3a ;  // fill counter
      N_3a_w=N_3a_w+newweight;
      

      // Mass cut on Z
      vector<candidateZ> Zcandisolmassvector;
      Zcandisolmassvector.clear();

      for (int index=0; index<Zcandvector.size();index++){
	if (!(Zcandvector.at(index).massvalue > 60 && Zcandvector.at(index).massvalue < 120)) continue;
	cout << "Z passing the 60 < mll < 120 cut with mass= " << Zcandvector.at(index).massvalue<< endl;
	Zcandisolmassvector.push_back(Zcandvector.at(index));
      };
      
      if (Zcandisolmassvector.size()<1) {
	cout << "No Z passing the mass cut"<< endl;
	continue;
      }

      cout << "Number of Z passing the isolation and the 60 < mll < 120 cut is= " << Zcandisolmassvector.size() << endl;

      ++N_3b ;  // fill counter
      N_3b_w=N_3b_w+newweight;

     
      cout << "Starting weight + pileup + LineShape + efficiency= " << newweight << endl;
//      if(debug) cout << "Efficiency Weight for the Z1: " << eff_weight_3 << " Final weight for Z1= " << newweight << endl;

      int issamesign = 0;

      //if( debug ) cout  << "\nStep 4: Number of good leptons: " << N_good << endl;

      int N_Z2_pairs = 0;

      int i2 = -1; //index of the first lepton (from Z1)
      int j2 = -1; //index of the second lepton (from Z1)
      int pi2 = -1; 
      int pj2 = -1; 
      
      bool has_FSR_Z2 = 0;


      
      
      // PT,20/10 for any di-lepton
      vector<candidateZ> firstpTcleanedgoodZ;    
      vector<float> leptonspTcleaned;
      
      for (int l=0;l<Zcandisolmassvector.size();l++){
        leptonspTcleaned.clear();
        leptonspTcleaned.push_back((Zcandisolmassvector.at(l)).pt1);
        leptonspTcleaned.push_back((Zcandisolmassvector.at(l)).pt2);
        std::sort(leptonspTcleaned.rbegin(),leptonspTcleaned.rend());

        if (leptonspTcleaned.at(0)>30. && leptonspTcleaned.at(1)>20.) {
          firstpTcleanedgoodZ.push_back(Zcandisolmassvector.at(l));
        }
      }

      vector<candidateZ> pTcleanedgoodZ;
      pTcleanedgoodZ=firstpTcleanedgoodZ;
            
     if (pTcleanedgoodZ.size()<1) {
        cout << "No Z(PT) passing the mass cut"<< endl;
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
      
      // Choice of Z1 as the closest to the Z mass
      for (int i=0;i<pTcleanedgoodZ.size();++i){
	
	if( fabs(pTcleanedgoodZ.at(i).massvalue - Zmass) < fabs(massZ1 - Zmass) ){
	  
	  massZ1 = pTcleanedgoodZ.at(i).massvalue;
	  indexZ1=i;
	  
	  pxZ1 = pTcleanedgoodZ.at(i).pxZ;
	  pyZ1 = pTcleanedgoodZ.at(i).pyZ;
	  pzZ1 = pTcleanedgoodZ.at(i).pzZ;
	  EZ1  = pTcleanedgoodZ.at(i).EZ;
	  
	  ptZ1 = sqrt( pxZ1*pxZ1 + pyZ1*pyZ1 );
	  sum_ptZ1 = pTcleanedgoodZ.at(i).pt1+pTcleanedgoodZ.at(i).pt2;
	  
	  Y_Z1 = 0.5 * log ( (EZ1 + pzZ1)/(EZ1 - pzZ1) );
	  indexlep1Z1=pTcleanedgoodZ.at(i).ilept1;
	  indexlep2Z1=pTcleanedgoodZ.at(i).ilept2;
	  Z1tag=pTcleanedgoodZ.at(i).tag;
	}
      }
      
      if (massZ1 < 40.) {
	cout << "The mass of Z1 is < 40 GeV...exiting" << endl;
	continue;
      } 
      
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

      
      
      
      ++N_4b ;  // fill counter
      N_4b_w=N_4b_w+newweight;

      
      // **** Step 5:

      
       // Execute Efficiency Reweighting
 
      int z1lept[2]={indexlep1Z1,indexlep2Z1};

      Double_t eff_weight = 1.;
     
      bool BB=false;
      bool EB=false;
      bool EE=false;   
 
      if (Z1tag==1){ 
          if(abs(RECOMU_ETA[ z1lept[0] ])<1.5 && abs(RECOMU_ETA[ z1lept[1] ])<1.5) BB=true;
          else if(abs(RECOMU_ETA[ z1lept[0] ])>1.5 && abs(RECOMU_ETA[ z1lept[1] ])>1.5) EE=true;
          else EB=true;	

        bool nomatch=false;
        bool noleg17=true;
        for(int i = 0; i < 2; ++i){
	  Double_t Pt = RECOMU_PT[ z1lept[i] ]; 
	  Double_t Eta = RECOMU_ETA[ z1lept[i] ]; 
	  if( MC_type == "Spring16" && DATA_type == "NO"){
//	    Int_t biny = mu_scale_2016->GetYaxis()->FindBin(Pt);
//	    Int_t binx = mu_scale_2016->GetXaxis()->FindBin(Eta);
//           cout << "xbin= " << binx <<" ybin="<<  biny << endl;
//            cout << "eff_weight test = " <<mu_scale_2016->GetBinContent(binx,biny)<< endl;
//	    if (mu_scale_2016->GetBinContent(binx,biny)>0.) eff_weight*=mu_scale_2016->GetBinContent(binx,biny); 


          int biny1 = mu_scale_factors_id_p1->GetYaxis()->FindBin(Pt);
          int binx1 = mu_scale_factors_id_p1->GetXaxis()->FindBin(abs(Eta));

          double sf_id=mu_scale_factors_id_p1->GetBinContent(binx1,biny1); //just for GH
         if (sf_id>0.) eff_weight*=sf_id; 
  
          int biny2 = mu_scale_factors_iso_p1->GetYaxis()->FindBin(Pt);
          int binx2 = mu_scale_factors_iso_p1->GetXaxis()->FindBin(abs(Eta)); 

          double sf_iso = mu_scale_factors_iso_p1->GetBinContent(binx2,biny2);  //just for GH
          if (sf_iso>0.)  eff_weight*=sf_iso; 
 
          double tk_sf = mu_scale_factors_tk->Eval(Eta);    
          if(mu_scale_factors_tk->Eval(Eta)>0) eff_weight*=tk_sf;

          cout << " id weight = " << sf_id << "\n iso weight = " << sf_iso << "\n tk weight = " << tk_sf << endl;

          if(RECOMU_dm_MuHLTMatch[ z1lept[i] ]== 2){
            int biny3 = mu_scale_factors_hlt_p1->GetYaxis()->FindBin(Pt);
            int binx3 = mu_scale_factors_hlt_p1->GetXaxis()->FindBin(abs(Eta));

//            int biny32 = mu_scale_factors_hlt_p2->GetYaxis()->FindBin(Pt);
//            int binx32 = mu_scale_factors_hlt_p2->GetXaxis()->FindBin(abs(Eta));
//            double sf_hlt = mu_scale_factors_hlt_p1->GetBinContent(binx3,biny3)*19.666/35.812+mu_scale_factors_hlt_p2->GetBinContent(binx32,biny32)*16.146/35.812;
            double sf_hlt = mu_scale_factors_hlt_p1->GetBinContent(binx3,biny3);
            if (sf_hlt>0.) eff_weight*=sf_hlt;
            cout << "l1trigger matching mu17 leg weight = " << sf_hlt << endl;
           }
         else{
            int biny32 = mu_scale_factors_hlt_p2->GetYaxis()->FindBin(Pt);
            int binx32 = mu_scale_factors_hlt_p2->GetXaxis()->FindBin(abs(Eta));

            double sf_hlt = mu_scale_factors_hlt_p2->GetBinContent(binx32,biny32);
            if (sf_hlt>0.) eff_weight*=sf_hlt;
            cout << "l1trigger matching mu8 leg weight = " << sf_hlt << endl;
          }
	  }
//tigger matching require (qier) 
  //for matching (qier)
          if(RECOMU_dm_MuHLTMatch[ z1lept[i] ]<0) nomatch=true;
          if(RECOMU_dm_MuHLTMatch[ z1lept[i] ]==2) noleg17=false;
         }
        if(nomatch||noleg17) continue;
  //sm  
//        int lead;
//        if(RECOMU_PT[ z1lept[0] ]>=RECOMU_PT[ z1lept[1] ]){lead=z1lept[0];}
//        else lead=z1lept[1]; 
//        if(RECOMU_sm_MuHLTMatch[lead]) nomatch=false;
//        if(nomatch) continue;
//        if( MC_type == "Spring16" && DATA_type == "NO"){
//            int biny33 = mu_scale_factors_hlt->GetYaxis()->FindBin(RECOMU_PT[lead]);
//            int binx33 = mu_scale_factors_hlt->GetXaxis()->FindBin(abs(RECOMU_ETA[lead]));
//            if (mu_scale_factors_hlt->GetBinContent(binx33,biny33)>0.) eff_weight*=mu_scale_factors_hlt->GetBinContent(binx33,biny33);
//            cout << "l1trigger matching weight = " << mu_scale_factors_hlt->GetBinContent(binx33,biny33) << endl;
//        }
      }
      else if (Z1tag==2){
          if(abs(RECOELE_ETA[ z1lept[0] ])<1.5 && abs(RECOELE_ETA[ z1lept[1] ])<1.5) BB=true;
          else if(abs(RECOELE_ETA[ z1lept[0] ])>1.5 && abs(RECOELE_ETA[ z1lept[1] ])>1.5) EE=true;
          else EB=true;

        bool nomatch=false;
        bool noleg1=true;
 
	for(int i = 0; i < 2; ++i){
	  Double_t Pt = RECOELE_PT[ z1lept[i] ]; 
	  Double_t Eta = RECOELE_scl_Eta[ z1lept[i] ]; 
	  
	  if( MC_type == "Spring16" && DATA_type == "NO"){
  
            cout << "Pt= " << Pt << " Eta= " << Eta << endl;
            int biny4 = ele_scale_factors_reco->GetYaxis()->FindBin(Pt);
            int binx4 = ele_scale_factors_reco->GetXaxis()->FindBin(Eta);
            if (ele_scale_factors_reco->GetBinContent(binx4,biny4)>0.) eff_weight*=ele_scale_factors_reco->GetBinContent(binx4,biny4); 
            cout << "ele reco sf = " << ele_scale_factors_reco->GetBinContent(binx4,biny4) << endl;
            
            int biny5 = ele_scale_factors_wp90->GetYaxis()->FindBin(Pt);
            int binx5 = ele_scale_factors_wp90->GetXaxis()->FindBin(Eta);
            if (ele_scale_factors_wp90->GetBinContent(binx5,biny5)>0.) eff_weight*=ele_scale_factors_wp90->GetBinContent(binx5,biny5);
            cout << "ele wp90 sf = " << ele_scale_factors_wp90->GetBinContent(binx5,biny5) << endl;
           
          if(RECOELE_de_EleHLTMatch[ z1lept[i] ]== 2){
            int biny6 = ele_scale_factors_leg1->GetYaxis()->FindBin(Pt);
            int binx6 = ele_scale_factors_leg1->GetXaxis()->FindBin(Eta);

//            int biny32 = mu_scale_factors_hlt_p2->GetYaxis()->FindBin(Pt);
//            int binx32 = mu_scale_factors_hlt_p2->GetXaxis()->FindBin(abs(Eta));
//            double sf_hlt = mu_scale_factors_hlt_p1->GetBinContent(binx3,biny3)*19.666/35.812+mu_scale_factors_hlt_p2->GetBinContent(binx32,biny32)*16.146/35.812;
            double sf_hlt = ele_scale_factors_leg1->GetBinContent(binx6,biny6);
            if (sf_hlt>0.) eff_weight*=sf_hlt;
            cout << "l1trigger matching ele leg1 weight = " << sf_hlt << endl;
           }
         else{
            int biny62 = ele_scale_factors_leg2->GetYaxis()->FindBin(Pt);
            int binx62 = ele_scale_factors_leg2->GetXaxis()->FindBin(Eta);

            double sf_hlt = ele_scale_factors_leg2->GetBinContent(binx62,biny62);
            if (sf_hlt>0.) eff_weight*=sf_hlt;
            cout << "l1trigger matching ele leg2 weight = " << sf_hlt << endl;
          }   
         
	  }
          if(RECOELE_de_EleHLTMatch[ z1lept[i] ]<=0) nomatch=true;
          if(RECOELE_de_EleHLTMatch[ z1lept[i] ]==2) noleg1=false;
	}
        if(nomatch||noleg1) continue;
      }

      cout << "eff_weight" << eff_weight << endl;
      
      // // Changing the weight for pileup and efficiency
      if (eff_weight>0.) newweight=newweight*eff_weight;
      
      cout << "Starting weight + pileup + efficiency= " << newweight << endl;
      if(debug) cout << "Efficiency Weight for the 4l: " << eff_weight << " Final weight= " << newweight << endl;
      
      //only keep 2mu (qier) 
      if(Z1tag!=2) continue;

      TLorentzVector Z1P4;
      Z1P4.SetPxPyPzE(pxZ1,pyZ1,pzZ1,EZ1);

      if( MC_type == "Spring16" && DATA_type == "NO" && zptweight){
        if(ptZ1<5) newweight=newweight*1.04;
        else if(ptZ1<10) newweight=newweight*1.05;
        else if(ptZ1>=10&&ptZ1<15) newweight=newweight*1.00;
        else if(ptZ1>=15&&ptZ1<20) newweight=newweight*0.96;
        else if(ptZ1>=20&&ptZ1<25) newweight=newweight*0.92;
        else if(ptZ1>=25&&ptZ1<30) newweight=newweight*0.91;
        else if(ptZ1>=30&&ptZ1<40) newweight=newweight*0.92;
        else if(ptZ1>=40&&ptZ1<45) newweight=newweight*0.93;
        else if(ptZ1>=45&&ptZ1<50) newweight=newweight*0.95;
        else if(ptZ1>=50&&ptZ1<55) newweight=newweight*0.95;
        else if(ptZ1>=55&&ptZ1<60) newweight=newweight*0.96;
        else if(ptZ1>=60&&ptZ1<65) newweight=newweight*0.98;
        else if(ptZ1>=60&&ptZ1<70) newweight=newweight*0.96;
      }

      if(MC_type == "Spring16" && DATA_type == "NO" && etaweight){
        if(abs(RECOMU_ETA[z1lept[0] ])>1.7) newweight=newweight*1.05;
      }


      // sort index by pt (kinematics not corrected for FSR)
      int ipt[2] ;
      double tmp_pt[2];
      int tmp_type[2];
      int lep_type[2];
  

      int indexleptonfinal[2]={indexlep1Z1,indexlep2Z1};
      //cout << "PTs= " << RECOMU_PT[indexleptonfinal[0]] << " " << RECOMU_PT[indexleptonfinal[1]] << " " <<  RECOMU_PT[indexleptonfinal[2]] << " " << RECOMU_PT[indexleptonfinal[3]]<< endl;

      for(int i = 0; i < 2; ++i){ 
	if (Z1tag==1) {tmp_pt[i] =  RECOMU_PT[indexleptonfinal[i]]; tmp_type[i]=1;}
	if (Z1tag==2) {tmp_pt[i] =  RECOELE_PT[indexleptonfinal[i]]; tmp_type[i]=2;}
	cout << tmp_pt[i] << endl;
      }

      float sortedpT[2];
 
      for(int i = 0; i < 2; ++i){		
        double tmp_max_pt = 0;
      	int jj = i;
        int type = 0;
        for(int j = 0; j < 2; ++j){
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
     

      // Format lepton syncronization                                                                                                                                      
      // {run}:{lumi}:{event}:{pdgId}:{pT:.2f}:{eta:.2f}:{phi:.2f}{SIP:.2f}:{PFChargedHadIso:.2f}:{PFNeutralHadIso:.2f}:{PFPhotonIso:.2f}:{PUCorr:.2f}:{combRelIsoPF:.3f}:{eleBDT:.3f}:{photpT:.2f}:{photDR:.2f}:{photRelIso:.2f}          
      
      Char_t leptformat[20000];

      for(int i = 0; i < 2; ++i){      
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
	  
       
      }//end fill leptons
      
      

      // N.B. Do NOT Update the Isolation values and correct the 4 momenta of leptons for FSR
      for(int i = 0; i < N_good; ++i){
	int flagFSR=0;
	int pfsr=-999;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if (iLp[p]==-1) continue;
	  if (iLp_l[p]==-1) continue;
	  
	  cout << "Index of lepton with photon ISR= " << iLp_l[ p ] << " and final lepton index= " << iL[i] << endl;
	  if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )  {
	    cout << "Muon with pT= " << RECOMU_PT[iL[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	    
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


     
     // Leptons PT, ETA, Phi, Isol corrected for FSR
     
     
     ++N_8_PFMET;
     N_8_PFMET_w=N_8_PFMET_w+newweight;
     
     
     //Basic cuts to jets AND delta R section
     int njets_pass=0;
     int nbtag_pass=0;
     int nhb_pass=0,nlb_pass=0;
     bool lead_jet=false;
     TLorentzVector JET1,JET2,BOT1,BOT2;
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
       
       // JET smearing
       double scaleFactor_1=1;
       double scaleFactor_1_up=1;
       double scaleFactor_1_dow=1;

       double Pt=RECO_PFJET_PT[i];
       double Eta=RECO_PFJET_ETA[i];

       if(RECO_PFJET_PT[i]>10. && fabs(RECO_PFJET_ETA[i])<2.4){
           hPtJet_7->Fill(RECO_PFJET_PT[i],newweight);
              if( MC_type == "Spring16" && DATA_type == "NO"){

           if(RECOBOT_MatchingMCTruth[i+1]==1){
             hPtJet_7_b->Fill(RECO_PFJET_PT[i],newweight);
           }
           else if(RECOBOT_MatchingMCTruth[i+1]==2){
             hPtJet_7_c->Fill(RECO_PFJET_PT[i],newweight);
           }
           else if(RECOBOT_MatchingMCTruth[i+1]==3){
             hPtJet_7_l->Fill(RECO_PFJET_PT[i],newweight);
           }
           else if(RECOBOT_MatchingMCTruth[i+1]==5){
             hPtJet_7_o->Fill(RECO_PFJET_PT[i],newweight);
           }


              if(cSV_BTagJet_DISCR[i]> 0.5426) {

                 if(RECOBOT_MatchingMCTruth[i+1]==2){
                    scaleFactor_1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, Eta, Pt);
                    scaleFactor_1_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_C, Eta, Pt);
                    scaleFactor_1_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_C, Eta, Pt);
                    hPtBot_7_c_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1_up);
                    hPtBot_7_c_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1_dow);
                    hPtBot_7_c->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1);
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==1){
                    scaleFactor_1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, Eta, Pt);
                    scaleFactor_1_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_B, Eta, Pt);
                    scaleFactor_1_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_B, Eta, Pt);
                    hPtBot_7_b_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1_up);
                    hPtBot_7_b_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1_dow);
                    hPtBot_7_b->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1);

                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==3){
                    scaleFactor_1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_1_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_1_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    hPtBot_7_l_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1_up);
                    hPtBot_7_l_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1_dow);
                    hPtBot_7_l->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1);
                 }
                 else if(RECOBOT_MatchingMCTruth[i+1]==5){
                    scaleFactor_1 = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_1_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_1_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    hPtBot_7_o->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1);
                    hPtBot_7_o_dow->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1_dow);
                    hPtBot_7_o_up->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1_up);
                 }

                }

              }
          if(cSV_BTagJet_DISCR[i]> 0.5426) {
          hPtBot_7->Fill(RECO_PFJET_PT[i],newweight*scaleFactor_1);}
       }


       double jercorr = 1.0; double jercorrup = 1.0; double jercorrdn = 1.0;

       bool goodjet = RECO_PFJET_NHF[i] < 0.99 &&
                      RECO_PFJET_NEF[i] < 0.99 &&
                      RECO_PFJET_CHF[i] < 0.99 &&
                      RECO_PFJET_CEF[i] < 0.99 &&
                      RECO_PFJET_nconstituents[i] > 1 &&
                      RECO_PFJET_NCH[i] > 0;
       
 
       if(RECO_PFJET_PT[i]>10. && fabs(RECO_PFJET_ETA[i])<2.4){
       
      	 //Check that jet has deltaR>0.4 away from any tight lepton corrected for FSR
	 for(int mu = 0; mu < N_good; ++mu){
//	   if (fabs(RECOMU_SIP[iL[mu]])>=4.) continue;  //qier
    	   if (RECOMU_PFX_dB_new[iL[mu]]>=0.20) continue;
	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOMU_PHI[iL[mu]]),2) + pow(RECO_PFJET_ETA[i] - RECOMU_ETA[iL[mu]],2));
	   cout << "1st lepton muon: " << " pT=" << RECOMU_PT[iL[mu]] <<" deltaR "<< deltaR <<endl;	   
	   if (deltaR<0.4){
	     jetfail[i]=1;
     	     cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
      	 for(int ele = 0; ele < Ne_good; ++ele){
      	//   if (fabs(RECOELE_SIP[iLe[ele]])>=4.) continue;
	   if (RECOELE_PFX_rho_new[iLe[ele]]>=0.35) continue;
      	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOELE_PHI[iLe[ele]]),2) + pow(RECO_PFJET_ETA[i] - RECOELE_ETA[iLe[ele]],2));
     	   cout << "1st lepton electron: " << " pT=" << RECOELE_PT[iLe[ele]] <<" deltaR "<< deltaR <<endl;
	   if (deltaR<0.4){
     	     jetfail[i]=1;
     	     cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
	 // cleaning w.r.t FSR photons attached to leptons
	 for(int j=0.;j<Nphotons;j++) {
           if (iLp_l[j]!=-1 && (iLp_tagEM[j]==0 || iLp_tagEM[j]==1) ) {
	     if (iLp_tagEM[j]==0) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[j]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[j]] << endl;
	     if (iLp_tagEM[j]==1) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[j]] << " attached to a electron with pT= " << RECOELE_PT[iLp_l[j]] << endl;
	     double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOPFPHOT_PHI[iLp[j]]),2) + pow(RECO_PFJET_ETA[i] - RECOPFPHOT_ETA[iLp[j]],2));
	     if (deltaR<0.4){
	       jetfail[i]=1;
	       cout << " jetfail " << jetfail[i] <<endl;
	       break;
	     }
	   }
         }
	 // 


	 if (jetfail[i]==0){
	   cout<< " PASS jet " <<i<<" PT= "<<RECO_PFJET_PT[i]<<" ETA= "<<RECO_PFJET_ETA[i]<<" PUID= "<<RECO_PFJET_PUID[i]<<endl;
	   njets_pass++;
           if(RECO_PFJET_PT[i]>50) lead_jet = true;
           if(RECOBOT_MatchingMCTruth[i+1]==1){
             hPtJet_8_b->Fill(RECO_PFJET_PT[i],newweight);
             hYJet_8_b->Fill(RECO_PFJET_ETA[i],newweight);
           }
           else if(RECOBOT_MatchingMCTruth[i+1]==2){
             hPtJet_8_c->Fill(RECO_PFJET_PT[i],newweight);
             hYJet_8_c->Fill(RECO_PFJET_ETA[i],newweight);
           }
           else if(RECOBOT_MatchingMCTruth[i+1]==3){
             hPtJet_8_l->Fill(RECO_PFJET_PT[i],newweight);
             hYJet_8_l->Fill(RECO_PFJET_ETA[i],newweight);
           }
           else if(RECOBOT_MatchingMCTruth[i+1]==5){
             hPtJet_8_o->Fill(RECO_PFJET_PT[i],newweight);
             hYJet_8_o->Fill(RECO_PFJET_ETA[i],newweight);
           }

           double scaleFactor=1;
           double scaleFactor_up=1;
           double scaleFactor_dow=1;

           //b-tagging
           if(cSV_BTagJet_DISCR[i]> 0.5426) {
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
                    scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_up = reader.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
                    scaleFactor_dow = reader.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, abs(Eta), Pt);
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
                 if(Eta<-0.6&&Eta>-1.8&&botsf) scaleFactor*=0.91;
                 if(Eta>0.6&&Eta<1.0&&botsf) scaleFactor*=0.95;
              }
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
     hNbjets_8->Fill(nbtag_pass,newweight*nbjetweight);
     hNbjets_8_h->Fill(nhb_pass,newweight*nbjetweight);
     hNbjets_8_h_up->Fill(nhb_pass,newweight*nbjetweight*(1+sqrt(scale_h_up)));
     hNbjets_8_h_dow->Fill(nhb_pass,newweight*nbjetweight*(1-sqrt(scale_h_dow)));
     hNbjets_8_l->Fill(nlb_pass,newweight*nbjetweight);
     hNbjets_8_l_up->Fill(nlb_pass,newweight*nbjetweight*(1+sqrt(scale_l_up)));
     hNbjets_8_l_dow->Fill(nlb_pass,newweight*nbjetweight*(1-sqrt(scale_l_dow)));

     hNjets_9->Fill(njets_pass,newweight*nbjetweight);

     if(nbtag_pass<2) continue;
  
     ++N_8; 
     double Mbb = (BOT1+BOT2).M();
     double PTbb = (BOT1+BOT2).Pt();

     Mbb_6->Fill(Mbb,newweight*nbjetweight);
     ptbb_6->Fill(PTbb,newweight*nbjetweight);
     double wu1=scaleFactor1_up/scaleFactor1-1;
     double wu2=scaleFactor2_up/scaleFactor2-1;
     double wd1=1-scaleFactor1_dow/scaleFactor1;
     double wd2=1-scaleFactor2_dow/scaleFactor2;
     double w_up=sqrt(wu1*wu1+wu2*wu2);//scaleFactor1_up*scaleFactor2_up/(scaleFactor1*scaleFactor2);
     double w_dow=sqrt(wd1*wd1+wd2*wd2);//scaleFactor1_dow*scaleFactor2_dow/(scaleFactor1*scaleFactor2);

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

   } // end loop on entries

   // write on output txt file:

   // write on output root file:
   _filePU->Close();
   theFile->cd();
   theFile->Write();
   theFile->Close();
  // finaltree->Write();
} // end main

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
