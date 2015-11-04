//_____________________________________________________________________________
// Standard Headers:
#include <fstream>
#include <cmath>
#include <iostream>
#include <unistd.h>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include <TVector3.h>
#include "TLorentzVector.h"
#include "TRandom.h"

#include <string.h>
#include "KinLine.C"
#include "Kstream.C"

#include"/Users/michaelkunkel/WORK/GIT_HUB/G12_Corrections/g12_corrections.hpp"


using namespace std;

int tot_events = 0;
int Bad_events = 0;
int Good_events = 0;


Double_t OneBoostComparison(TLorentzVector daughter1, TLorentzVector parent){
  
  daughter1.Boost(-parent.BoostVector());
  
  Double_t CosTheta = daughter1.CosTheta();
  
  if(TMath::IsNaN(CosTheta)){return -1000.;}
  else{return CosTheta;}
  
}

int main(){
  // Initialization of Momentum Correction Parameters
  clog << "Reading in momentum correction parameters...\n";
  clas::g12::MomentumCorrection pcor;
  //extern int optind;
  //gROOT->Reset();
  //TROOT troot();
  
  TChain *chain = new TChain("lepTree");
  chain->Add("IsLepG7_g12.root");
  
  Int_t in_Nphoton;
  Int_t N_in_timePhoton;
  Int_t in_photon_stat[100];
  Float_t in_photon_E[100];
  Float_t in_photon_tpho[100];
  Float_t in_photon_ttag[100];
  
  chain->SetBranchAddress("in_Nphoton",&in_Nphoton);
  chain->SetBranchAddress("N_in_timePhoton",&N_in_timePhoton);
  chain->SetBranchAddress("in_photon_stat",&in_photon_stat);
  chain->SetBranchAddress("in_photon_E",&in_photon_E);
  chain->SetBranchAddress("in_photon_tpho",&in_photon_tpho);
  chain->SetBranchAddress("in_photon_ttag",&in_photon_ttag);
  
  
  Float_t prot_P, prot_Theta, prot_Phi, egam, tpho;
  Float_t pip_P, pip_Theta, pip_Phi;
  Float_t pim_P, pim_Theta, pim_Phi;
  Float_t prot_mat[5][5], pim_mat[5][5], pip_mat[5][5];
  
  Int_t event, run, trigBits;
  chain->SetBranchAddress("event",&event);
  chain->SetBranchAddress("run",&run);
  chain->SetBranchAddress("trigBits",&trigBits);
  
  
  chain->SetBranchAddress("egam",&egam);
  chain->SetBranchAddress("tpho",&tpho);
  
  chain->SetBranchAddress("prot_P",&prot_P);
  chain->SetBranchAddress("prot_Theta",&prot_Theta);
  chain->SetBranchAddress("prot_Phi",&prot_Phi);
  
  chain->SetBranchAddress("pip_P",&pip_P);
  chain->SetBranchAddress("pip_Theta",&pip_Theta);
  chain->SetBranchAddress("pip_Phi",&pip_Phi);
  
  chain->SetBranchAddress("pim_P",&pim_P);
  chain->SetBranchAddress("pim_Theta",&pim_Theta);
  chain->SetBranchAddress("pim_Phi",&pim_Phi);
  
  chain->SetBranchAddress("prot_mat",&prot_mat);
  chain->SetBranchAddress("pim_mat",&pim_mat);
  chain->SetBranchAddress("pip_mat",&pip_mat);
  
  Int_t prot_sec, pip_sec, pim_sec;
  chain->SetBranchAddress("prot_sec",&prot_sec);
  chain->SetBranchAddress("pip_sec",&pip_sec);
  chain->SetBranchAddress("pim_sec",&pim_sec);
  
  Int_t pip_EC_hit, pim_EC_hit;
  chain->SetBranchAddress("pip_EC_hit",&pip_EC_hit);
  chain->SetBranchAddress("pim_EC_hit",&pim_EC_hit);
  
  float pim_ECx;
  float pim_ECy;
  float pim_ECz;
  float pim_ECin, pim_ECout;
  
  chain->SetBranchAddress("pim_ECx",&pim_ECx);
  chain->SetBranchAddress("pim_ECy",&pim_ECy);
  chain->SetBranchAddress("pim_ECz",&pim_ECz);
  chain->SetBranchAddress("pim_ECin",&pim_ECin);
  chain->SetBranchAddress("pim_ECout",&pim_ECout);
  
  
  float pip_ECx;
  float pip_ECy;
  float pip_ECz;
  float pip_ECin, pip_ECout;
  
  chain->SetBranchAddress("pip_ECx",&pip_ECx);
  chain->SetBranchAddress("pip_ECy",&pip_ECy);
  chain->SetBranchAddress("pip_ECz",&pip_ECz);
  chain->SetBranchAddress("pip_ECin",&pip_ECin);
  chain->SetBranchAddress("pip_ECout",&pip_ECout);
  
  float prot_vz;
  float prot_vy;
  float prot_vx;
  
  chain->SetBranchAddress("prot_vz",&prot_vz);
  chain->SetBranchAddress("prot_vy",&prot_vy);
  chain->SetBranchAddress("prot_vx",&prot_vx);
  
  float pip_vz;
  float pip_vy;
  float pip_vx;
  
  chain->SetBranchAddress("pip_vz",&pip_vz);
  chain->SetBranchAddress("pip_vy",&pip_vy);
  chain->SetBranchAddress("pip_vx",&pip_vx);
  
  float pim_vz;
  float pim_vy;
  float pim_vx;
  
  chain->SetBranchAddress("pim_vz",&pim_vz);
  chain->SetBranchAddress("pim_vy",&pim_vy);
  chain->SetBranchAddress("pim_vx",&pim_vx);
  
  float IVEpEm, MM2PEpEm, MM2P, MEPEpEm;
  
  
  
  chain->SetBranchAddress("IVEpEm",&IVEpEm);
  chain->SetBranchAddress("MM2PEpEm",&MM2PEpEm);
  chain->SetBranchAddress("MM2P",&MM2P);
  chain->SetBranchAddress("MEPEpEm",&MEPEpEm);
  
  Float_t vx, vy, vz;
  chain->SetBranchAddress("vx",&vx);
  chain->SetBranchAddress("vy",&vy);
  chain->SetBranchAddress("vz",&vz);
  
  
  Int_t NProt, NPip, NPim;
  chain->SetBranchAddress("NProt",&NProt);
  chain->SetBranchAddress("NPip",&NPip);
  chain->SetBranchAddress("NPim",&NPim);
  
  
  float pip_dTOF, Timing_pip, pip_Tprop;
  chain->SetBranchAddress("pip_dTOF",&pip_dTOF);
  chain->SetBranchAddress("Timing_pip",&Timing_pip);
  chain->SetBranchAddress("pip_Tprop",&pip_Tprop);
  float pip_scLen, pip_SC_time, pip_ST_time, pip_STSCpl;
  chain->SetBranchAddress("pip_scLen",&pip_scLen);
  chain->SetBranchAddress("pip_SC_time",&pip_SC_time);
  chain->SetBranchAddress("pip_ST_time",&pip_ST_time);
  chain->SetBranchAddress("pip_STSCpl",&pip_STSCpl);
  Int_t pip_SC_paddle;
  chain->SetBranchAddress("pip_SC_paddle",&pip_SC_paddle);
  
  float pim_dTOF, Timing_pim, pim_Tprop;
  chain->SetBranchAddress("pim_dTOF",&pim_dTOF);
  chain->SetBranchAddress("Timing_pim",&Timing_pim);
  chain->SetBranchAddress("pim_Tprop",&pim_Tprop);
  float pim_scLen, pim_SC_time, pim_ST_time, pim_STSCpl;
  chain->SetBranchAddress("pim_scLen",&pim_scLen);
  chain->SetBranchAddress("pim_SC_time",&pim_SC_time);
  chain->SetBranchAddress("pim_ST_time",&pim_ST_time);
  chain->SetBranchAddress("pim_STSCpl",&pim_STSCpl);
  Int_t pim_SC_paddle;
  chain->SetBranchAddress("pim_SC_paddle",&pim_SC_paddle);
  
  float prot_dTOF, Timing_prot, prot_Tprop;
  chain->SetBranchAddress("prot_dTOF",&prot_dTOF);
  chain->SetBranchAddress("Timing_prot",&Timing_prot);
  chain->SetBranchAddress("prot_Tprop",&prot_Tprop);
  float prot_scLen, prot_SC_time, prot_ST_time, prot_STSCpl;
  chain->SetBranchAddress("prot_scLen",&prot_scLen);
  chain->SetBranchAddress("prot_SC_time",&prot_SC_time);
  chain->SetBranchAddress("prot_ST_time",&prot_ST_time);
  chain->SetBranchAddress("prot_STSCpl",&prot_STSCpl);
  Int_t prot_SC_paddle;
  chain->SetBranchAddress("prot_SC_paddle",&prot_SC_paddle);
  
  TLorentzVector P4pho,P4pro,P4em,P4ep,P4target;
  TLorentzVector P4pho_gamfit_diLep,P4pro_gamfit_diLep,P4em_gamfit_diLep,P4ep_gamfit_diLep;
  TLorentzVector P4pho_pi0fit_diLep, P4pro_pi0fit_diLep, P4em_pi0fit_diLep, P4ep_pi0fit_diLep;
  TLorentzVector P4pho_nothingfit_diLep, P4pro_nothingfit_diLep, P4em_nothingfit_diLep, P4ep_nothingfit_diLep;
  
  P4target.SetPxPyPzE(0.0,0.0,0.0,0.93828);
  
  TVector3 V3pro,V3ep,V3em;
  
  //Pion fitted 4vectors
  TLorentzVector P4Pim,P4Pip;
  TLorentzVector P4pho_gamfit_diPion,P4pro_gamfit_diPion,P4em_gamfit_diPion,P4ep_gamfit_diPion;
  TLorentzVector P4pho_pi0fit_diPion, P4pro_pi0fit_diPion, P4em_pi0fit_diPion, P4ep_pi0fit_diPion;
  TLorentzVector P4pho_nothingfit_diPion, P4pro_nothingfit_diPion, P4em_nothingfit_diPion, P4ep_nothingfit_diPion;
  
  
  Int_t nentries = (Int_t)chain->GetEntries();
  
  tot_events = tot_events + nentries;
  
  cout<<"############################################"<<endl;
  cout<<tot_events<<" TOTAL NUMBER PROCESSED"<<endl;
  cout<<"############################################"<<endl;
  
  
  TFile outFile("Lep_KinFit.root","recreate");
  TTree *t4 = new TTree("LepTree","LepTree");
  
  Double_t mE_PEpEm, mm2_PEpEm, mm2_P, E_g, IV_EpEm, IV_EpEm_P, mandelstam_t, CM_Theta, open_angle;
  
  t4->Branch("mE_PEpEm",&mE_PEpEm,"mE_PEpEm/D");
  t4->Branch("mm2_PEpEm",&mm2_PEpEm,"mm2_PEpEm/D");
  t4->Branch("mm2_P",&mm2_P,"mm2_P/D");
  t4->Branch("IV_EpEm",&IV_EpEm,"IV_EpEm/D");
  t4->Branch("IV_EpEm_P",&IV_EpEm_P,"IV_EpEm_P/D");
  t4->Branch("mandelstam_t",&mandelstam_t,"mandelstam_t/D");
  t4->Branch("CM_Theta",&CM_Theta,"CM_Theta/D");
  t4->Branch("open_angle",&open_angle,"open_angle/D");
  
  Double_t P_px, P_py, P_pz, P_p, P_Theta, P_Phi;
  t4->Branch("P_px",&P_px,"P_px/D");
  t4->Branch("P_py",&P_py,"P_py/D");
  t4->Branch("P_pz",&P_pz,"P_pz/D");
  t4->Branch("P_p",&P_p,"P_p/D");
  t4->Branch("P_Theta",&P_Theta,"P_Theta/D");
  t4->Branch("P_Phi",&P_Phi,"P_Phi/D");
  
  Double_t Ep_px, Ep_py, Ep_pz, Ep_p, Ep_Theta, Ep_Phi;
  t4->Branch("Ep_px",&Ep_px,"Ep_px/D");
  t4->Branch("Ep_py",&Ep_py,"Ep_py/D");
  t4->Branch("Ep_pz",&Ep_pz,"Ep_pz/D");
  t4->Branch("Ep_p",&Ep_p,"Ep_p/D");
  t4->Branch("Ep_Theta",&Ep_Theta,"Ep_Theta/D");
  t4->Branch("Ep_Phi",&Ep_Phi,"Ep_Phi/D");
  
  Double_t Em_px, Em_py, Em_pz, Em_p, Em_Theta, Em_Phi;
  t4->Branch("Em_px",&Em_px,"Em_px/D");
  t4->Branch("Em_py",&Em_py,"Em_py/D");
  t4->Branch("Em_pz",&Em_pz,"Em_pz/D");
  t4->Branch("Em_p",&Em_p,"Em_p/D");
  t4->Branch("Em_Theta",&Em_Theta,"Em_Theta/D");
  t4->Branch("Em_Phi",&Em_Phi,"Em_Phi/D");
  
  Int_t P_sec, Ep_sec, Em_sec;
  t4->Branch("P_sec",&P_sec,"P_sec/I");
  t4->Branch("Ep_sec",&Ep_sec,"Ep_sec/I");
  t4->Branch("Em_sec",&Em_sec,"Em_sec/I");
  
  
  Int_t nPim, nPip, nProt;
  t4->Branch("nPim",&nPim,"nPim/I");
  t4->Branch("nPip",&nPip,"nPip/I");
  t4->Branch("nProt",&nProt,"nProt/I");
  
  t4->Branch("E_g",&E_g,"E_g/D");
  
  Double_t P_dTOF, Timing_P, P_Tprop;
  t4->Branch("P_dTOF",&P_dTOF,"P_dTOF/D");
  t4->Branch("Timing_P",&Timing_P,"Timing_P/D");
  t4->Branch("P_Tprop",&P_Tprop,"P_Tprop/D");
  
  Double_t Ep_dTOF, Timing_Ep, Ep_Tprop;
  t4->Branch("Ep_dTOF",&Ep_dTOF,"Ep_dTOF/D");
  t4->Branch("Timing_Ep",&Timing_Ep,"Timing_Ep/D");
  t4->Branch("Ep_Tprop",&Ep_Tprop,"Ep_Tprop/D");
  
  Double_t Em_dTOF, Timing_Em, Em_Tprop;
  t4->Branch("Em_dTOF",&Em_dTOF,"Em_dTOF/D");
  t4->Branch("Timing_Em",&Timing_Em,"Timing_Em/D");
  t4->Branch("Em_Tprop",&Em_Tprop,"Em_Tprop/D");
  
  Double_t Vx, Vy, Vz;
  t4->Branch("Vx",&Vx,"Vx/D");
  t4->Branch("Vy",&Vy,"Vy/D");
  t4->Branch("Vz",&Vz,"Vz/D");
  
  Int_t trip_pass, Em_tofpass, Ep_tofpass, P_tofpass, Ep_EC_pass, Em_EC_pass, Em_geofid, Ep_geofid, P_geofid, Pass_all;
  t4->Branch("trip_pass",&trip_pass,"trip_pass/I");
  t4->Branch("Em_tofpass",&Em_tofpass,"Em_tofpass/I");
  t4->Branch("Ep_tofpass",&Ep_tofpass,"Ep_tofpass/I");
  t4->Branch("P_tofpass",&P_tofpass,"P_tofpass/I");
  t4->Branch("Ep_EC_pass",&Ep_EC_pass,"Ep_EC_pass/I");
  t4->Branch("Em_EC_pass",&Em_EC_pass,"Em_EC_pass/I");
  t4->Branch("Em_geofid",&Em_geofid,"Em_geofid/I");
  t4->Branch("Ep_geofid",&Ep_geofid,"Ep_geofid/I");
  t4->Branch("P_geofid",&P_geofid,"P_geofid/I");
  t4->Branch("Pass_all",&Pass_all,"Pass_all/I");
  
  Int_t Run, Event, TrigBits;
  t4->Branch("Run",&Run,"Run/I");
  t4->Branch("Event",&Event,"Event/I");
  t4->Branch("TrigBits",&TrigBits,"TrigBits/I");
  /*############################################################################*
   *############################################################################*
   *################## Lets Start with e+e- kinematic fitting ##################*
   *##################                3 fits                  ##################*
   *##################              e+e-(gamma)  1-C          ##################*
   *##################              e+e-(pi0)    1-C          ##################*
   *##################              e+e-(0)      4-C          ##################*
   *############################################################################*
   *############################################################################*/
  //Example of test fit
  Double_t test_Pull_Chi, test_Pull_Prob;
  
  t4->Branch("test_Pull_Chi",&test_Pull_Chi,"test_Pull_Chi/D");
  t4->Branch("test_Pull_Prob",&test_Pull_Prob,"test_Pull_Prob/D");
  
  //e+e-(gamma) fit
  Double_t Pull_Chi_gamfit_diLep, Pull_Prob_gamfit_diLep;
  
  t4->Branch("Pull_Chi_gamfit_diLep",&Pull_Chi_gamfit_diLep,"Pull_Chi_gamfit_diLep/D");
  t4->Branch("Pull_Prob_gamfit_diLep",&Pull_Prob_gamfit_diLep,"Pull_Prob_gamfit_diLep/D");
  
  Double_t mE_PEpEm_gamfit_diLep, mm2_PEpEm_gamfit_diLep, mm_PEpEm_gamfit_diLep, mm_P_gamfit_diLep, mm2_P_gamfit_diLep, Eg_gamfit_diLep, IV_EpEm_gamfit_diLep, IV_EpEm_gamfit_diLep_P;
  
  t4->Branch("mE_PEpEm_gamfit_diLep",&mE_PEpEm_gamfit_diLep,"mE_PEpEm_gamfit_diLep/D");
  t4->Branch("mm2_PEpEm_gamfit_diLep",&mm2_PEpEm_gamfit_diLep,"mm2_PEpEm_gamfit_diLep/D");
  t4->Branch("mm_PEpEm_gamfit_diLep",&mm_PEpEm_gamfit_diLep,"mm_PEpEm_gamfit_diLep/D");
  t4->Branch("mm_P_gamfit_diLep",&mm_P_gamfit_diLep,"mm_P_gamfit_diLep/D");
  t4->Branch("mm2_P_gamfit_diLep",&mm2_P_gamfit_diLep,"mm2_P_gamfit_diLep/D");
  t4->Branch("IV_EpEm_gamfit_diLep",&IV_EpEm_gamfit_diLep,"IV_EpEm_gamfit_diLep/D");
  t4->Branch("IV_EpEm_gamfit_diLep_P",&IV_EpEm_gamfit_diLep_P,"IV_EpEm_gamfit_diLep_P/D");
  
  t4->Branch("Eg_gamfit_diLep",&Eg_gamfit_diLep,"Eg_gamfit_diLep/D");
  
  Double_t open_angle_gamfit_diLep;
  t4->Branch("open_angle_gamfit_diLep",&open_angle_gamfit_diLep,"open_angle_gamfit_diLep/D");
  
  Double_t P_px_gamfit_diLep, P_py_gamfit_diLep, P_pz_gamfit_diLep, P_Theta_gamfit_diLep, P_Phi_gamfit_diLep, P_Ptot_gamfit_diLep;
  Double_t Ep_px_gamfit_diLep, Ep_py_gamfit_diLep, Ep_pz_gamfit_diLep, Ep_Theta_gamfit_diLep, Ep_Phi_gamfit_diLep, Ep_Ptot_gamfit_diLep;
  Double_t Em_px_gamfit_diLep, Em_py_gamfit_diLep, Em_pz_gamfit_diLep, Em_Theta_gamfit_diLep, Em_Phi_gamfit_diLep, Em_Ptot_gamfit_diLep;
  
  t4->Branch("P_px_gamfit_diLep",&P_px_gamfit_diLep,"P_px_gamfit_diLep/D");
  t4->Branch("P_py_gamfit_diLep",&P_py_gamfit_diLep,"P_py_gamfit_diLep/D");
  t4->Branch("P_pz_gamfit_diLep",&P_pz_gamfit_diLep,"P_pz_gamfit_diLep/D");
  t4->Branch("P_Theta_gamfit_diLep",&P_Theta_gamfit_diLep,"P_Theta_gamfit_diLep/D");
  t4->Branch("P_Phi_gamfit_diLep",&P_Phi_gamfit_diLep,"P_Phi_gamfit_diLep/D");
  t4->Branch("P_Ptot_gamfit_diLep",&P_Ptot_gamfit_diLep,"P_Ptot_gamfit_diLep/D");
  
  t4->Branch("Ep_px_gamfit_diLep",&Ep_px_gamfit_diLep,"Ep_px_gamfit_diLep/D");
  t4->Branch("Ep_py_gamfit_diLep",&Ep_py_gamfit_diLep,"Ep_py_gamfit_diLep/D");
  t4->Branch("Ep_pz_gamfit_diLep",&Ep_pz_gamfit_diLep,"Ep_pz_gamfit_diLep/D");
  t4->Branch("Ep_Theta_gamfit_diLep",&Ep_Theta_gamfit_diLep,"Ep_Theta_gamfit_diLep/D");
  t4->Branch("Ep_Phi_gamfit_diLep",&Ep_Phi_gamfit_diLep,"Ep_Phi_gamfit_diLep/D");
  t4->Branch("Ep_Ptot_gamfit_diLep",&Ep_Ptot_gamfit_diLep,"Ep_Ptot_gamfit_diLep/D");
  
  t4->Branch("Em_px_gamfit_diLep",&Em_px_gamfit_diLep,"Em_px_gamfit_diLep/D");
  t4->Branch("Em_py_gamfit_diLep",&Em_py_gamfit_diLep,"Em_py_gamfit_diLep/D");
  t4->Branch("Em_pz_gamfit_diLep",&Em_pz_gamfit_diLep,"Em_pz_gamfit_diLep/D");
  t4->Branch("Em_Theta_gamfit_diLep",&Em_Theta_gamfit_diLep,"Em_Theta_gamfit_diLep/D");
  t4->Branch("Em_Phi_gamfit_diLep",&Em_Phi_gamfit_diLep,"Em_Phi_gamfit_diLep/D");
  t4->Branch("Em_Ptot_gamfit_diLep",&Em_Ptot_gamfit_diLep,"Em_Ptot_gamfit_diLep/D");
  
  Double_t CM_Theta_gamfit_diLep;
  t4->Branch("CM_Theta_gamfit_diLep",&CM_Theta_gamfit_diLep,"CM_Theta_gamfit_diLep/D");
  
  Double_t mandelstam_t_gamfit_diLep;
  t4->Branch("mandelstam_t_gamfit_diLep",&mandelstam_t_gamfit_diLep,"mandelstam_t_gamfit_diLep/D");
  
  
  //e+e-(pi0) fit
  Double_t Pull_Chi_pi0fit_diLep, Pull_Prob_pi0fit_diLep;
  
  t4->Branch("Pull_Chi_pi0fit_diLep",&Pull_Chi_pi0fit_diLep,"Pull_Chi_pi0fit_diLep/D");
  t4->Branch("Pull_Prob_pi0fit_diLep",&Pull_Prob_pi0fit_diLep,"Pull_Prob_pi0fit_diLep/D");
  
  Double_t mE_PEpEm_pi0fit_diLep, mm2_PEpEm_pi0fit_diLep, mm_PEpEm_pi0fit_diLep, mm_P_pi0fit_diLep, mm2_P_pi0fit_diLep, Eg_pi0fit_diLep, IV_EpEm_pi0fit_diLep, IV_EpEm_pi0fit_diLep_P;
  
  t4->Branch("mE_PEpEm_pi0fit_diLep",&mE_PEpEm_pi0fit_diLep,"mE_PEpEm_pi0fit_diLep/D");
  t4->Branch("mm2_PEpEm_pi0fit_diLep",&mm2_PEpEm_pi0fit_diLep,"mm2_PEpEm_pi0fit_diLep/D");
  t4->Branch("mm_PEpEm_pi0fit_diLep",&mm_PEpEm_pi0fit_diLep,"mm_PEpEm_pi0fit_diLep/D");
  t4->Branch("mm_P_pi0fit_diLep",&mm_P_pi0fit_diLep,"mm_P_pi0fit_diLep/D");
  t4->Branch("mm2_P_pi0fit_diLep",&mm2_P_pi0fit_diLep,"mm2_P_pi0fit_diLep/D");
  t4->Branch("IV_EpEm_pi0fit_diLep",&IV_EpEm_pi0fit_diLep,"IV_EpEm_pi0fit_diLep/D");
  t4->Branch("IV_EpEm_pi0fit_diLep_P",&IV_EpEm_pi0fit_diLep_P,"IV_EpEm_pi0fit_diLep_P/D");
  
  t4->Branch("Eg_pi0fit_diLep",&Eg_pi0fit_diLep,"Eg_pi0fit_diLep/D");
  
  Double_t open_angle_pi0fit_diLep;
  t4->Branch("open_angle_pi0fit_diLep",&open_angle_pi0fit_diLep,"open_angle_pi0fit_diLep/D");
  
  Double_t P_px_pi0fit_diLep, P_py_pi0fit_diLep, P_pz_pi0fit_diLep, P_Theta_pi0fit_diLep, P_Phi_pi0fit_diLep, P_Ptot_pi0fit_diLep;
  Double_t Ep_px_pi0fit_diLep, Ep_py_pi0fit_diLep, Ep_pz_pi0fit_diLep, Ep_Theta_pi0fit_diLep, Ep_Phi_pi0fit_diLep, Ep_Ptot_pi0fit_diLep;
  Double_t Em_px_pi0fit_diLep, Em_py_pi0fit_diLep, Em_pz_pi0fit_diLep, Em_Theta_pi0fit_diLep, Em_Phi_pi0fit_diLep, Em_Ptot_pi0fit_diLep;
  
  t4->Branch("P_px_pi0fit_diLep",&P_px_pi0fit_diLep,"P_px_pi0fit_diLep/D");
  t4->Branch("P_py_pi0fit_diLep",&P_py_pi0fit_diLep,"P_py_pi0fit_diLep/D");
  t4->Branch("P_pz_pi0fit_diLep",&P_pz_pi0fit_diLep,"P_pz_pi0fit_diLep/D");
  t4->Branch("P_Theta_pi0fit_diLep",&P_Theta_pi0fit_diLep,"P_Theta_pi0fit_diLep/D");
  t4->Branch("P_Phi_pi0fit_diLep",&P_Phi_pi0fit_diLep,"P_Phi_pi0fit_diLep/D");
  t4->Branch("P_Ptot_pi0fit_diLep",&P_Ptot_pi0fit_diLep,"P_Ptot_pi0fit_diLep/D");
  
  t4->Branch("Ep_px_pi0fit_diLep",&Ep_px_pi0fit_diLep,"Ep_px_pi0fit_diLep/D");
  t4->Branch("Ep_py_pi0fit_diLep",&Ep_py_pi0fit_diLep,"Ep_py_pi0fit_diLep/D");
  t4->Branch("Ep_pz_pi0fit_diLep",&Ep_pz_pi0fit_diLep,"Ep_pz_pi0fit_diLep/D");
  t4->Branch("Ep_Theta_pi0fit_diLep",&Ep_Theta_pi0fit_diLep,"Ep_Theta_pi0fit_diLep/D");
  t4->Branch("Ep_Phi_pi0fit_diLep",&Ep_Phi_pi0fit_diLep,"Ep_Phi_pi0fit_diLep/D");
  t4->Branch("Ep_Ptot_pi0fit_diLep",&Ep_Ptot_pi0fit_diLep,"Ep_Ptot_pi0fit_diLep/D");
  
  t4->Branch("Em_px_pi0fit_diLep",&Em_px_pi0fit_diLep,"Em_px_pi0fit_diLep/D");
  t4->Branch("Em_py_pi0fit_diLep",&Em_py_pi0fit_diLep,"Em_py_pi0fit_diLep/D");
  t4->Branch("Em_pz_pi0fit_diLep",&Em_pz_pi0fit_diLep,"Em_pz_pi0fit_diLep/D");
  t4->Branch("Em_Theta_pi0fit_diLep",&Em_Theta_pi0fit_diLep,"Em_Theta_pi0fit_diLep/D");
  t4->Branch("Em_Phi_pi0fit_diLep",&Em_Phi_pi0fit_diLep,"Em_Phi_pi0fit_diLep/D");
  t4->Branch("Em_Ptot_pi0fit_diLep",&Em_Ptot_pi0fit_diLep,"Em_Ptot_pi0fit_diLep/D");
  
  Double_t CM_Theta_pi0fit_diLep;
  t4->Branch("CM_Theta_pi0fit_diLep",&CM_Theta_pi0fit_diLep,"CM_Theta_pi0fit_diLep/D");
  
  Double_t mandelstam_t_pi0fit_diLep;
  t4->Branch("mandelstam_t_pi0fit_diLep",&mandelstam_t_pi0fit_diLep,"mandelstam_t_pi0fit_diLep/D");
  
  //e+e-(0) fit
  Double_t Pull_Chi_nothingfit_diLep, Pull_Prob_nothingfit_diLep;
  
  t4->Branch("Pull_Chi_nothingfit_diLep",&Pull_Chi_nothingfit_diLep,"Pull_Chi_nothingfit_diLep/D");
  t4->Branch("Pull_Prob_nothingfit_diLep",&Pull_Prob_nothingfit_diLep,"Pull_Prob_nothingfit_diLep/D");
  
  Double_t Pull_Zero, Pull_One, Pull_Two, Pull_Three, Pull_Four, Pull_Five, Pull_Six, Pull_Seven, Pull_Eight, Pull_Nine;
  t4->Branch("Pull_Zero",&Pull_Zero,"Pull_Zero/D");
  t4->Branch("Pull_One",&Pull_One,"Pull_One/D");
  t4->Branch("Pull_Two",&Pull_Two,"Pull_Two/D");
  t4->Branch("Pull_Three",&Pull_Three,"Pull_Three/D");
  t4->Branch("Pull_Four",&Pull_Four,"Pull_Four/D");
  t4->Branch("Pull_Five",&Pull_Five,"Pull_Five/D");
  t4->Branch("Pull_Six",&Pull_Six,"Pull_Six/D");
  t4->Branch("Pull_Seven",&Pull_Seven,"Pull_Seven/D");
  t4->Branch("Pull_Eight",&Pull_Eight,"Pull_Eight/D");
  t4->Branch("Pull_Nine",&Pull_Nine,"Pull_Nine/D");
  
  Double_t mE_PEpEm_nothingfit_diLep, mm2_PEpEm_nothingfit_diLep, mm_PEpEm_nothingfit_diLep, mm_P_nothingfit_diLep, mm2_P_nothingfit_diLep, Eg_nothingfit_diLep, IV_EpEm_nothingfit_diLep, IV_EpEm_nothingfit_diLep_P;
  
  t4->Branch("mE_PEpEm_nothingfit_diLep",&mE_PEpEm_nothingfit_diLep,"mE_PEpEm_nothingfit_diLep/D");
  t4->Branch("mm2_PEpEm_nothingfit_diLep",&mm2_PEpEm_nothingfit_diLep,"mm2_PEpEm_nothingfit_diLep/D");
  t4->Branch("mm_PEpEm_nothingfit_diLep",&mm_PEpEm_nothingfit_diLep,"mm_PEpEm_nothingfit_diLep/D");
  t4->Branch("mm_P_nothingfit_diLep",&mm_P_nothingfit_diLep,"mm_P_nothingfit_diLep/D");
  t4->Branch("mm2_P_nothingfit_diLep",&mm2_P_nothingfit_diLep,"mm2_P_nothingfit_diLep/D");
  t4->Branch("IV_EpEm_nothingfit_diLep",&IV_EpEm_nothingfit_diLep,"IV_EpEm_nothingfit_diLep/D");
  t4->Branch("IV_EpEm_nothingfit_diLep_P",&IV_EpEm_nothingfit_diLep_P,"IV_EpEm_nothingfit_diLep_P/D");
  
  t4->Branch("Eg_nothingfit_diLep",&Eg_nothingfit_diLep,"Eg_nothingfit_diLep/D");
  
  Double_t open_angle_nothingfit_diLep;
  t4->Branch("open_angle_nothingfit_diLep",&open_angle_nothingfit_diLep,"open_angle_nothingfit_diLep/D");
  
  Double_t P_px_nothingfit_diLep, P_py_nothingfit_diLep, P_pz_nothingfit_diLep, P_Theta_nothingfit_diLep, P_Phi_nothingfit_diLep, P_Ptot_nothingfit_diLep;
  Double_t Ep_px_nothingfit_diLep, Ep_py_nothingfit_diLep, Ep_pz_nothingfit_diLep, Ep_Theta_nothingfit_diLep, Ep_Phi_nothingfit_diLep, Ep_Ptot_nothingfit_diLep;
  Double_t Em_px_nothingfit_diLep, Em_py_nothingfit_diLep, Em_pz_nothingfit_diLep, Em_Theta_nothingfit_diLep, Em_Phi_nothingfit_diLep, Em_Ptot_nothingfit_diLep;
  
  t4->Branch("P_px_nothingfit_diLep",&P_px_nothingfit_diLep,"P_px_nothingfit_diLep/D");
  t4->Branch("P_py_nothingfit_diLep",&P_py_nothingfit_diLep,"P_py_nothingfit_diLep/D");
  t4->Branch("P_pz_nothingfit_diLep",&P_pz_nothingfit_diLep,"P_pz_nothingfit_diLep/D");
  t4->Branch("P_Theta_nothingfit_diLep",&P_Theta_nothingfit_diLep,"P_Theta_nothingfit_diLep/D");
  t4->Branch("P_Phi_nothingfit_diLep",&P_Phi_nothingfit_diLep,"P_Phi_nothingfit_diLep/D");
  t4->Branch("P_Ptot_nothingfit_diLep",&P_Ptot_nothingfit_diLep,"P_Ptot_nothingfit_diLep/D");
  
  t4->Branch("Ep_px_nothingfit_diLep",&Ep_px_nothingfit_diLep,"Ep_px_nothingfit_diLep/D");
  t4->Branch("Ep_py_nothingfit_diLep",&Ep_py_nothingfit_diLep,"Ep_py_nothingfit_diLep/D");
  t4->Branch("Ep_pz_nothingfit_diLep",&Ep_pz_nothingfit_diLep,"Ep_pz_nothingfit_diLep/D");
  t4->Branch("Ep_Theta_nothingfit_diLep",&Ep_Theta_nothingfit_diLep,"Ep_Theta_nothingfit_diLep/D");
  t4->Branch("Ep_Phi_nothingfit_diLep",&Ep_Phi_nothingfit_diLep,"Ep_Phi_nothingfit_diLep/D");
  t4->Branch("Ep_Ptot_nothingfit_diLep",&Ep_Ptot_nothingfit_diLep,"Ep_Ptot_nothingfit_diLep/D");
  
  t4->Branch("Em_px_nothingfit_diLep",&Em_px_nothingfit_diLep,"Em_px_nothingfit_diLep/D");
  t4->Branch("Em_py_nothingfit_diLep",&Em_py_nothingfit_diLep,"Em_py_nothingfit_diLep/D");
  t4->Branch("Em_pz_nothingfit_diLep",&Em_pz_nothingfit_diLep,"Em_pz_nothingfit_diLep/D");
  t4->Branch("Em_Theta_nothingfit_diLep",&Em_Theta_nothingfit_diLep,"Em_Theta_nothingfit_diLep/D");
  t4->Branch("Em_Phi_nothingfit_diLep",&Em_Phi_nothingfit_diLep,"Em_Phi_nothingfit_diLep/D");
  t4->Branch("Em_Ptot_nothingfit_diLep",&Em_Ptot_nothingfit_diLep,"Em_Ptot_nothingfit_diLep/D");
  
  Double_t CM_Theta_nothingfit_diLep;
  t4->Branch("CM_Theta_nothingfit_diLep",&CM_Theta_nothingfit_diLep,"CM_Theta_nothingfit_diLep/D");
  
  Double_t mandelstam_t_nothingfit_diLep;
  t4->Branch("mandelstam_t_nothingfit_diLep",&mandelstam_t_nothingfit_diLep,"mandelstam_t_nothingfit_diLep/D");
  
  
  Double_t M_P = 0.938272;   //Proton
  Double_t M_Pi = 0.139570;  //Pion
  Double_t M_PiZ = 0.1349766;  //Pion Zero
  Double_t M_Eta = 0.547853;  //Eta
  Double_t M_Eta_Prime = 0.95778;  //Eta_Prime
  Double_t Melectron = 0.000510999; //Electron
  Double_t c = 29.9792458; //units m/s/e7
  Double_t pi = TMath::Pi();
  Double_t DegToRad = pi/180.0;
  
  
  for (Int_t j=0;j<=nentries;j++) {//nentries
    chain->GetEntry(j);
    
    if(!(j%9500)) std::cout << "\r done " << j << " out of " << nentries << " ==> " << double(j)*100.0/double(nentries) << "%" << flush;
    if(j== nentries) std::cout << " DONE" << endl;
    
    Pass_all = 0;
    
    if (NPip ==1 && NPim ==1 && NProt ==1) {
      //Lets Loop through each in-time photon
      for (Int_t jphoton = 0; jphoton<N_in_timePhoton; jphoton++){
        
        //####################### Beam Corrections #######################
        Double_t egam_chosen = in_photon_E[jphoton];
        
        Int_t run_input;
        if(run == 56400){run_input = 56401;}
        else if(run == 57314){run_input = 57315;}
        else{run_input = run;}
        double egam_chosen_corrected = clas::g12::corrected_beam_energy(run_input, egam_chosen);
        //####################### End Beam Corrections #######################
        
        
        P_p = prot_P + pcor.pcor(prot_Phi*DegToRad,14); //Momemtum corrected
        P_Theta = prot_Theta;
        P_Phi = prot_Phi;
        P_px = (P_p*sin(DegToRad*P_Theta)*cos(DegToRad*P_Phi));
        P_py = (P_p*sin(DegToRad*P_Theta)*sin(DegToRad*P_Phi));
        P_pz = (P_p*cos(DegToRad*P_Theta));
        
        Em_p = pim_P + pcor.pcor(pim_Phi*DegToRad,9); //Momemtum corrected
        Em_Theta = pim_Theta;
        Em_Phi = pim_Phi;
        Em_px = Em_p*sin(DegToRad*Em_Theta)*cos(DegToRad*Em_Phi);
        Em_py = Em_p*sin(DegToRad*Em_Theta)*sin(DegToRad*Em_Phi);
        Em_pz = Em_p*cos(DegToRad*Em_Theta);
        
        
        Ep_p = pip_P + pcor.pcor(pip_Phi*DegToRad,8); //Momemtum corrected
        Ep_Theta = pip_Theta;
        Ep_Phi = pip_Phi;
        Ep_px = Ep_p*sin(DegToRad*Ep_Theta)*cos(DegToRad*Ep_Phi);
        Ep_py = Ep_p*sin(DegToRad*Ep_Theta)*sin(DegToRad*Ep_Phi);
        Ep_pz = Ep_p*cos(DegToRad*Ep_Theta);
        
        P4pho.SetPxPyPzE(0.0,0.0,egam_chosen_corrected,egam_chosen_corrected);
        P4em.SetPxPyPzE(Em_px,Em_py,Em_pz,sqrt(Em_p*Em_p + Melectron*Melectron));
        P4pro.SetPxPyPzE(P_px,P_py,P_pz,sqrt(P_p*P_p + M_P*M_P));
        P4ep.SetPxPyPzE(Ep_px,Ep_py,Ep_pz,sqrt(Ep_p*Ep_p + Melectron*Melectron));
        V3pro.SetXYZ(prot_vx,prot_vy,prot_vz);
        V3ep.SetXYZ(pip_vx,pip_vy,pip_vz);
        V3em.SetXYZ(pim_vx,pim_vy,pim_vz);
        
        TLorentzVector P4tot, P4mis;
        
        P4tot = P4pho + P4target;
        P4mis = P4tot - ( P4pro + P4em + P4ep );
        
        if(abs((P4tot - ( P4em + P4ep )).M()-0.93828)< 0.055 && abs(P4mis.M())< 0.035){
          
          //Lets start with fiducial and TOF cuts
          
          trip_pass = clas::g12::is_good(run, event);
          
          TVector3 pimUVW = clas::g12::g12_ECxyz_2uvw(pim_ECx, pim_ECy, pim_ECz);// EC UVW
          TVector3 pipUVW = clas::g12::g12_ECxyz_2uvw(pip_ECx, pip_ECy, pip_ECz);// EC UVW
          float pim_u = pimUVW.X();
          float pim_v = pimUVW.Y();
          float pim_w = pimUVW.Z();
          float pip_u = pipUVW.X();
          float pip_v = pipUVW.Y();
          float pip_w = pipUVW.Z();
          
          Ep_EC_pass = clas::g12::pass_g12_ec_knockout(pip_ECin, pip_ECout, pip_u, pip_v, pip_w, pip_sec);// EC Knockout
          Em_EC_pass = clas::g12::pass_g12_ec_knockout(pim_ECin, pim_ECout, pim_u, pim_v, pim_w, pim_sec);// EC Knockout
          
          Em_tofpass = clas::g12::pass_g12_TOFKO(pim_sec, pim_SC_paddle, 1);// TOF Knockout
          Ep_tofpass = clas::g12::pass_g12_TOFKO(pip_sec, pip_SC_paddle, 1);// TOF Knockout
          P_tofpass = clas::g12::pass_g12_TOFKO(prot_sec, prot_SC_paddle, 1);// TOF Knockout
          
          Ep_geofid = clas::g12::g12_PosParticle_fiducial_cuts(Em_p, Em_Theta, Em_Phi,"nominal");//Geometric Fiducial Cut
          Em_geofid = clas::g12::g12_NegParticle_fiducial_cuts(Ep_p, Ep_Theta, Ep_Phi,"nominal");//Geometric Fiducial Cut
          P_geofid = clas::g12::g12_PosParticle_fiducial_cuts(P_p, P_Theta, P_Phi,"nominal");//Geometric Fiducial Cut
          
          
          if (Em_tofpass && Ep_tofpass && P_tofpass && Ep_EC_pass && Em_EC_pass && Ep_geofid && Em_geofid && P_geofid)
          {
            Good_events++;
            Pass_all = 1;
            
          }
          else{
            Bad_events++;
            Pass_all = -1;
          }
          //End Fiducial and TOF cuts
          
          Run = run;
          Event = event;
          TrigBits = trigBits;
          
          P_sec = prot_sec;
          Ep_sec = pip_sec;
          Em_sec = pip_sec;
          
          
          P_dTOF = prot_dTOF;
          Timing_P = Timing_prot;
          P_Tprop = prot_Tprop;
          
          
          Ep_dTOF = pip_dTOF;
          Timing_Ep = Timing_pip;
          Ep_Tprop = pip_Tprop;
          
          
          Em_dTOF = pim_dTOF;
          Timing_Em = Timing_pim;
          Em_Tprop = pim_Tprop;
          
          
          mE_PEpEm = (double)MEPEpEm; mm2_PEpEm = (double)MM2PEpEm;
          mm2_P = (double)MM2P; E_g = (double)egam_chosen_corrected; IV_EpEm = (double)IVEpEm;
          Vx = (double)vx; Vy = (double)vy; Vz = (double)vz;
          
          
          open_angle = P4ep.Vect().Angle(P4em.Vect());
          TLorentzVector IV_EpEmVec_nofit = P4em + P4ep;
          IV_EpEm_P= sqrt(pow(IV_EpEmVec_nofit.Px(),2) + pow(IV_EpEmVec_nofit.Py(),2) + pow(IV_EpEmVec_nofit.Pz(),2));
          
          CM_Theta = OneBoostComparison((P4pho + P4target - P4pro), (P4pho + P4target));
          mandelstam_t = (P4target - P4pro).M2();
          
          
          TMatrixD covTrack(10,10);
          Double_t c00 = pow((0.001*5.715),2)/3;
          covTrack(0,0)=c00;
          covTrack(1,1)= prot_mat[0][0]*pow(P_p,4);
          covTrack(1,2)=-prot_mat[0][1]*pow(P_p,2);
          covTrack(1,3)= -prot_mat[0][2]*pow(P_p,2);
          covTrack(2,1)= -prot_mat[0][1]*pow(P_p,2);
          covTrack(2,2)= prot_mat[1][1];
          covTrack(2,3)=prot_mat[1][2];
          covTrack(3,1)= -prot_mat[0][2]*pow(P_p,2);
          covTrack(3,2)= prot_mat[1][2];
          covTrack(3,3)=prot_mat[2][2];
          
          covTrack(4,4)= pip_mat[0][0]*pow(Ep_p,4);
          covTrack(4,5)= -pip_mat[0][1]*pow(Ep_p,2);
          covTrack(4,6)= -pip_mat[0][2]*pow(Ep_p,2);
          covTrack(5,4)= -pip_mat[0][1]*pow(Ep_p,2);
          covTrack(5,5)= pip_mat[1][1];
          covTrack(5,6)= pip_mat[1][2];
          covTrack(6,4)= -pip_mat[0][2]*pow(Ep_p,2);
          covTrack(6,5)= pip_mat[1][2];
          covTrack(6,6)= pip_mat[2][2];
          
          covTrack(7,7)= pim_mat[0][0]*pow(Em_p,4);
          covTrack(7,8)= pim_mat[0][1]*pow(Em_p,2);
          covTrack(7,9)= pim_mat[0][2]*pow(Em_p,2);
          covTrack(8,7)= pim_mat[0][1]*pow(Em_p,2);
          covTrack(8,8)= pim_mat[1][1];
          covTrack(8,9)= pim_mat[1][2];
          covTrack(9,7)= pim_mat[0][2]*pow(Em_p,2);
          covTrack(9,8)= pim_mat[1][2];
          covTrack(9,9)= pim_mat[2][2];
          
          const int num_parts = 3;
          
          std::vector<TLorentzVector> p4(num_parts);
          std::vector<TLorentzVector> ptest(num_parts);
          std::vector<TLorentzVector> p4pi0(num_parts);
          std::vector<TLorentzVector> p4nothing(num_parts);
          
          
          std::vector<TVector3> vert(num_parts);
          std::vector<string> particles(num_parts);
          
          bool multi = false;
          bool is_mc = false;
          
          std::vector<bool> set(num_parts);
          std::vector<bool> settest(num_parts);
          
          set[0] = false;     settest[0] = false;
          set[1] = false;     settest[1] = true;
          set[2] = false;     settest[2] = true;
          
          Double_t m_targ = 0.93828;
          Double_t e_gamma = egam_chosen_corrected;
          
          p4[0] = P4pro; ptest[0] = P4pro;  p4pi0[0] = P4pro; p4nothing[0] = P4pro;
          p4[1] = P4ep;  ptest[1] = P4ep;   p4pi0[1] = P4ep;  p4nothing[1] = P4ep;
          p4[2] = P4em;  ptest[2] = P4em;   p4pi0[2] = P4em;  p4nothing[2] = P4em;
          
          
          vert[0] = V3pro;
          vert[1] = V3ep;
          vert[2] = V3em;
          
          particles[0] = "p";
          particles[1] = "e+";
          particles[2] = "e-";
          
          string experiment = "g12";
          
          //This is an example of a 2-C fit
          TMatrixD covtestMatrix(10,10);
          
          covtestMatrix = CorrectCLAS_V(covTrack,particles,ptest,vert,multi,is_mc,experiment);
          Kstream testfit;
          testfit.StringNames(particles);
          testfit.FitInput(e_gamma,ptest,covtestMatrix,m_targ);
          bool include = true; //Include missing particle in constraint no (false).
          
          testfit.Fit("gamma",settest,include,M_Eta);
          
          
          test_Pull_Chi = testfit.Chi2();
          test_Pull_Prob = testfit.Prob();
          
          //For e+e-(gamma) fitting
          TMatrixD covMatrix(10,10);
          
          covMatrix = CorrectCLAS_V(covTrack,particles,p4,vert,multi,is_mc,experiment);
          Kstream gamfit_diLep;
          gamfit_diLep.StringNames(particles);
          gamfit_diLep.FitInput(e_gamma,p4,covMatrix,m_targ);
          gamfit_diLep.Fit("gamma");
          
          Pull_Chi_gamfit_diLep = gamfit_diLep.Chi2();
          Pull_Prob_gamfit_diLep = gamfit_diLep.Prob();
          
          Eg_gamfit_diLep = gamfit_diLep.FitPhotonEnergy();
          
          for(int i = 0; i < 3; i++) p4[i] = gamfit_diLep.FitP4(i);
          
          P4pho_gamfit_diLep.SetPxPyPzE(0.0,0.0,Eg_gamfit_diLep,Eg_gamfit_diLep);
          P4pro_gamfit_diLep = p4[0];
          P4ep_gamfit_diLep = p4[1];
          P4em_gamfit_diLep = p4[2];
          
          TLorentzVector MM_PVec = (P4pho_gamfit_diLep + P4target) - P4pro_gamfit_diLep;
          TLorentzVector MM_PEpEmVec = (P4pho_gamfit_diLep + P4target) - (P4pro_gamfit_diLep + P4ep_gamfit_diLep + P4em_gamfit_diLep);
          TLorentzVector IV_EpEmVec = P4ep_gamfit_diLep + P4em_gamfit_diLep;
          open_angle_gamfit_diLep = P4ep_gamfit_diLep.Vect().Angle(P4em_gamfit_diLep.Vect());
          IV_EpEm_gamfit_diLep_P = sqrt(pow(IV_EpEmVec.Px(),2) + pow(IV_EpEmVec.Py(),2) + pow(IV_EpEmVec.Pz(),2));
          
          CM_Theta_gamfit_diLep = OneBoostComparison(MM_PVec, (P4pho_gamfit_diLep + P4target));
          mandelstam_t_gamfit_diLep = (P4target - P4pro_gamfit_diLep).M2();
          
          
          mE_PEpEm_gamfit_diLep = MM_PEpEmVec.E();
          mm2_PEpEm_gamfit_diLep = MM_PEpEmVec.M2();
          mm_PEpEm_gamfit_diLep = MM_PEpEmVec.M();
          
          mm_P_gamfit_diLep = MM_PVec.M();
          mm2_P_gamfit_diLep = MM_PVec.M2();
          IV_EpEm_gamfit_diLep = IV_EpEmVec.M();
          
          P_px_gamfit_diLep = P4pro_gamfit_diLep.Px();
          P_py_gamfit_diLep = P4pro_gamfit_diLep.Py();
          P_pz_gamfit_diLep = P4pro_gamfit_diLep.Pz();
          P_Theta_gamfit_diLep = P4pro_gamfit_diLep.Theta()*180./TMath::Pi();
          P_Phi_gamfit_diLep = P4pro_gamfit_diLep.Phi()*180./TMath::Pi();
          
          Ep_px_gamfit_diLep = P4ep_gamfit_diLep.Px();
          Ep_py_gamfit_diLep = P4ep_gamfit_diLep.Py();
          Ep_pz_gamfit_diLep = P4ep_gamfit_diLep.Pz();
          Ep_Theta_gamfit_diLep = P4ep_gamfit_diLep.Theta()*180./TMath::Pi();
          Ep_Phi_gamfit_diLep = P4ep_gamfit_diLep.Phi()*180./TMath::Pi();
          
          Em_px_gamfit_diLep = P4em_gamfit_diLep.Px();
          Em_py_gamfit_diLep = P4em_gamfit_diLep.Py();
          Em_pz_gamfit_diLep = P4em_gamfit_diLep.Pz();
          Em_Theta_gamfit_diLep = P4em_gamfit_diLep.Theta()*180./TMath::Pi();
          Em_Phi_gamfit_diLep = P4em_gamfit_diLep.Phi()*180./TMath::Pi();
          
          P_Ptot_gamfit_diLep = sqrt(pow(P_px_gamfit_diLep,2) + pow(P_py_gamfit_diLep,2) + pow(P_pz_gamfit_diLep,2));
          Ep_Ptot_gamfit_diLep = sqrt(pow(Ep_px_gamfit_diLep,2) + pow(Ep_py_gamfit_diLep,2) + pow(Ep_pz_gamfit_diLep,2));
          Em_Ptot_gamfit_diLep = sqrt(pow(Em_px_gamfit_diLep,2) + pow(Em_py_gamfit_diLep,2) + pow(Em_pz_gamfit_diLep,2));
          
          //For e+e-(pi0) fitting
          TMatrixD covMatrixpi0(10,10);
          
          covMatrixpi0 = CorrectCLAS_V(covTrack,particles,p4pi0,vert,multi,is_mc,experiment);
          Kstream pi0fit_diLep;
          pi0fit_diLep.StringNames(particles);
          pi0fit_diLep.FitInput(e_gamma,p4pi0,covMatrixpi0,m_targ);
          pi0fit_diLep.Fit("pi0");
          
          
          Pull_Chi_pi0fit_diLep = pi0fit_diLep.Chi2();
          Pull_Prob_pi0fit_diLep = pi0fit_diLep.Prob();
          
          Eg_pi0fit_diLep = pi0fit_diLep.FitPhotonEnergy();
          
          for(int i = 0; i < 3; i++) p4pi0[i] = pi0fit_diLep.FitP4(i);
          
          P4pho_pi0fit_diLep.SetPxPyPzE(0.0,0.0,Eg_pi0fit_diLep,Eg_pi0fit_diLep);
          P4pro_pi0fit_diLep = p4pi0[0];
          P4ep_pi0fit_diLep = p4pi0[1];
          P4em_pi0fit_diLep = p4pi0[2];
          
          TLorentzVector MM_PVec_pi0fit_diLep = (P4pho_pi0fit_diLep + P4target) - P4pro_pi0fit_diLep;
          TLorentzVector MM_PEpEmVec_pi0fit_diLep = (P4pho_pi0fit_diLep + P4target) - (P4pro_pi0fit_diLep + P4ep_pi0fit_diLep + P4em_pi0fit_diLep);
          TLorentzVector IV_EpEmVec_pi0fit_diLep = P4ep_pi0fit_diLep + P4em_pi0fit_diLep;
          open_angle_pi0fit_diLep = P4ep_pi0fit_diLep.Vect().Angle(P4em_pi0fit_diLep.Vect());
          IV_EpEm_pi0fit_diLep_P = sqrt(pow(IV_EpEmVec_pi0fit_diLep.Px(),2) + pow(IV_EpEmVec_pi0fit_diLep.Py(),2) + pow(IV_EpEmVec_pi0fit_diLep.Pz(),2));
          
          CM_Theta_pi0fit_diLep = OneBoostComparison(MM_PVec_pi0fit_diLep, (P4pho_pi0fit_diLep + P4target));
          mandelstam_t_pi0fit_diLep = (P4target - P4pro_pi0fit_diLep).M2();
          
          
          mE_PEpEm_pi0fit_diLep = MM_PEpEmVec_pi0fit_diLep.E();
          mm2_PEpEm_pi0fit_diLep = MM_PEpEmVec_pi0fit_diLep.M2();
          mm_PEpEm_pi0fit_diLep = MM_PEpEmVec_pi0fit_diLep.M();
          
          mm_P_pi0fit_diLep = MM_PVec_pi0fit_diLep.M();
          mm2_P_pi0fit_diLep = MM_PVec_pi0fit_diLep.M2();
          IV_EpEm_pi0fit_diLep = IV_EpEmVec_pi0fit_diLep.M();
          
          P_px_pi0fit_diLep = P4pro_pi0fit_diLep.Px();
          P_py_pi0fit_diLep = P4pro_pi0fit_diLep.Py();
          P_pz_pi0fit_diLep = P4pro_pi0fit_diLep.Pz();
          P_Theta_pi0fit_diLep = P4pro_pi0fit_diLep.Theta()*180./TMath::Pi();
          P_Phi_pi0fit_diLep = P4pro_pi0fit_diLep.Phi()*180./TMath::Pi();
          
          Ep_px_pi0fit_diLep = P4ep_pi0fit_diLep.Px();
          Ep_py_pi0fit_diLep = P4ep_pi0fit_diLep.Py();
          Ep_pz_pi0fit_diLep = P4ep_pi0fit_diLep.Pz();
          Ep_Theta_pi0fit_diLep = P4ep_pi0fit_diLep.Theta()*180./TMath::Pi();
          Ep_Phi_pi0fit_diLep = P4ep_pi0fit_diLep.Phi()*180./TMath::Pi();
          
          Em_px_pi0fit_diLep = P4em_pi0fit_diLep.Px();
          Em_py_pi0fit_diLep = P4em_pi0fit_diLep.Py();
          Em_pz_pi0fit_diLep = P4em_pi0fit_diLep.Pz();
          Em_Theta_pi0fit_diLep = P4em_pi0fit_diLep.Theta()*180./TMath::Pi();
          Em_Phi_pi0fit_diLep = P4em_pi0fit_diLep.Phi()*180./TMath::Pi();
          
          P_Ptot_pi0fit_diLep = sqrt(pow(P_px_pi0fit_diLep,2) + pow(P_py_pi0fit_diLep,2) + pow(P_pz_pi0fit_diLep,2));
          Ep_Ptot_pi0fit_diLep = sqrt(pow(Ep_px_pi0fit_diLep,2) + pow(Ep_py_pi0fit_diLep,2) + pow(Ep_pz_pi0fit_diLep,2));
          Em_Ptot_pi0fit_diLep = sqrt(pow(Em_px_pi0fit_diLep,2) + pow(Em_py_pi0fit_diLep,2) + pow(Em_pz_pi0fit_diLep,2));
          
          //For e+e-(0) fitting
          TMatrixD covMatrixnothing(10,10);
          
          covMatrixnothing = CorrectCLAS_V(covTrack,particles,p4nothing,vert,multi,is_mc,experiment);
          Kstream nothingfit_diLep;
          nothingfit_diLep.StringNames(particles);
          nothingfit_diLep.FitInput(e_gamma,p4nothing,covMatrixnothing,m_targ);
          nothingfit_diLep.Fit();
          
          
          Pull_Chi_nothingfit_diLep = nothingfit_diLep.Chi2();
          Pull_Prob_nothingfit_diLep = nothingfit_diLep.Prob();
          
          Pull_Zero = nothingfit_diLep.GetPull(0);
          Pull_One = nothingfit_diLep.GetPull(1);
          Pull_Two = nothingfit_diLep.GetPull(2);
          Pull_Three = nothingfit_diLep.GetPull(3);
          Pull_Four = nothingfit_diLep.GetPull(4);
          Pull_Five = nothingfit_diLep.GetPull(5);
          Pull_Six = nothingfit_diLep.GetPull(6);
          Pull_Seven = nothingfit_diLep.GetPull(7);
          Pull_Eight = nothingfit_diLep.GetPull(8);
          Pull_Nine = nothingfit_diLep.GetPull(9);
          
          Eg_nothingfit_diLep = nothingfit_diLep.FitPhotonEnergy();
          
          for(int i = 0; i < 3; i++) p4nothing[i] = nothingfit_diLep.FitP4(i);
          
          P4pho_nothingfit_diLep.SetPxPyPzE(0.0,0.0,Eg_nothingfit_diLep,Eg_nothingfit_diLep);
          P4pro_nothingfit_diLep = p4nothing[0];
          P4ep_nothingfit_diLep = p4nothing[1];
          P4em_nothingfit_diLep = p4nothing[2];
          
          TLorentzVector MM_PVec_nothingfit_diLep = (P4pho_nothingfit_diLep + P4target) - P4pro_nothingfit_diLep;
          TLorentzVector MM_PEpEmVec_nothingfit_diLep = (P4pho_nothingfit_diLep + P4target) - (P4pro_nothingfit_diLep + P4ep_nothingfit_diLep + P4em_nothingfit_diLep);
          TLorentzVector IV_EpEmVec_nothingfit_diLep = P4ep_nothingfit_diLep + P4em_nothingfit_diLep;
          open_angle_nothingfit_diLep = P4ep_nothingfit_diLep.Vect().Angle(P4em_nothingfit_diLep.Vect());
          IV_EpEm_nothingfit_diLep_P = sqrt(pow(IV_EpEmVec_nothingfit_diLep.Px(),2) + pow(IV_EpEmVec_nothingfit_diLep.Py(),2) + pow(IV_EpEmVec_nothingfit_diLep.Pz(),2));
          
          CM_Theta_nothingfit_diLep = OneBoostComparison(MM_PVec_nothingfit_diLep, (P4pho_nothingfit_diLep + P4target));
          mandelstam_t_nothingfit_diLep = (P4target - P4pro_nothingfit_diLep).M2();
          
          
          mE_PEpEm_nothingfit_diLep = MM_PEpEmVec_nothingfit_diLep.E();
          mm2_PEpEm_nothingfit_diLep = MM_PEpEmVec_nothingfit_diLep.M2();
          mm_PEpEm_nothingfit_diLep = MM_PEpEmVec_nothingfit_diLep.M();
          
          mm_P_nothingfit_diLep = MM_PVec_nothingfit_diLep.M();
          mm2_P_nothingfit_diLep = MM_PVec_nothingfit_diLep.M2();
          IV_EpEm_nothingfit_diLep = IV_EpEmVec_nothingfit_diLep.M();
          
          P_px_nothingfit_diLep = P4pro_nothingfit_diLep.Px();
          P_py_nothingfit_diLep = P4pro_nothingfit_diLep.Py();
          P_pz_nothingfit_diLep = P4pro_nothingfit_diLep.Pz();
          P_Theta_nothingfit_diLep = P4pro_nothingfit_diLep.Theta()*180./TMath::Pi();
          P_Phi_nothingfit_diLep = P4pro_nothingfit_diLep.Phi()*180./TMath::Pi();
          
          Ep_px_nothingfit_diLep = P4ep_nothingfit_diLep.Px();
          Ep_py_nothingfit_diLep = P4ep_nothingfit_diLep.Py();
          Ep_pz_nothingfit_diLep = P4ep_nothingfit_diLep.Pz();
          Ep_Theta_nothingfit_diLep = P4ep_nothingfit_diLep.Theta()*180./TMath::Pi();
          Ep_Phi_nothingfit_diLep = P4ep_nothingfit_diLep.Phi()*180./TMath::Pi();
          
          Em_px_nothingfit_diLep = P4em_nothingfit_diLep.Px();
          Em_py_nothingfit_diLep = P4em_nothingfit_diLep.Py();
          Em_pz_nothingfit_diLep = P4em_nothingfit_diLep.Pz();
          Em_Theta_nothingfit_diLep = P4em_nothingfit_diLep.Theta()*180./TMath::Pi();
          Em_Phi_nothingfit_diLep = P4em_nothingfit_diLep.Phi()*180./TMath::Pi();
          
          P_Ptot_nothingfit_diLep = sqrt(pow(P_px_nothingfit_diLep,2) + pow(P_py_nothingfit_diLep,2) + pow(P_pz_nothingfit_diLep,2));
          Ep_Ptot_nothingfit_diLep = sqrt(pow(Ep_px_nothingfit_diLep,2) + pow(Ep_py_nothingfit_diLep,2) + pow(Ep_pz_nothingfit_diLep,2));
          Em_Ptot_nothingfit_diLep = sqrt(pow(Em_px_nothingfit_diLep,2) + pow(Em_py_nothingfit_diLep,2) + pow(Em_pz_nothingfit_diLep,2));
          
          
          t4->Fill();
        }//end of exclusive cut loop
      }
    }
  }
  t4->Write();
  
  cout<<"######-------------------------------------#####"<<endl;
  cout<<Good_events<<" are GOOD "<<Bad_events<<" are BAD "<<endl;
  cout<<"######-------------------------------------#####"<<endl;
  
  outFile.Write(); // write to the output file
  outFile.Close(); // close the output file
  
}//end of main


