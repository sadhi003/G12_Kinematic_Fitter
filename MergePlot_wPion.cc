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

#include"/Volumes/Mac_Storage/Work_Data/G12_NECCESSITIES/g12_corrections/g12_corrections/g12_corrections.hpp"


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

  cout<<"here"<<endl;
  //extern int optind;
  //gROOT->Reset();
  //TROOT troot();
  cout<<"here"<<endl;

  Float_t prot_P, prot_Theta, prot_Phi, egam, tpho;
  Float_t pip_P, pip_Theta, pip_Phi;
  Float_t pim_P, pim_Theta, pim_Phi;
  Float_t real_pip_P, real_pip_Theta, real_pip_Phi;
  Float_t real_pim_P, real_pim_Theta, real_pim_Phi;
  Float_t prot_mat[5][5], pim_mat[5][5], pip_mat[5][5];

  
  TChain *chain = new TChain("lepTree");
  chain->Add("IsLepG7_g12.root");
  
  Int_t event, run, trigBits;
  chain->SetBranchAddress("event",&event);
  chain->SetBranchAddress("run",&run);
  chain->SetBranchAddress("trigBits",&trigBits);

  
  chain->SetBranchAddress("egam",&egam);
  chain->SetBranchAddress("tpho",&tpho);
  
  chain->SetBranchAddress("prot_P",&prot_P);
  chain->SetBranchAddress("prot_Theta",&prot_Theta);
  chain->SetBranchAddress("prot_Phi",&prot_Phi);
  
  //for electron/positron 3-momenta
  
  chain->SetBranchAddress("pip_P",&pip_P);
  chain->SetBranchAddress("pip_Theta",&pip_Theta);
  chain->SetBranchAddress("pip_Phi",&pip_Phi);
  
  chain->SetBranchAddress("pim_P",&pim_P);
  chain->SetBranchAddress("pim_Theta",&pim_Theta);
  chain->SetBranchAddress("pim_Phi",&pim_Phi);
  //end electron/positron 3-momenta
  
  //for pion 3-momenta
  chain->SetBranchAddress("real_pip_P",&real_pip_P);
  chain->SetBranchAddress("real_pip_Theta",&real_pip_Theta);
  chain->SetBranchAddress("real_pip_Phi",&real_pip_Phi);
  
  chain->SetBranchAddress("real_pim_P",&real_pim_P);
  chain->SetBranchAddress("real_pim_Theta",&real_pim_Theta);
  chain->SetBranchAddress("real_pim_Phi",&real_pim_Phi);
  //End pion 3-momenta
  
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
  
  Double_t mE_PEpEm, mm2_PEpEm, mm2_P, E_g, IV_EpEm, IV_EpEm_P, mandelstam_t, CM_Theta, open_angle, open_angle_diPion, IV_EpEm_P_diPion;
  
  t4->Branch("mE_PEpEm",&mE_PEpEm,"mE_PEpEm/D");
  t4->Branch("mm2_PEpEm",&mm2_PEpEm,"mm2_PEpEm/D");
  t4->Branch("mm2_P",&mm2_P,"mm2_P/D");
  t4->Branch("IV_EpEm",&IV_EpEm,"IV_EpEm/D");
  t4->Branch("IV_EpEm_P",&IV_EpEm_P,"IV_EpEm_P/D");
  t4->Branch("mandelstam_t",&mandelstam_t,"mandelstam_t/D");
  t4->Branch("CM_Theta",&CM_Theta,"CM_Theta/D");
  t4->Branch("open_angle",&open_angle,"open_angle/D");

  t4->Branch("open_angle_diPion",&open_angle_diPion,"open_angle_diPion/D");
  t4->Branch("IV_EpEm_P_diPion",&IV_EpEm_P_diPion,"IV_EpEm_P_diPion/D");
  
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
  
  Int_t Pip_geofid, Pim_geofid;
  t4->Branch("Pim_geofid",&Pim_geofid,"Pim_geofid/I");
  t4->Branch("Pip_geofid",&Pip_geofid,"Pip_geofid/I");
  
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
  
  /*############################################################################*
   *############################################################################*
   *################## Lets Start with pi+pi- kinematic fitting ################*
   *##################                3 fits                    ################*
   *##################              pi+pi-(gamma)  1-C          ################*
   *##################              pi+pi-(pi0)    1-C          ################*
   *##################              pi+pi-(0)      4-C          ################*
   *############################################################################*
   *############################################################################*/
  //pi+pi-(gamma) fit
  Double_t Pull_Chi_gamfit_diPion, Pull_Prob_gamfit_diPion;
  
  t4->Branch("Pull_Chi_gamfit_diPion",&Pull_Chi_gamfit_diPion,"Pull_Chi_gamfit_diPion/D");
  t4->Branch("Pull_Prob_gamfit_diPion",&Pull_Prob_gamfit_diPion,"Pull_Prob_gamfit_diPion/D");
  
  Double_t mE_PEpEm_gamfit_diPion, mm2_PEpEm_gamfit_diPion, mm_PEpEm_gamfit_diPion, mm_P_gamfit_diPion, mm2_P_gamfit_diPion, Eg_gamfit_diPion, IV_EpEm_gamfit_diPion, IV_EpEm_gamfit_diPion_P;
  
  t4->Branch("mE_PEpEm_gamfit_diPion",&mE_PEpEm_gamfit_diPion,"mE_PEpEm_gamfit_diPion/D");
  t4->Branch("mm2_PEpEm_gamfit_diPion",&mm2_PEpEm_gamfit_diPion,"mm2_PEpEm_gamfit_diPion/D");
  t4->Branch("mm_PEpEm_gamfit_diPion",&mm_PEpEm_gamfit_diPion,"mm_PEpEm_gamfit_diPion/D");
  t4->Branch("mm_P_gamfit_diPion",&mm_P_gamfit_diPion,"mm_P_gamfit_diPion/D");
  t4->Branch("mm2_P_gamfit_diPion",&mm2_P_gamfit_diPion,"mm2_P_gamfit_diPion/D");
  t4->Branch("IV_EpEm_gamfit_diPion",&IV_EpEm_gamfit_diPion,"IV_EpEm_gamfit_diPion/D");
  t4->Branch("IV_EpEm_gamfit_diPion_P",&IV_EpEm_gamfit_diPion_P,"IV_EpEm_gamfit_diPion_P/D");
  
  t4->Branch("Eg_gamfit_diPion",&Eg_gamfit_diPion,"Eg_gamfit_diPion/D");
  
  Double_t open_angle_gamfit_diPion;
  t4->Branch("open_angle_gamfit_diPion",&open_angle_gamfit_diPion,"open_angle_gamfit_diPion/D");
  
  Double_t P_px_gamfit_diPion, P_py_gamfit_diPion, P_pz_gamfit_diPion, P_Theta_gamfit_diPion, P_Phi_gamfit_diPion, P_Ptot_gamfit_diPion;
  Double_t Ep_px_gamfit_diPion, Ep_py_gamfit_diPion, Ep_pz_gamfit_diPion, Ep_Theta_gamfit_diPion, Ep_Phi_gamfit_diPion, Ep_Ptot_gamfit_diPion;
  Double_t Em_px_gamfit_diPion, Em_py_gamfit_diPion, Em_pz_gamfit_diPion, Em_Theta_gamfit_diPion, Em_Phi_gamfit_diPion, Em_Ptot_gamfit_diPion;
  
  t4->Branch("P_px_gamfit_diPion",&P_px_gamfit_diPion,"P_px_gamfit_diPion/D");
  t4->Branch("P_py_gamfit_diPion",&P_py_gamfit_diPion,"P_py_gamfit_diPion/D");
  t4->Branch("P_pz_gamfit_diPion",&P_pz_gamfit_diPion,"P_pz_gamfit_diPion/D");
  t4->Branch("P_Theta_gamfit_diPion",&P_Theta_gamfit_diPion,"P_Theta_gamfit_diPion/D");
  t4->Branch("P_Phi_gamfit_diPion",&P_Phi_gamfit_diPion,"P_Phi_gamfit_diPion/D");
  t4->Branch("P_Ptot_gamfit_diPion",&P_Ptot_gamfit_diPion,"P_Ptot_gamfit_diPion/D");
  
  t4->Branch("Ep_px_gamfit_diPion",&Ep_px_gamfit_diPion,"Ep_px_gamfit_diPion/D");
  t4->Branch("Ep_py_gamfit_diPion",&Ep_py_gamfit_diPion,"Ep_py_gamfit_diPion/D");
  t4->Branch("Ep_pz_gamfit_diPion",&Ep_pz_gamfit_diPion,"Ep_pz_gamfit_diPion/D");
  t4->Branch("Ep_Theta_gamfit_diPion",&Ep_Theta_gamfit_diPion,"Ep_Theta_gamfit_diPion/D");
  t4->Branch("Ep_Phi_gamfit_diPion",&Ep_Phi_gamfit_diPion,"Ep_Phi_gamfit_diPion/D");
  t4->Branch("Ep_Ptot_gamfit_diPion",&Ep_Ptot_gamfit_diPion,"Ep_Ptot_gamfit_diPion/D");
  
  t4->Branch("Em_px_gamfit_diPion",&Em_px_gamfit_diPion,"Em_px_gamfit_diPion/D");
  t4->Branch("Em_py_gamfit_diPion",&Em_py_gamfit_diPion,"Em_py_gamfit_diPion/D");
  t4->Branch("Em_pz_gamfit_diPion",&Em_pz_gamfit_diPion,"Em_pz_gamfit_diPion/D");
  t4->Branch("Em_Theta_gamfit_diPion",&Em_Theta_gamfit_diPion,"Em_Theta_gamfit_diPion/D");
  t4->Branch("Em_Phi_gamfit_diPion",&Em_Phi_gamfit_diPion,"Em_Phi_gamfit_diPion/D");
  t4->Branch("Em_Ptot_gamfit_diPion",&Em_Ptot_gamfit_diPion,"Em_Ptot_gamfit_diPion/D");
  
  Double_t CM_Theta_gamfit_diPion;
  t4->Branch("CM_Theta_gamfit_diPion",&CM_Theta_gamfit_diPion,"CM_Theta_gamfit_diPion/D");
  
  Double_t mandelstam_t_gamfit_diPion;
  t4->Branch("mandelstam_t_gamfit_diPion",&mandelstam_t_gamfit_diPion,"mandelstam_t_gamfit_diPion/D");
  
  
  //pi+pi-(pi0) fit
  Double_t Pull_Chi_pi0fit_diPion, Pull_Prob_pi0fit_diPion;
  
  t4->Branch("Pull_Chi_pi0fit_diPion",&Pull_Chi_pi0fit_diPion,"Pull_Chi_pi0fit_diPion/D");
  t4->Branch("Pull_Prob_pi0fit_diPion",&Pull_Prob_pi0fit_diPion,"Pull_Prob_pi0fit_diPion/D");
  
  Double_t mE_PEpEm_pi0fit_diPion, mm2_PEpEm_pi0fit_diPion, mm_PEpEm_pi0fit_diPion, mm_P_pi0fit_diPion, mm2_P_pi0fit_diPion, Eg_pi0fit_diPion, IV_EpEm_pi0fit_diPion, IV_EpEm_pi0fit_diPion_P;
  
  t4->Branch("mE_PEpEm_pi0fit_diPion",&mE_PEpEm_pi0fit_diPion,"mE_PEpEm_pi0fit_diPion/D");
  t4->Branch("mm2_PEpEm_pi0fit_diPion",&mm2_PEpEm_pi0fit_diPion,"mm2_PEpEm_pi0fit_diPion/D");
  t4->Branch("mm_PEpEm_pi0fit_diPion",&mm_PEpEm_pi0fit_diPion,"mm_PEpEm_pi0fit_diPion/D");
  t4->Branch("mm_P_pi0fit_diPion",&mm_P_pi0fit_diPion,"mm_P_pi0fit_diPion/D");
  t4->Branch("mm2_P_pi0fit_diPion",&mm2_P_pi0fit_diPion,"mm2_P_pi0fit_diPion/D");
  t4->Branch("IV_EpEm_pi0fit_diPion",&IV_EpEm_pi0fit_diPion,"IV_EpEm_pi0fit_diPion/D");
  t4->Branch("IV_EpEm_pi0fit_diPion_P",&IV_EpEm_pi0fit_diPion_P,"IV_EpEm_pi0fit_diPion_P/D");
  
  t4->Branch("Eg_pi0fit_diPion",&Eg_pi0fit_diPion,"Eg_pi0fit_diPion/D");
  
  Double_t open_angle_pi0fit_diPion;
  t4->Branch("open_angle_pi0fit_diPion",&open_angle_pi0fit_diPion,"open_angle_pi0fit_diPion/D");
  
  Double_t P_px_pi0fit_diPion, P_py_pi0fit_diPion, P_pz_pi0fit_diPion, P_Theta_pi0fit_diPion, P_Phi_pi0fit_diPion, P_Ptot_pi0fit_diPion;
  Double_t Ep_px_pi0fit_diPion, Ep_py_pi0fit_diPion, Ep_pz_pi0fit_diPion, Ep_Theta_pi0fit_diPion, Ep_Phi_pi0fit_diPion, Ep_Ptot_pi0fit_diPion;
  Double_t Em_px_pi0fit_diPion, Em_py_pi0fit_diPion, Em_pz_pi0fit_diPion, Em_Theta_pi0fit_diPion, Em_Phi_pi0fit_diPion, Em_Ptot_pi0fit_diPion;
  
  t4->Branch("P_px_pi0fit_diPion",&P_px_pi0fit_diPion,"P_px_pi0fit_diPion/D");
  t4->Branch("P_py_pi0fit_diPion",&P_py_pi0fit_diPion,"P_py_pi0fit_diPion/D");
  t4->Branch("P_pz_pi0fit_diPion",&P_pz_pi0fit_diPion,"P_pz_pi0fit_diPion/D");
  t4->Branch("P_Theta_pi0fit_diPion",&P_Theta_pi0fit_diPion,"P_Theta_pi0fit_diPion/D");
  t4->Branch("P_Phi_pi0fit_diPion",&P_Phi_pi0fit_diPion,"P_Phi_pi0fit_diPion/D");
  t4->Branch("P_Ptot_pi0fit_diPion",&P_Ptot_pi0fit_diPion,"P_Ptot_pi0fit_diPion/D");
  
  t4->Branch("Ep_px_pi0fit_diPion",&Ep_px_pi0fit_diPion,"Ep_px_pi0fit_diPion/D");
  t4->Branch("Ep_py_pi0fit_diPion",&Ep_py_pi0fit_diPion,"Ep_py_pi0fit_diPion/D");
  t4->Branch("Ep_pz_pi0fit_diPion",&Ep_pz_pi0fit_diPion,"Ep_pz_pi0fit_diPion/D");
  t4->Branch("Ep_Theta_pi0fit_diPion",&Ep_Theta_pi0fit_diPion,"Ep_Theta_pi0fit_diPion/D");
  t4->Branch("Ep_Phi_pi0fit_diPion",&Ep_Phi_pi0fit_diPion,"Ep_Phi_pi0fit_diPion/D");
  t4->Branch("Ep_Ptot_pi0fit_diPion",&Ep_Ptot_pi0fit_diPion,"Ep_Ptot_pi0fit_diPion/D");
  
  t4->Branch("Em_px_pi0fit_diPion",&Em_px_pi0fit_diPion,"Em_px_pi0fit_diPion/D");
  t4->Branch("Em_py_pi0fit_diPion",&Em_py_pi0fit_diPion,"Em_py_pi0fit_diPion/D");
  t4->Branch("Em_pz_pi0fit_diPion",&Em_pz_pi0fit_diPion,"Em_pz_pi0fit_diPion/D");
  t4->Branch("Em_Theta_pi0fit_diPion",&Em_Theta_pi0fit_diPion,"Em_Theta_pi0fit_diPion/D");
  t4->Branch("Em_Phi_pi0fit_diPion",&Em_Phi_pi0fit_diPion,"Em_Phi_pi0fit_diPion/D");
  t4->Branch("Em_Ptot_pi0fit_diPion",&Em_Ptot_pi0fit_diPion,"Em_Ptot_pi0fit_diPion/D");
  
  Double_t CM_Theta_pi0fit_diPion;
  t4->Branch("CM_Theta_pi0fit_diPion",&CM_Theta_pi0fit_diPion,"CM_Theta_pi0fit_diPion/D");
  
  Double_t mandelstam_t_pi0fit_diPion;
  t4->Branch("mandelstam_t_pi0fit_diPion",&mandelstam_t_pi0fit_diPion,"mandelstam_t_pi0fit_diPion/D");
  
  //pi+pi-(0) fit
  Double_t Pull_Chi_nothingfit_diPion, Pull_Prob_nothingfit_diPion;
  
  t4->Branch("Pull_Chi_nothingfit_diPion",&Pull_Chi_nothingfit_diPion,"Pull_Chi_nothingfit_diPion/D");
  t4->Branch("Pull_Prob_nothingfit_diPion",&Pull_Prob_nothingfit_diPion,"Pull_Prob_nothingfit_diPion/D");
  
  Double_t mE_PEpEm_nothingfit_diPion, mm2_PEpEm_nothingfit_diPion, mm_PEpEm_nothingfit_diPion, mm_P_nothingfit_diPion, mm2_P_nothingfit_diPion, Eg_nothingfit_diPion, IV_EpEm_nothingfit_diPion, IV_EpEm_nothingfit_diPion_P;
  
  t4->Branch("mE_PEpEm_nothingfit_diPion",&mE_PEpEm_nothingfit_diPion,"mE_PEpEm_nothingfit_diPion/D");
  t4->Branch("mm2_PEpEm_nothingfit_diPion",&mm2_PEpEm_nothingfit_diPion,"mm2_PEpEm_nothingfit_diPion/D");
  t4->Branch("mm_PEpEm_nothingfit_diPion",&mm_PEpEm_nothingfit_diPion,"mm_PEpEm_nothingfit_diPion/D");
  t4->Branch("mm_P_nothingfit_diPion",&mm_P_nothingfit_diPion,"mm_P_nothingfit_diPion/D");
  t4->Branch("mm2_P_nothingfit_diPion",&mm2_P_nothingfit_diPion,"mm2_P_nothingfit_diPion/D");
  t4->Branch("IV_EpEm_nothingfit_diPion",&IV_EpEm_nothingfit_diPion,"IV_EpEm_nothingfit_diPion/D");
  t4->Branch("IV_EpEm_nothingfit_diPion_P",&IV_EpEm_nothingfit_diPion_P,"IV_EpEm_nothingfit_diPion_P/D");
  
  t4->Branch("Eg_nothingfit_diPion",&Eg_nothingfit_diPion,"Eg_nothingfit_diPion/D");
  
  Double_t open_angle_nothingfit_diPion;
  t4->Branch("open_angle_nothingfit_diPion",&open_angle_nothingfit_diPion,"open_angle_nothingfit_diPion/D");
  
  Double_t P_px_nothingfit_diPion, P_py_nothingfit_diPion, P_pz_nothingfit_diPion, P_Theta_nothingfit_diPion, P_Phi_nothingfit_diPion, P_Ptot_nothingfit_diPion;
  Double_t Ep_px_nothingfit_diPion, Ep_py_nothingfit_diPion, Ep_pz_nothingfit_diPion, Ep_Theta_nothingfit_diPion, Ep_Phi_nothingfit_diPion, Ep_Ptot_nothingfit_diPion;
  Double_t Em_px_nothingfit_diPion, Em_py_nothingfit_diPion, Em_pz_nothingfit_diPion, Em_Theta_nothingfit_diPion, Em_Phi_nothingfit_diPion, Em_Ptot_nothingfit_diPion;
  
  t4->Branch("P_px_nothingfit_diPion",&P_px_nothingfit_diPion,"P_px_nothingfit_diPion/D");
  t4->Branch("P_py_nothingfit_diPion",&P_py_nothingfit_diPion,"P_py_nothingfit_diPion/D");
  t4->Branch("P_pz_nothingfit_diPion",&P_pz_nothingfit_diPion,"P_pz_nothingfit_diPion/D");
  t4->Branch("P_Theta_nothingfit_diPion",&P_Theta_nothingfit_diPion,"P_Theta_nothingfit_diPion/D");
  t4->Branch("P_Phi_nothingfit_diPion",&P_Phi_nothingfit_diPion,"P_Phi_nothingfit_diPion/D");
  t4->Branch("P_Ptot_nothingfit_diPion",&P_Ptot_nothingfit_diPion,"P_Ptot_nothingfit_diPion/D");
  
  t4->Branch("Ep_px_nothingfit_diPion",&Ep_px_nothingfit_diPion,"Ep_px_nothingfit_diPion/D");
  t4->Branch("Ep_py_nothingfit_diPion",&Ep_py_nothingfit_diPion,"Ep_py_nothingfit_diPion/D");
  t4->Branch("Ep_pz_nothingfit_diPion",&Ep_pz_nothingfit_diPion,"Ep_pz_nothingfit_diPion/D");
  t4->Branch("Ep_Theta_nothingfit_diPion",&Ep_Theta_nothingfit_diPion,"Ep_Theta_nothingfit_diPion/D");
  t4->Branch("Ep_Phi_nothingfit_diPion",&Ep_Phi_nothingfit_diPion,"Ep_Phi_nothingfit_diPion/D");
  t4->Branch("Ep_Ptot_nothingfit_diPion",&Ep_Ptot_nothingfit_diPion,"Ep_Ptot_nothingfit_diPion/D");
  
  t4->Branch("Em_px_nothingfit_diPion",&Em_px_nothingfit_diPion,"Em_px_nothingfit_diPion/D");
  t4->Branch("Em_py_nothingfit_diPion",&Em_py_nothingfit_diPion,"Em_py_nothingfit_diPion/D");
  t4->Branch("Em_pz_nothingfit_diPion",&Em_pz_nothingfit_diPion,"Em_pz_nothingfit_diPion/D");
  t4->Branch("Em_Theta_nothingfit_diPion",&Em_Theta_nothingfit_diPion,"Em_Theta_nothingfit_diPion/D");
  t4->Branch("Em_Phi_nothingfit_diPion",&Em_Phi_nothingfit_diPion,"Em_Phi_nothingfit_diPion/D");
  t4->Branch("Em_Ptot_nothingfit_diPion",&Em_Ptot_nothingfit_diPion,"Em_Ptot_nothingfit_diPion/D");
  
  Double_t CM_Theta_nothingfit_diPion;
  t4->Branch("CM_Theta_nothingfit_diPion",&CM_Theta_nothingfit_diPion,"CM_Theta_nothingfit_diPion/D");
  
  Double_t mandelstam_t_nothingfit_diPion;
  t4->Branch("mandelstam_t_nothingfit_diPion",&mandelstam_t_nothingfit_diPion,"mandelstam_t_nothingfit_diPion/D");

  
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

      Ep_geofid = clas::g12::g12_PosParticle_fiducial_cuts(pip_P, pip_Theta, pip_Phi,"nominal");//Geometric Fiducial Cut
      Em_geofid = clas::g12::g12_NegParticle_fiducial_cuts(pim_P, pim_Theta, pim_Phi,"nominal");//Geometric Fiducial Cut
      Pip_geofid = clas::g12::g12_PosParticle_fiducial_cuts(real_pip_P, real_pip_Theta, real_pip_Phi,"nominal");//Geometric Fiducial Cut
      Pim_geofid = clas::g12::g12_NegParticle_fiducial_cuts(real_pim_P, real_pim_Theta, real_pim_Phi,"nominal");//Geometric Fiducial Cut
      P_geofid = clas::g12::g12_PosParticle_fiducial_cuts(prot_P, prot_Theta, prot_Phi,"nominal");//Geometric Fiducial Cut
      
      
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
      mm2_P = (double)MM2P; E_g = (double)egam; IV_EpEm = (double)IVEpEm;
      Vx = (double)vx; Vy = (double)vy; Vz = (double)vz;
      

      P_px = (prot_P*sin(DegToRad*prot_Theta)*cos(DegToRad*prot_Phi));
      P_py = (prot_P*sin(DegToRad*prot_Theta)*sin(DegToRad*prot_Phi));
      P_pz = (prot_P*cos(DegToRad*prot_Theta));
      
      Double_t pip_px = pip_P*sin(DegToRad*pip_Theta)*cos(DegToRad*pip_Phi);
      Double_t pip_py = pip_P*sin(DegToRad*pip_Theta)*sin(DegToRad*pip_Phi);
      Double_t pip_pz = pip_P*cos(DegToRad*pip_Theta);
      
      Double_t pim_px = pim_P*sin(DegToRad*pim_Theta)*cos(DegToRad*pim_Phi);
      Double_t pim_py = pim_P*sin(DegToRad*pim_Theta)*sin(DegToRad*pim_Phi);
      Double_t pim_pz = pim_P*cos(DegToRad*pim_Theta);
      
      Ep_px = pip_px; Ep_py = pip_py; Ep_pz = pip_pz;
      Em_px = pim_px; Em_py = pim_py; Em_pz = pim_pz;
      P_p = prot_P;
      Ep_p = pip_P;
      Em_p = pim_P;
      
      P_Theta = prot_Theta;
      Ep_Theta = pip_Theta;
      Em_Theta = pim_Theta;
      
      P_Phi = prot_Phi;
      Ep_Phi = pip_Phi;
      Em_Phi = pim_Phi;
      
      
      P4pho.SetPxPyPzE(0.0,0.0,egam,egam);
      P4em.SetPxPyPzE(pim_px,pim_py,pim_pz,sqrt(pim_P*pim_P+0.000510998*0.000510998));
      P4pro.SetPxPyPzE(P_px,P_py,P_pz,sqrt(prot_P*prot_P+0.93828*0.93828));
      P4ep.SetPxPyPzE(pip_px,pip_py,pip_pz,sqrt(pip_P*pip_P+0.000510998*0.000510998));
      V3pro.SetXYZ(prot_vx,prot_vy,prot_vz);
      V3ep.SetXYZ(pip_vx,pip_vy,pip_vz);
      V3em.SetXYZ(pim_vx,pim_vy,pim_vz);
      
      P4Pim.SetPxPyPzE(pim_px,pim_py,pim_pz,sqrt(pim_P*pim_P + M_Pi*M_Pi));
      P4Pip.SetPxPyPzE(pip_px,pip_py,pip_pz,sqrt(pip_P*pip_P + M_Pi*M_Pi));
      
      open_angle = P4ep.Vect().Angle(P4em.Vect());
      TLorentzVector IV_EpEmVec_nofit = P4em + P4ep;
      IV_EpEm_P= sqrt(pow(IV_EpEmVec_nofit.Px(),2) + pow(IV_EpEmVec_nofit.Py(),2) + pow(IV_EpEmVec_nofit.Pz(),2));
      
      CM_Theta = OneBoostComparison((P4pho + P4target - P4pro), (P4pho + P4target));
      mandelstam_t = (P4target - P4pro).M2();
      
      
      TMatrixD covTrack(10,10);
      Double_t c00 = pow((0.001*5.715),2)/3;
      covTrack(0,0)=c00;
      covTrack(1,1)= prot_mat[0][0]*pow(prot_P,4);
      covTrack(1,2)=-prot_mat[0][1]*pow(prot_P,2);
      covTrack(1,3)= -prot_mat[0][2]*pow(prot_P,2);
      covTrack(2,1)= -prot_mat[0][1]*pow(prot_P,2);
      covTrack(2,2)= prot_mat[1][1];
      covTrack(2,3)=prot_mat[1][2];
      covTrack(3,1)= -prot_mat[0][2]*pow(prot_P,2);
      covTrack(3,2)= prot_mat[1][2];
      covTrack(3,3)=prot_mat[2][2];
      
      covTrack(4,4)= pip_mat[0][0]*pow(pip_P,4);
      covTrack(4,5)= -pip_mat[0][1]*pow(pip_P,2);
      covTrack(4,6)= -pip_mat[0][2]*pow(pip_P,2);
      covTrack(5,4)= -pip_mat[0][1]*pow(pip_P,2);
      covTrack(5,5)= pip_mat[1][1];
      covTrack(5,6)= pip_mat[1][2];
      covTrack(6,4)= -pip_mat[0][2]*pow(pip_P,2);
      covTrack(6,5)= pip_mat[1][2];
      covTrack(6,6)= pip_mat[2][2];
      
      covTrack(7,7)= pim_mat[0][0]*pow(pim_P,4);
      covTrack(7,8)= pim_mat[0][1]*pow(pim_P,2);
      covTrack(7,9)= pim_mat[0][2]*pow(pim_P,2);
      covTrack(8,7)= pim_mat[0][1]*pow(pim_P,2);
      covTrack(8,8)= pim_mat[1][1];
      covTrack(8,9)= pim_mat[1][2];
      covTrack(9,7)= pim_mat[0][2]*pow(pim_P,2);
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
      Double_t e_gamma = egam;
      
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
      
      //For pi+pi-(gamma) fitting
      P4Pim.SetPxPyPzE(pim_px,pim_py,pim_pz,sqrt(pim_P*pim_P + M_Pi*M_Pi));
      P4Pip.SetPxPyPzE(pip_px,pip_py,pip_pz,sqrt(pip_P*pip_P + M_Pi*M_Pi));
      
      open_angle_diPion = P4Pip.Vect().Angle(P4Pim.Vect());
      TLorentzVector IV_PipPimVec_nofit = P4Pim + P4Pip;
      IV_EpEm_P_diPion= sqrt(pow(IV_PipPimVec_nofit.Px(),2) + pow(IV_PipPimVec_nofit.Py(),2) + pow(IV_PipPimVec_nofit.Pz(),2));
      
      
      TMatrixD picovTrack(10,10);
      picovTrack(0,0)=c00;
      picovTrack(1,1)= prot_mat[0][0]*pow(prot_P,4);
      picovTrack(1,2)=-prot_mat[0][1]*pow(prot_P,2);
      picovTrack(1,3)= -prot_mat[0][2]*pow(prot_P,2);
      picovTrack(2,1)= -prot_mat[0][1]*pow(prot_P,2);
      picovTrack(2,2)= prot_mat[1][1];
      picovTrack(2,3)=prot_mat[1][2];
      picovTrack(3,1)= -prot_mat[0][2]*pow(prot_P,2);
      picovTrack(3,2)= prot_mat[1][2];
      picovTrack(3,3)=prot_mat[2][2];
      
      picovTrack(4,4)= pip_mat[0][0]*pow(real_pip_P,4);
      picovTrack(4,5)= -pip_mat[0][1]*pow(real_pip_P,2);
      picovTrack(4,6)= -pip_mat[0][2]*pow(real_pip_P,2);
      picovTrack(5,4)= -pip_mat[0][1]*pow(real_pip_P,2);
      picovTrack(5,5)= pip_mat[1][1];
      picovTrack(5,6)= pip_mat[1][2];
      picovTrack(6,4)= -pip_mat[0][2]*pow(real_pip_P,2);
      picovTrack(6,5)= pip_mat[1][2];
      picovTrack(6,6)= pip_mat[2][2];
      
      picovTrack(7,7)= pim_mat[0][0]*pow(real_pim_P,4);
      picovTrack(7,8)= pim_mat[0][1]*pow(real_pim_P,2);
      picovTrack(7,9)= pim_mat[0][2]*pow(real_pim_P,2);
      picovTrack(8,7)= pim_mat[0][1]*pow(real_pim_P,2);
      picovTrack(8,8)= pim_mat[1][1];
      picovTrack(8,9)= pim_mat[1][2];
      picovTrack(9,7)= pim_mat[0][2]*pow(real_pim_P,2);
      picovTrack(9,8)= pim_mat[1][2];
      picovTrack(9,9)= pim_mat[2][2];
      
      std::vector<TLorentzVector> pip4(num_parts);
      std::vector<TLorentzVector> pip4pi0(num_parts);
      std::vector<TLorentzVector> pip4nothing(num_parts);
      
      std::vector<string> piparticles(num_parts);
  
      pip4[0] = P4pro;  pip4pi0[0] = P4pro;  pip4nothing[0] = P4pro;
      pip4[1] = P4Pip;  pip4pi0[1] = P4Pip;  pip4nothing[1] = P4Pip;
      pip4[2] = P4Pim;  pip4pi0[2] = P4Pim;  pip4nothing[2] = P4Pim;
      
      piparticles[0] = "p";
      piparticles[1] = "pi+";
      piparticles[2] = "pi-";
      
      TMatrixD picovMatrix(10,10);
      
      picovMatrix = CorrectCLAS_V(picovTrack,piparticles,pip4,vert,multi,is_mc,experiment);
      Kstream gamfit_diPion;
      gamfit_diPion.StringNames(piparticles);
      gamfit_diPion.FitInput(e_gamma,pip4,picovMatrix,m_targ);
      gamfit_diPion.Fit("gamma");
      
      Pull_Chi_gamfit_diPion = gamfit_diPion.Chi2();
      Pull_Prob_gamfit_diPion = gamfit_diPion.Prob();
      
      Eg_gamfit_diPion = gamfit_diPion.FitPhotonEnergy();
      
      for(int i = 0; i < 3; i++) p4[i] = gamfit_diPion.FitP4(i);
      
      P4pho_gamfit_diPion.SetPxPyPzE(0.0,0.0,Eg_gamfit_diPion,Eg_gamfit_diPion);
      P4pro_gamfit_diPion = p4[0];
      P4ep_gamfit_diPion = p4[1];
      P4em_gamfit_diPion = p4[2];
      
      TLorentzVector MM_PVec_diPion = (P4pho_gamfit_diPion + P4target) - P4pro_gamfit_diPion;
      TLorentzVector MM_PEpEmVec_diPion = (P4pho_gamfit_diPion + P4target) - (P4pro_gamfit_diPion + P4ep_gamfit_diPion + P4em_gamfit_diPion);
      TLorentzVector IV_EpEmVec_diPion = P4ep_gamfit_diPion + P4em_gamfit_diPion;
      open_angle_gamfit_diPion = P4ep_gamfit_diPion.Vect().Angle(P4em_gamfit_diPion.Vect());
      IV_EpEm_gamfit_diPion_P = sqrt(pow(IV_EpEmVec.Px(),2) + pow(IV_EpEmVec.Py(),2) + pow(IV_EpEmVec.Pz(),2));
      
      CM_Theta_gamfit_diPion = OneBoostComparison(MM_PVec, (P4pho_gamfit_diPion + P4target));
      mandelstam_t_gamfit_diPion = (P4target - P4pro_gamfit_diPion).M2();
      
      
      mE_PEpEm_gamfit_diPion = MM_PEpEmVec_diPion.E();
      mm2_PEpEm_gamfit_diPion = MM_PEpEmVec_diPion.M2();
      mm_PEpEm_gamfit_diPion = MM_PEpEmVec_diPion.M();
      
      mm_P_gamfit_diPion = MM_PVec_diPion.M();
      mm2_P_gamfit_diPion = MM_PVec_diPion.M2();
      IV_EpEm_gamfit_diPion = IV_EpEmVec_diPion.M();
      
      P_px_gamfit_diPion = P4pro_gamfit_diPion.Px();
      P_py_gamfit_diPion = P4pro_gamfit_diPion.Py();
      P_pz_gamfit_diPion = P4pro_gamfit_diPion.Pz();
      P_Theta_gamfit_diPion = P4pro_gamfit_diPion.Theta()*180./TMath::Pi();
      P_Phi_gamfit_diPion = P4pro_gamfit_diPion.Phi()*180./TMath::Pi();
      
      Ep_px_gamfit_diPion = P4ep_gamfit_diPion.Px();
      Ep_py_gamfit_diPion = P4ep_gamfit_diPion.Py();
      Ep_pz_gamfit_diPion = P4ep_gamfit_diPion.Pz();
      Ep_Theta_gamfit_diPion = P4ep_gamfit_diPion.Theta()*180./TMath::Pi();
      Ep_Phi_gamfit_diPion = P4ep_gamfit_diPion.Phi()*180./TMath::Pi();
      
      Em_px_gamfit_diPion = P4em_gamfit_diPion.Px();
      Em_py_gamfit_diPion = P4em_gamfit_diPion.Py();
      Em_pz_gamfit_diPion = P4em_gamfit_diPion.Pz();
      Em_Theta_gamfit_diPion = P4em_gamfit_diPion.Theta()*180./TMath::Pi();
      Em_Phi_gamfit_diPion = P4em_gamfit_diPion.Phi()*180./TMath::Pi();
      
      P_Ptot_gamfit_diPion = sqrt(pow(P_px_gamfit_diPion,2) + pow(P_py_gamfit_diPion,2) + pow(P_pz_gamfit_diPion,2));
      Ep_Ptot_gamfit_diPion = sqrt(pow(Ep_px_gamfit_diPion,2) + pow(Ep_py_gamfit_diPion,2) + pow(Ep_pz_gamfit_diPion,2));
      Em_Ptot_gamfit_diPion = sqrt(pow(Em_px_gamfit_diPion,2) + pow(Em_py_gamfit_diPion,2) + pow(Em_pz_gamfit_diPion,2));
      
      //For pi+pi-(pi0) fitting
      TMatrixD picovMatrixpi0(10,10);
      
      picovMatrixpi0 = CorrectCLAS_V(picovTrack,piparticles,pip4pi0,vert,multi,is_mc,experiment);
      Kstream pi0fit_diPion;
      pi0fit_diPion.StringNames(piparticles);
      pi0fit_diPion.FitInput(e_gamma,pip4pi0,picovMatrixpi0,m_targ);
      pi0fit_diPion.Fit("pi0");
      
      
      Pull_Chi_pi0fit_diPion = pi0fit_diPion.Chi2();
      Pull_Prob_pi0fit_diPion = pi0fit_diPion.Prob();
      
      Eg_pi0fit_diPion = pi0fit_diPion.FitPhotonEnergy();
      
      for(int i = 0; i < 3; i++) p4pi0[i] = pi0fit_diPion.FitP4(i);
      
      P4pho_pi0fit_diPion.SetPxPyPzE(0.0,0.0,Eg_pi0fit_diPion,Eg_pi0fit_diPion);
      P4pro_pi0fit_diPion = p4pi0[0];
      P4ep_pi0fit_diPion = p4pi0[1];
      P4em_pi0fit_diPion = p4pi0[2];
      
      TLorentzVector MM_PVec_pi0fit_diPion = (P4pho_pi0fit_diPion + P4target) - P4pro_pi0fit_diPion;
      TLorentzVector MM_PEpEmVec_pi0fit_diPion = (P4pho_pi0fit_diPion + P4target) - (P4pro_pi0fit_diPion + P4ep_pi0fit_diPion + P4em_pi0fit_diPion);
      TLorentzVector IV_EpEmVec_pi0fit_diPion = P4ep_pi0fit_diPion + P4em_pi0fit_diPion;
      open_angle_pi0fit_diPion = P4ep_pi0fit_diPion.Vect().Angle(P4em_pi0fit_diPion.Vect());
      IV_EpEm_pi0fit_diPion_P = sqrt(pow(IV_EpEmVec_pi0fit_diPion.Px(),2) + pow(IV_EpEmVec_pi0fit_diPion.Py(),2) + pow(IV_EpEmVec_pi0fit_diPion.Pz(),2));
      
      CM_Theta_pi0fit_diPion = OneBoostComparison(MM_PVec_pi0fit_diPion, (P4pho_pi0fit_diPion + P4target));
      mandelstam_t_pi0fit_diPion = (P4target - P4pro_pi0fit_diPion).M2();
      
      
      mE_PEpEm_pi0fit_diPion = MM_PEpEmVec_pi0fit_diPion.E();
      mm2_PEpEm_pi0fit_diPion = MM_PEpEmVec_pi0fit_diPion.M2();
      mm_PEpEm_pi0fit_diPion = MM_PEpEmVec_pi0fit_diPion.M();
      
      mm_P_pi0fit_diPion = MM_PVec_pi0fit_diPion.M();
      mm2_P_pi0fit_diPion = MM_PVec_pi0fit_diPion.M2();
      IV_EpEm_pi0fit_diPion = IV_EpEmVec_pi0fit_diPion.M();
      
      P_px_pi0fit_diPion = P4pro_pi0fit_diPion.Px();
      P_py_pi0fit_diPion = P4pro_pi0fit_diPion.Py();
      P_pz_pi0fit_diPion = P4pro_pi0fit_diPion.Pz();
      P_Theta_pi0fit_diPion = P4pro_pi0fit_diPion.Theta()*180./TMath::Pi();
      P_Phi_pi0fit_diPion = P4pro_pi0fit_diPion.Phi()*180./TMath::Pi();
      
      Ep_px_pi0fit_diPion = P4ep_pi0fit_diPion.Px();
      Ep_py_pi0fit_diPion = P4ep_pi0fit_diPion.Py();
      Ep_pz_pi0fit_diPion = P4ep_pi0fit_diPion.Pz();
      Ep_Theta_pi0fit_diPion = P4ep_pi0fit_diPion.Theta()*180./TMath::Pi();
      Ep_Phi_pi0fit_diPion = P4ep_pi0fit_diPion.Phi()*180./TMath::Pi();
      
      Em_px_pi0fit_diPion = P4em_pi0fit_diPion.Px();
      Em_py_pi0fit_diPion = P4em_pi0fit_diPion.Py();
      Em_pz_pi0fit_diPion = P4em_pi0fit_diPion.Pz();
      Em_Theta_pi0fit_diPion = P4em_pi0fit_diPion.Theta()*180./TMath::Pi();
      Em_Phi_pi0fit_diPion = P4em_pi0fit_diPion.Phi()*180./TMath::Pi();
      
      P_Ptot_pi0fit_diPion = sqrt(pow(P_px_pi0fit_diPion,2) + pow(P_py_pi0fit_diPion,2) + pow(P_pz_pi0fit_diPion,2));
      Ep_Ptot_pi0fit_diPion = sqrt(pow(Ep_px_pi0fit_diPion,2) + pow(Ep_py_pi0fit_diPion,2) + pow(Ep_pz_pi0fit_diPion,2));
      Em_Ptot_pi0fit_diPion = sqrt(pow(Em_px_pi0fit_diPion,2) + pow(Em_py_pi0fit_diPion,2) + pow(Em_pz_pi0fit_diPion,2));
      
      //For pi+pi-(0) fitting
      TMatrixD picovMatrixnothing(10,10);
      
      picovMatrixnothing = CorrectCLAS_V(picovTrack,piparticles,pip4nothing,vert,multi,is_mc,experiment);
      Kstream nothingfit_diPion;
      nothingfit_diPion.StringNames(piparticles);
      nothingfit_diPion.FitInput(e_gamma,pip4nothing,picovMatrix,m_targ);
      nothingfit_diPion.Fit();
      
      
      Pull_Chi_nothingfit_diPion = nothingfit_diPion.Chi2();
      Pull_Prob_nothingfit_diPion = nothingfit_diPion.Prob();
      
      Eg_nothingfit_diPion = nothingfit_diPion.FitPhotonEnergy();
      
      for(int i = 0; i < 3; i++) p4nothing[i] = nothingfit_diPion.FitP4(i);
      
      P4pho_nothingfit_diPion.SetPxPyPzE(0.0,0.0,Eg_nothingfit_diPion,Eg_nothingfit_diPion);
      P4pro_nothingfit_diPion = p4nothing[0];
      P4ep_nothingfit_diPion = p4nothing[1];
      P4em_nothingfit_diPion = p4nothing[2];
      
      TLorentzVector MM_PVec_nothingfit_diPion = (P4pho_nothingfit_diPion + P4target) - P4pro_nothingfit_diPion;
      TLorentzVector MM_PEpEmVec_nothingfit_diPion = (P4pho_nothingfit_diPion + P4target) - (P4pro_nothingfit_diPion + P4ep_nothingfit_diPion + P4em_nothingfit_diPion);
      TLorentzVector IV_EpEmVec_nothingfit_diPion = P4ep_nothingfit_diPion + P4em_nothingfit_diPion;
      open_angle_nothingfit_diPion = P4ep_nothingfit_diPion.Vect().Angle(P4em_nothingfit_diPion.Vect());
      IV_EpEm_nothingfit_diPion_P = sqrt(pow(IV_EpEmVec_nothingfit_diPion.Px(),2) + pow(IV_EpEmVec_nothingfit_diPion.Py(),2) + pow(IV_EpEmVec_nothingfit_diPion.Pz(),2));
      
      CM_Theta_nothingfit_diPion = OneBoostComparison(MM_PVec_nothingfit_diPion, (P4pho_nothingfit_diPion + P4target));
      mandelstam_t_nothingfit_diPion = (P4target - P4pro_nothingfit_diPion).M2();
      
      
      mE_PEpEm_nothingfit_diPion = MM_PEpEmVec_nothingfit_diPion.E();
      mm2_PEpEm_nothingfit_diPion = MM_PEpEmVec_nothingfit_diPion.M2();
      mm_PEpEm_nothingfit_diPion = MM_PEpEmVec_nothingfit_diPion.M();
      
      mm_P_nothingfit_diPion = MM_PVec_nothingfit_diPion.M();
      mm2_P_nothingfit_diPion = MM_PVec_nothingfit_diPion.M2();
      IV_EpEm_nothingfit_diPion = IV_EpEmVec_nothingfit_diPion.M();
      
      P_px_nothingfit_diPion = P4pro_nothingfit_diPion.Px();
      P_py_nothingfit_diPion = P4pro_nothingfit_diPion.Py();
      P_pz_nothingfit_diPion = P4pro_nothingfit_diPion.Pz();
      P_Theta_nothingfit_diPion = P4pro_nothingfit_diPion.Theta()*180./TMath::Pi();
      P_Phi_nothingfit_diPion = P4pro_nothingfit_diPion.Phi()*180./TMath::Pi();
      
      Ep_px_nothingfit_diPion = P4ep_nothingfit_diPion.Px();
      Ep_py_nothingfit_diPion = P4ep_nothingfit_diPion.Py();
      Ep_pz_nothingfit_diPion = P4ep_nothingfit_diPion.Pz();
      Ep_Theta_nothingfit_diPion = P4ep_nothingfit_diPion.Theta()*180./TMath::Pi();
      Ep_Phi_nothingfit_diPion = P4ep_nothingfit_diPion.Phi()*180./TMath::Pi();
      
      Em_px_nothingfit_diPion = P4em_nothingfit_diPion.Px();
      Em_py_nothingfit_diPion = P4em_nothingfit_diPion.Py();
      Em_pz_nothingfit_diPion = P4em_nothingfit_diPion.Pz();
      Em_Theta_nothingfit_diPion = P4em_nothingfit_diPion.Theta()*180./TMath::Pi();
      Em_Phi_nothingfit_diPion = P4em_nothingfit_diPion.Phi()*180./TMath::Pi();
      
      P_Ptot_nothingfit_diPion = sqrt(pow(P_px_nothingfit_diPion,2) + pow(P_py_nothingfit_diPion,2) + pow(P_pz_nothingfit_diPion,2));
      Ep_Ptot_nothingfit_diPion = sqrt(pow(Ep_px_nothingfit_diPion,2) + pow(Ep_py_nothingfit_diPion,2) + pow(Ep_pz_nothingfit_diPion,2));
      Em_Ptot_nothingfit_diPion = sqrt(pow(Em_px_nothingfit_diPion,2) + pow(Em_py_nothingfit_diPion,2) + pow(Em_pz_nothingfit_diPion,2));
      
      t4->Fill();
    }
    //t4->Fill();
    
  }
  t4->Write();
  
  cout<<"######-------------------------------------#####"<<endl;
  cout<<Good_events<<" are GOOD "<<Bad_events<<" are BAD "<<endl;
  cout<<"######-------------------------------------#####"<<endl;
  
  outFile.Write(); // write to the output file
  outFile.Close(); // close the output file
  
}//end of main


