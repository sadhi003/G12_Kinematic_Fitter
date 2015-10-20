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

#include"/Volumes/Mac_Storage/Work_Data/G12_NECCESSITIES/trip.hpp"

#include "/Volumes/Mac_Storage/Work_Data/G12_NECCESSITIES/g12_corrections/MK_ORIGINAL/g12_TOF_knockout.hpp"

#include "/Volumes/Mac_Storage/Work_Data/G12_NECCESSITIES/g12_corrections/MK_ORIGINAL/g12_EC_knockout.hpp"
#include "/Volumes/Mac_Storage/Work_Data/G12_NECCESSITIES/g12_corrections/MK_ORIGINAL/g12_ECxyz_2uvw.hpp"
#include "/Volumes/Mac_Storage/Work_Data/G12_NECCESSITIES/g12_corrections/MK_ORIGINAL/g12_EC_fiducial.hpp"

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

//int main(int __argc,char *__argv[]){
//int main(){
void Fiducial_Clone_Tree(TString file_in, TString file_out) {
  

  //char *outFileName = (char *) "Fitted_Pi0_Conversion.root";
  extern int optind;
  gROOT->Reset();
  TROOT troot();
  
  Float_t prot_P, prot_Theta, prot_Phi, egam, tpho;
  Float_t pip_P, pip_Theta, pip_Phi;
  Float_t pim_P, pim_Theta, pim_Phi;
  
  Float_t prot_mat[5][5], pim_mat[5][5], pip_mat[5][5];
  
  //  TFile inFile(__argv[n_arg]); // open the input file
  //
  //  TTree *Lep = (TTree*)inFile.Get("Lep")
  
  TChain *chain = new TChain("lepTree");
  
  //chain->Add(__argv[1]);
  chain->Add(file_in);
  
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
  
  Int_t pip_EC_hit, pip_CC_hit, pim_EC_hit, pim_CC_hit, pip_IsLepG7, pim_IsLepG7, pip_IsLepG7ect, pim_IsLepG7ect, pip_IsLepG7cc, pim_IsLepG7cc;
  
  chain->SetBranchAddress("pip_EC_hit",&pip_EC_hit);
  chain->SetBranchAddress("pip_CC_hit",&pip_CC_hit);
  chain->SetBranchAddress("pim_EC_hit",&pim_EC_hit);
  chain->SetBranchAddress("pim_CC_hit",&pim_CC_hit);
  chain->SetBranchAddress("pip_IsLepG7",&pip_IsLepG7);
  chain->SetBranchAddress("pim_IsLepG7",&pim_IsLepG7);
  chain->SetBranchAddress("pip_IsLepG7ect",&pip_IsLepG7ect);
  chain->SetBranchAddress("pip_IsLepG7cc",&pip_IsLepG7cc);
  chain->SetBranchAddress("pim_IsLepG7ect",&pim_IsLepG7ect);
  chain->SetBranchAddress("pim_IsLepG7cc",&pim_IsLepG7cc);
  
  float pim_ECx;
  float pim_ECy;
  float pim_ECz;
  float pim_ECin, pim_ECout;
  int pim_CChit_stat;
  
  chain->SetBranchAddress("pim_ECx",&pim_ECx);
  chain->SetBranchAddress("pim_ECy",&pim_ECy);
  chain->SetBranchAddress("pim_ECz",&pim_ECz);
  chain->SetBranchAddress("pim_ECin",&pim_ECin);
  chain->SetBranchAddress("pim_ECout",&pim_ECout);
  chain->SetBranchAddress("pim_CChit_stat",&pim_CChit_stat);
  
  
  float pip_ECx;
  float pip_ECy;
  float pip_ECz;
  float pip_ECin, pip_ECout;
  int pip_CChit_stat;
  
  chain->SetBranchAddress("pip_ECx",&pip_ECx);
  chain->SetBranchAddress("pip_ECy",&pip_ECy);
  chain->SetBranchAddress("pip_ECz",&pip_ECz);
  chain->SetBranchAddress("pip_ECin",&pip_ECin);
  chain->SetBranchAddress("pip_ECout",&pip_ECout);
  chain->SetBranchAddress("pip_CChit_stat",&pip_CChit_stat);
  
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
  
  float pim_dTOF, Timing_pim, pim_Tprop;
  chain->SetBranchAddress("pim_dTOF",&pim_dTOF);
  chain->SetBranchAddress("Timing_pim",&Timing_pim);
  chain->SetBranchAddress("pim_Tprop",&pim_Tprop);
  float pim_scLen, pim_SC_time, pim_ST_time, pim_STSCpl;
  chain->SetBranchAddress("pim_scLen",&pim_scLen);
  chain->SetBranchAddress("pim_SC_time",&pim_SC_time);
  chain->SetBranchAddress("pim_ST_time",&pim_ST_time);
  chain->SetBranchAddress("pim_STSCpl",&pim_STSCpl);
  
  float prot_dTOF, Timing_prot, prot_Tprop;
  chain->SetBranchAddress("prot_dTOF",&prot_dTOF);
  chain->SetBranchAddress("Timing_prot",&Timing_prot);
  chain->SetBranchAddress("prot_Tprop",&prot_Tprop);
  float prot_scLen, prot_SC_time, prot_ST_time, prot_STSCpl;
  chain->SetBranchAddress("prot_scLen",&prot_scLen);
  chain->SetBranchAddress("prot_SC_time",&prot_SC_time);
  chain->SetBranchAddress("prot_ST_time",&prot_ST_time);
  chain->SetBranchAddress("prot_STSCpl",&prot_STSCpl);
  
  TLorentzVector P4pho,P4pro,P4em,P4ep,P4target;
  TLorentzVector P4pho_fit,P4pro_fit,P4pim_fit,P4pip_fit;
  
  P4target.SetPxPyPzE(0.0,0.0,0.0,0.93828);
  
  //TLorentzVector vEg_in, vEp_in, vEm_in, vP_in;
  
  TVector3 V3pro,V3ep,V3em;
  
  //Pion fitted 4vectors
  TLorentzVector P4Pim,P4Pip;
  TLorentzVector P4pho_PiGam_fit,P4pro_PiGam_fit,P4em_PiGam_fit,P4ep_PiGam_fit;
  TLorentzVector P4pho_Pi_fit,P4pro_Pi_fit,P4em_Pi_fit,P4ep_Pi_fit;
  TLorentzVector P4pho_PiPi0_fit,P4pro_PiPi0_fit,P4em_PiPi0_fit,P4ep_PiPi0_fit;
  
  
  Int_t nentries = (Int_t)chain->GetEntries();
  
  tot_events = tot_events + nentries;
  
  cout<<"############################################"<<endl;
  cout<<tot_events<<" TOTAL NUMBER PROCESSED"<<endl;
  cout<<"############################################"<<endl;

  
  TFile outFile(file_out,"recreate");
  TTree *t4 = new TTree("LepTree","LepTree");
  
  Double_t Pull_Chi, Pull_Prob;
  
  t4->Branch("Pull_Chi",&Pull_Chi,"Pull_Chi/D");
  t4->Branch("Pull_Prob",&Pull_Prob,"Pull_Prob/D");
  
  Double_t mE_PEpEm, mm2_PEpEm, mm2_P, E_g, IV_EpEm, IV_EpEm_P;
  
  t4->Branch("mE_PEpEm",&mE_PEpEm,"mE_PEpEm/D");
  t4->Branch("mm2_PEpEm",&mm2_PEpEm,"mm2_PEpEm/D");
  t4->Branch("mm2_P",&mm2_P,"mm2_P/D");
  t4->Branch("IV_EpEm",&IV_EpEm,"IV_EpEm/D");
  t4->Branch("IV_EpEm_P",&IV_EpEm_P,"IV_EpEm_P/D");
  
  
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
  
  
  
  t4->Branch("E_g",&E_g,"E_g/D");

  
  Double_t mE_PEpEm_fit, mm2_PEpEm_fit, mm_PEpEm_fit, mm_P_fit, mm2_P_fit, Eg_fit, IV_EpEm_fit, IV_EpEm_fit_P;
  
  t4->Branch("mE_PEpEm_fit",&mE_PEpEm_fit,"mE_PEpEm_fit/D");
  t4->Branch("mm2_PEpEm_fit",&mm2_PEpEm_fit,"mm2_PEpEm_fit/D");
  t4->Branch("mm_PEpEm_fit",&mm_PEpEm_fit,"mm_PEpEm_fit/D");
  t4->Branch("mm_P_fit",&mm_P_fit,"mm_P_fit/D");
  t4->Branch("mm2_P_fit",&mm2_P_fit,"mm2_P_fit/D");
  t4->Branch("IV_EpEm_fit",&IV_EpEm_fit,"IV_EpEm_fit/D");
  t4->Branch("IV_EpEm_fit_P",&IV_EpEm_fit_P,"IV_EpEm_fit_P/D");
  
  
  t4->Branch("Eg_fit",&Eg_fit,"Eg_fit/D");
  
  Int_t Em_CChit, Em_EChit, Ep_CChit, Ep_EChit;
  t4->Branch("Em_CChit",&Em_CChit,"Em_CChit/I");
  t4->Branch("Em_EChit",&Em_EChit,"Em_EChit/I");
  t4->Branch("Ep_CChit",&Ep_CChit,"Ep_CChit/I");
  t4->Branch("Ep_EChit",&Ep_EChit,"Ep_EChit/I");
  
  
  Int_t Ep_IsLepG7, Em_IsLepG7, Ep_IsLepG7ect, Em_IsLepG7ect, Ep_IsLepG7cc, Em_IsLepG7cc;
  t4->Branch("Ep_IsLepG7",&Ep_IsLepG7,"Ep_IsLepG7/I");
  t4->Branch("Em_IsLepG7",&Em_IsLepG7,"Em_IsLepG7/I");
  t4->Branch("Ep_IsLepG7ect",&Ep_IsLepG7ect,"Ep_IsLepG7ect/I");
  t4->Branch("Em_IsLepG7ect",&Em_IsLepG7ect,"Em_IsLepG7ect/I");
  t4->Branch("Ep_IsLepG7cc",&Ep_IsLepG7cc,"Ep_IsLepG7cc/I");
  t4->Branch("Em_IsLepG7cc",&Em_IsLepG7cc,"Em_IsLepG7cc/I");
  
  
  Int_t nPim, nPip, nProt;
  t4->Branch("nPim",&nPim,"nPim/I");
  t4->Branch("nPip",&nPip,"nPip/I");
  t4->Branch("nProt",&nProt,"nProt/I");
  
  Double_t test_Pull_Chi, test_Pull_Prob;
  
  t4->Branch("test_Pull_Chi",&test_Pull_Chi,"test_Pull_Chi/D");
  t4->Branch("test_Pull_Prob",&test_Pull_Prob,"test_Pull_Prob/D");
  
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
  
  Double_t P_dTOF_fit, Ep_dTOF_fit, Em_dTOF_fit;
  
  t4->Branch("P_dTOF_fit",&P_dTOF_fit,"P_dTOF_fit/D");
  t4->Branch("Ep_dTOF_fit",&Ep_dTOF_fit,"Ep_dTOF_fit/D");
  t4->Branch("Em_dTOF_fit",&Em_dTOF_fit,"Em_dTOF_fit/D");
  
  Double_t P_Timing_fit, Ep_Timing_fit, Em_Timing_fit;
  t4->Branch("P_Timing_fit",&P_Timing_fit,"P_Timing_fit/D");
  t4->Branch("Ep_Timing_fit",&Ep_Timing_fit,"Ep_Timing_fit/D");
  t4->Branch("Em_Timing_fit",&Em_Timing_fit,"Em_Timing_fit/D");
  
  Double_t open_angle_fit, open_angle;
  t4->Branch("open_angle_fit",&open_angle_fit,"open_angle_fit/D");
  t4->Branch("open_angle",&open_angle,"open_angle/D");
  
  Double_t Vx, Vy, Vz;
  t4->Branch("Vx",&Vx,"Vx/D");
  t4->Branch("Vy",&Vy,"Vy/D");
  t4->Branch("Vz",&Vz,"Vz/D");
  
  Double_t P_px_fit, P_py_fit, P_pz_fit, P_Theta_fit, P_Phi_fit, P_Ptot_fit;
  Double_t Ep_px_fit, Ep_py_fit, Ep_pz_fit, Ep_Theta_fit, Ep_Phi_fit, Ep_Ptot_fit;
  Double_t Em_px_fit, Em_py_fit, Em_pz_fit, Em_Theta_fit, Em_Phi_fit, Em_Ptot_fit;
  
  t4->Branch("P_px_fit",&P_px_fit,"P_px_fit/D");
  t4->Branch("P_py_fit",&P_py_fit,"P_py_fit/D");
  t4->Branch("P_pz_fit",&P_pz_fit,"P_pz_fit/D");
  t4->Branch("P_Theta_fit",&P_Theta_fit,"P_Theta_fit/D");
  t4->Branch("P_Phi_fit",&P_Phi_fit,"P_Phi_fit/D");
  t4->Branch("P_Ptot_fit",&P_Ptot_fit,"P_Ptot_fit/D");
  
  t4->Branch("Ep_px_fit",&Ep_px_fit,"Ep_px_fit/D");
  t4->Branch("Ep_py_fit",&Ep_py_fit,"Ep_py_fit/D");
  t4->Branch("Ep_pz_fit",&Ep_pz_fit,"Ep_pz_fit/D");
  t4->Branch("Ep_Theta_fit",&Ep_Theta_fit,"Ep_Theta_fit/D");
  t4->Branch("Ep_Phi_fit",&Ep_Phi_fit,"Ep_Phi_fit/D");
  t4->Branch("Ep_Ptot_fit",&Ep_Ptot_fit,"Ep_Ptot_fit/D");
  
  t4->Branch("Em_px_fit",&Em_px_fit,"Em_px_fit/D");
  t4->Branch("Em_py_fit",&Em_py_fit,"Em_py_fit/D");
  t4->Branch("Em_pz_fit",&Em_pz_fit,"Em_pz_fit/D");
  t4->Branch("Em_Theta_fit",&Em_Theta_fit,"Em_Theta_fit/D");
  t4->Branch("Em_Phi_fit",&Em_Phi_fit,"Em_Phi_fit/D");
  t4->Branch("Em_Ptot_fit",&Em_Ptot_fit,"Em_Ptot_fit/D");
  
  Double_t CM_Theta, CM_Theta_fit;
  t4->Branch("CM_Theta",&CM_Theta,"CM_Theta/D");
  t4->Branch("CM_Theta_fit",&CM_Theta_fit,"CM_Theta_fit/D");
  
  Double_t mandelstam_t, mandelstam_t_fit;
  t4->Branch("mandelstam_t",&mandelstam_t,"mandelstam_t/D");
  t4->Branch("mandelstam_t_fit",&mandelstam_t_fit,"mandelstam_t_fit/D");
  
  
  Int_t trip_pass, ECCC_req, pim_tofpass, pip_tofpass, prot_tofpass, Pass_all;
  t4->Branch("trip_pass",&trip_pass,"trip_pass/I");
  t4->Branch("ECCC_req",&ECCC_req,"ECCC_req/I");
  t4->Branch("pim_tofpass",&pim_tofpass,"pim_tofpass/I");
  t4->Branch("pip_tofpass",&pip_tofpass,"pip_tofpass/I");
  t4->Branch("prot_tofpass",&prot_tofpass,"prot_tofpass/I");
  t4->Branch("Pass_all",&Pass_all,"Pass_all/I");
  
  Int_t Run, Event, TrigBits;
  t4->Branch("Run",&Run,"Run/I");
  t4->Branch("Event",&Event,"Event/I");
  t4->Branch("TrigBits",&TrigBits,"TrigBits/I");
  
  
  //Pions with missing gamma
  Double_t PiGam_Pull_Chi, PiGam_Pull_Prob;
  
  t4->Branch("PiGam_Pull_Chi",&PiGam_Pull_Chi,"PiGam_Pull_Chi/D");
  t4->Branch("PiGam_Pull_Prob",&PiGam_Pull_Prob,"PiGam_Pull_Prob/D");
  
  Double_t mE_PEpEm_PiGam_fit, mm2_PEpEm_PiGam_fit, mm_PEpEm_PiGam_fit, mm_P_PiGam_fit, mm2_P_PiGam_fit, Eg_PiGam_fit, IV_EpEm_PiGam_fit;
  
  t4->Branch("mE_PEpEm_PiGam_fit",&mE_PEpEm_PiGam_fit,"mE_PEpEm_PiGam_fit/D");
  t4->Branch("mm2_PEpEm_PiGam_fit",&mm2_PEpEm_PiGam_fit,"mm2_PEpEm_PiGam_fit/D");
  t4->Branch("mm_PEpEm_PiGam_fit",&mm_PEpEm_PiGam_fit,"mm_PEpEm_PiGam_fit/D");
  t4->Branch("mm_P_PiGam_fit",&mm_P_PiGam_fit,"mm_P_PiGam_fit/D");
  t4->Branch("mm2_P_PiGam_fit",&mm2_P_PiGam_fit,"mm2_P_PiGam_fit/D");
  t4->Branch("IV_EpEm_PiGam_fit",&IV_EpEm_PiGam_fit,"IV_EpEm_PiGam_fit/D");
  
  t4->Branch("Eg_PiGam_fit",&Eg_PiGam_fit,"Eg_PiGam_fit/D");
  
  //Pions with nothing missing
  Double_t PiPull_Chi, PiPull_Prob;
  
  t4->Branch("PiPull_Chi",&PiPull_Chi,"PiPull_Chi/D");
  t4->Branch("PiPull_Prob",&PiPull_Prob,"PiPull_Prob/D");
  
  Double_t mE_PEpEm_Pi_fit, mm2_PEpEm_Pi_fit, mm_PEpEm_Pi_fit, mm_P_Pi_fit, mm2_P_Pi_fit, Eg_Pi_fit, IV_EpEm_Pi_fit;
  
  t4->Branch("mE_PEpEm_Pi_fit",&mE_PEpEm_Pi_fit,"mE_PEpEm_Pi_fit/D");
  t4->Branch("mm2_PEpEm_Pi_fit",&mm2_PEpEm_Pi_fit,"mm2_PEpEm_Pi_fit/D");
  t4->Branch("mm_PEpEm_Pi_fit",&mm_PEpEm_Pi_fit,"mm_PEpEm_Pi_fit/D");
  t4->Branch("mm_P_Pi_fit",&mm_P_Pi_fit,"mm_P_Pi_fit/D");
  t4->Branch("mm2_P_Pi_fit",&mm2_P_Pi_fit,"mm2_P_Pi_fit/D");
  t4->Branch("IV_EpEm_Pi_fit",&IV_EpEm_Pi_fit,"IV_EpEm_Pi_fit/D");
  
  t4->Branch("Eg_Pi_fit",&Eg_Pi_fit,"Eg_Pi_fit/D");
  
  
  //Pions with missing pi0
  Double_t PiPi0Pull_Chi, PiPi0Pull_Prob;
  
  t4->Branch("PiPi0Pull_Chi",&PiPi0Pull_Chi,"PiPi0Pull_Chi/D");
  t4->Branch("PiPi0Pull_Prob",&PiPi0Pull_Prob,"PiPi0Pull_Prob/D");
  
  Double_t mE_PEpEm_PiPi0_fit, mm2_PEpEm_PiPi0_fit, mm_PEpEm_PiPi0_fit, mm_P_PiPi0_fit, mm2_P_PiPi0_fit, Eg_PiPi0_fit, IV_EpEm_PiPi0_fit;
  
  t4->Branch("mE_PEpEm_PiPi0_fit",&mE_PEpEm_PiPi0_fit,"mE_PEpEm_PiPi0_fit/D");
  t4->Branch("mm2_PEpEm_PiPi0_fit",&mm2_PEpEm_PiPi0_fit,"mm2_PEpEm_PiPi0_fit/D");
  t4->Branch("mm_PEpEm_PiPi0_fit",&mm_PEpEm_PiPi0_fit,"mm_PEpEm_PiPi0_fit/D");
  t4->Branch("mm_P_PiPi0_fit",&mm_P_PiPi0_fit,"mm_P_PiPi0_fit/D");
  t4->Branch("mm2_P_PiPi0_fit",&mm2_P_PiPi0_fit,"mm2_P_PiPi0_fit/D");
  t4->Branch("IV_EpEm_PiPi0_fit",&IV_EpEm_PiPi0_fit,"IV_EpEm_PiPi0_fit/D");
  
  t4->Branch("Eg_PiPi0_fit",&Eg_PiPi0_fit,"Eg_PiPi0_fit/D");
  
  Double_t P_px_PiPi0_fit, P_py_PiPi0_fit, P_pz_PiPi0_fit, P_Theta_PiPi0_fit, P_Phi_PiPi0_fit;
  t4->Branch("P_px_PiPi0_fit",&P_px_PiPi0_fit,"P_px_PiPi0_fit/D");
  t4->Branch("P_py_PiPi0_fit",&P_py_PiPi0_fit,"P_py_PiPi0_fit/D");
  t4->Branch("P_pz_PiPi0_fit",&P_pz_PiPi0_fit,"P_pz_PiPi0_fit/D");
  t4->Branch("P_Theta_PiPi0_fit",&P_Theta_PiPi0_fit,"P_Theta_PiPi0_fit/D");
  t4->Branch("P_Phi_PiPi0_fit",&P_Phi_PiPi0_fit,"P_Phi_PiPi0_fit/D");
  
  Double_t Pip_px_PiPi0_fit, Pip_py_PiPi0_fit, Pip_pz_PiPi0_fit, Pip_Theta_PiPi0_fit, Pip_Phi_PiPi0_fit;
  t4->Branch("Pip_px_PiPi0_fit",&Pip_px_PiPi0_fit,"Pip_px_PiPi0_fit/D");
  t4->Branch("Pip_py_PiPi0_fit",&Pip_py_PiPi0_fit,"Pip_py_PiPi0_fit/D");
  t4->Branch("Pip_pz_PiPi0_fit",&Pip_pz_PiPi0_fit,"Pip_pz_PiPi0_fit/D");
  t4->Branch("Pip_Theta_PiPi0_fit",&Pip_Theta_PiPi0_fit,"Pip_Theta_PiPi0_fit/D");
  t4->Branch("Pip_Phi_PiPi0_fit",&Pip_Phi_PiPi0_fit,"Pip_Phi_PiPi0_fit/D");
  
  Double_t Pim_px_PiPi0_fit, Pim_py_PiPi0_fit, Pim_pz_PiPi0_fit, Pim_Theta_PiPi0_fit, Pim_Phi_PiPi0_fit;
  t4->Branch("Pim_px_PiPi0_fit",&Pim_px_PiPi0_fit,"Pim_px_PiPi0_fit/D");
  t4->Branch("Pim_py_PiPi0_fit",&Pim_py_PiPi0_fit,"Pim_py_PiPi0_fit/D");
  t4->Branch("Pim_pz_PiPi0_fit",&Pim_pz_PiPi0_fit,"Pim_pz_PiPi0_fit/D");
  t4->Branch("Pim_Theta_PiPi0_fit",&Pim_Theta_PiPi0_fit,"Pim_Theta_PiPi0_fit/D");
  t4->Branch("Pim_Phi_PiPi0_fit",&Pim_Phi_PiPi0_fit,"Pim_Phi_PiPi0_fit/D");
  
  Double_t P_Ptot_PiPi0_fit, Pip_Ptot_PiPi0_fit, Pim_Ptot_PiPi0_fit;
  t4->Branch("P_Ptot_PiPi0_fit",&P_Ptot_PiPi0_fit,"P_Ptot_PiPi0_fit/D");
  t4->Branch("Pip_Ptot_PiPi0_fit",&Pip_Ptot_PiPi0_fit,"Pip_Ptot_PiPi0_fit/D");
  t4->Branch("Pim_Ptot_PiPi0_fit",&Pim_Ptot_PiPi0_fit,"Pim_Ptot_PiPi0_fit/D");
  
  Double_t P_dTOF_PiPi0_fit, Pip_dTOF_PiPi0_fit, Pim_dTOF_PiPi0_fit;
  t4->Branch("P_dTOF_PiPi0_fit",&P_dTOF_PiPi0_fit,"P_dTOF_PiPi0_fit/D");
  t4->Branch("Pip_dTOF_PiPi0_fit",&Pip_dTOF_PiPi0_fit,"Pip_dTOF_PiPi0_fit/D");
  t4->Branch("Pim_dTOF_PiPi0_fit",&Pim_dTOF_PiPi0_fit,"Pim_dTOF_PiPi0_fit/D");
  
  Double_t P_Timing_PiPi0_fit, Pip_Timing_PiPi0_fit, Pim_Timing_PiPi0_fit;
  t4->Branch("P_Timing_PiPi0_fit",&P_Timing_PiPi0_fit,"P_Timing_PiPi0_fit/D");
  t4->Branch("Pip_Timing_PiPi0_fit",&Pip_Timing_PiPi0_fit,"Pip_Timing_PiPi0_fit/D");
  t4->Branch("Pim_Timing_PiPi0_fit",&Pim_Timing_PiPi0_fit,"Pim_Timing_PiPi0_fit/D");
  
  Double_t Pip_Timing, Pim_Timing;
  t4->Branch("Pip_Timing",&Pip_Timing,"Pip_Timing/D");
  t4->Branch("Pim_Timing",&Pim_Timing,"Pim_Timing/D");
  
  Double_t Pip_dTOF, Pim_dTOF;
  t4->Branch("Pip_dTOF",&Pip_dTOF,"Pip_dTOF/D");
  t4->Branch("Pim_dTOF",&Pim_dTOF,"Pim_dTOF/D");
  
  
  
  //Eta test prob
  Double_t Eta_test_Pull_Chi, Eta_test_Pull_Prob;
  
  t4->Branch("Eta_test_Pull_Chi",&Eta_test_Pull_Chi,"Eta_test_Pull_Chi/D");
  t4->Branch("Eta_test_Pull_Prob",&Eta_test_Pull_Prob,"Eta_test_Pull_Prob/D");
  
  Double_t M_P = 0.938272;   //Proton
  Double_t M_Pi = 0.139570;  //Pion
  Double_t M_PiZ = 0.1349766;  //Pion Zero
  Double_t M_Eta = 0.547853;  //Eta
  Double_t M_Eta_Prime = 0.95778;  //Eta_Prime
  Double_t Melectron = 0.000510999; //Electron
  Double_t c = 29.9792458; //units m/s/e7
  Double_t pi = TMath::Pi();
  Double_t DegToRad = pi/180.0;
  
  string tripdir = "/Volumes/Mac_Storage/Work_Data/G12_NECCESSITIES/TRIPFILES";

  Int_t GOOD = 0;
  Int_t BAD = 0;
  
  for (Int_t j=0;j<=nentries;j++) {//nentries
    chain->GetEntry(j);
    
    if(!(j%95000)) std::cout << "\r done " << j << " out of " << nentries << " ==> " << double(j)*100.0/double(nentries) << "%" << flush;
    if(j== nentries) std::cout << " DONE" << endl;
    
    
    trip_pass = 0;
    ECCC_req = 0;
    pim_tofpass = 0;
    pip_tofpass = 0;
    prot_tofpass = 0;
    Pass_all = 0;
    
    if (NPip ==1 && NPim ==1 && NProt ==1 && abs(MM2PEpEm)<0.075) {
      
      Int_t trip_pass_here = 0;
      
      if (is_good(run, event, tripdir))
      {
        //cout<<"GOOD"<<endl;
        trip_pass = 1;
        trip_pass_here =1;
      }
      else{
        //cout<<"BAD"<<endl;
        trip_pass = -1;
      }
      
      TVector3 pimUVW = clas::g12::g12_ECxyz_2uvw(pim_ECx, pim_ECy, pim_ECz);// EC UVW
      TVector3 pipUVW = clas::g12::g12_ECxyz_2uvw(pip_ECx, pip_ECy, pip_ECz);// EC UVW
      float pim_u = pimUVW.X();
      float pim_v = pimUVW.Y();
      float pim_w = pimUVW.Z();
      float pip_u = pipUVW.X();
      float pip_v = pipUVW.Y();
      float pip_w = pipUVW.Z();
      
      int pip_ecpass = clas::g12::g12_ec_knockout(pip_ECin, pip_ECout, pip_u, pip_v, pip_w, pip_sec);// EC Knockout
      int pim_ecpass = clas::g12::g12_ec_knockout(pim_ECin, pim_ECout, pim_u, pim_v, pim_w, pim_sec);// EC Knockout
      
      Int_t ECCC_req_here = 0;
      
      if ((pip_EC_hit ==1 && pip_ecpass  && pip_CChit_stat ==1) || (pim_EC_hit ==1 && pim_ecpass && pim_CChit_stat ==1) || (pip_EC_hit ==1 && pip_ecpass && pim_CChit_stat ==1) || (pim_EC_hit ==1 && pim_ecpass && pip_CChit_stat ==1)) {
        ECCC_req = 1;
        ECCC_req_here = 1;
      }
      else{
        ECCC_req = -1;
      }
      
      int pim_tof = clas::g12::g12_TOF_knockout(pim_Theta, pim_Phi, pim_sec);// TOF Knockout
      int pip_tof = clas::g12::g12_TOF_knockout(pip_Theta, pip_Phi, pip_sec);// TOF Knockout
      int prot_tof = clas::g12::g12_TOF_knockout(prot_Theta, prot_Phi, prot_sec);// TOF Knockout
      
      if (pim_tof) {
        pim_tofpass = 1;
      }
      else{
        pim_tofpass = -1;
      }
      
      if (pip_tof) {
        pip_tofpass = 1;
      }
      else{
        pip_tofpass = -1;
      }
      
      if (prot_tof) {
        prot_tofpass = 1;
      }
      else{
        prot_tofpass = -1;
      }
      
      if (ECCC_req_here && pim_tof && pip_tof && prot_tof && trip_pass_here)
      {
        //cout<<"GOOD"<<endl;
        GOOD++;
        Good_events++;
        Pass_all = 1;
        
      }
      else{
        //cout<<"BAD"<<endl;
        BAD++;
        Bad_events++;
        Pass_all = -1;
      }
      
      
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
      
      
      
      Ep_EChit = pip_EC_hit;
      Ep_CChit = pip_CC_hit;
      Em_EChit = pim_EC_hit;
      Em_CChit = pim_CC_hit;
      
      Ep_IsLepG7 = pip_IsLepG7;
      Em_IsLepG7 = pim_IsLepG7;
      Ep_IsLepG7ect = pip_IsLepG7ect;
      Em_IsLepG7ect = pim_IsLepG7ect;
      Ep_IsLepG7cc = pip_IsLepG7cc;
      Em_IsLepG7cc = pim_IsLepG7cc;
      
      
      mE_PEpEm = (double)MEPEpEm; mm2_PEpEm = (double)MM2PEpEm;
      mm2_P = (double)MM2P; E_g = (double)egam; IV_EpEm = (double)IVEpEm;
      Vx = (double)vx; Vy = (double)vy; Vz = (double)vz;
      
      
      //if(ECCC_trigger){
      
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
      
      
      TMatrixD covTrack(10,10),covMatrix(10,10), covtestMatrix(10,10), covPiGamMatrix(10,10),covPiMatrix(10,10), covPiPi0Matrix(10,10), covEta_testMatrix(10,10);
      
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
      
      p4[0] = P4pro; ptest[0] = P4pro;
      p4[1] = P4ep;  ptest[1] = P4ep;
      p4[2] = P4em;  ptest[2] = P4em;
      
      vert[0] = V3pro;
      vert[1] = V3ep;
      vert[2] = V3em;
      
      particles[0] = "p";
      particles[1] = "e+";
      particles[2] = "e-";
      
      string experiment = "g12";
      
      covtestMatrix = CorrectCLAS_V(covTrack,particles,ptest,vert,multi,is_mc,experiment);
      Kstream testfit;
      testfit.StringNames(particles);
      testfit.FitInput(e_gamma,ptest,covtestMatrix,m_targ);
      bool include = true; //Include missing particle in constraint no (false).
      
      testfit.Fit(0,settest,include,0.1349);
      //testfit.Fit(0.1349,settest,include,0.782);

      
      test_Pull_Chi = testfit.Chi2();
      test_Pull_Prob = testfit.Prob();
      
      
      covMatrix = CorrectCLAS_V(covTrack,particles,p4,vert,multi,is_mc,experiment);
      
      
      Kstream gamfit;
      gamfit.StringNames(particles);
      gamfit.FitInput(e_gamma,p4,covMatrix,m_targ);
      //gamfit.Fit("pi0");
      gamfit.Fit("gamma");

      
      Pull_Chi = gamfit.Chi2();
      Pull_Prob = gamfit.Prob();
      
      
      Eg_fit = gamfit.FitPhotonEnergy();
      
      for(int i = 0; i < 3; i++) p4[i] = gamfit.FitP4(i);
      //for(int i = 0; i < 2; i++) p4[i] = pifit.FitP4(i);
      
      P4pho_fit.SetPxPyPzE(0.0,0.0,Eg_fit,Eg_fit);
      P4pro_fit = p4[0];
      P4pip_fit = p4[1];
      P4pim_fit = p4[2];
      
      TLorentzVector MM_PVec = (P4pho_fit + P4target) - P4pro_fit;
      TLorentzVector MM_PEpEmVec = (P4pho_fit + P4target) - (P4pro_fit + P4pip_fit + P4pim_fit);
      TLorentzVector IV_EpEmVec = P4pip_fit + P4pim_fit;
      open_angle_fit = P4pip_fit.Vect().Angle(P4pim_fit.Vect());
      IV_EpEm_fit_P = sqrt(pow(IV_EpEmVec.Px(),2) + pow(IV_EpEmVec.Py(),2) + pow(IV_EpEmVec.Pz(),2));
      
      CM_Theta_fit = OneBoostComparison(MM_PVec, (P4pho_fit + P4target));
      mandelstam_t_fit = (P4target - P4pro_fit).M2();

      
      mE_PEpEm_fit = MM_PEpEmVec.E();
      mm2_PEpEm_fit = MM_PEpEmVec.M2();
      mm_PEpEm_fit = MM_PEpEmVec.M();
      
      mm_P_fit = MM_PVec.M();
      mm2_P_fit = MM_PVec.M2();
      IV_EpEm_fit = IV_EpEmVec.M();
      
      
      
      P_px_fit = P4pro_fit.Px();
      P_py_fit = P4pro_fit.Py();
      P_pz_fit = P4pro_fit.Pz();
      P_Theta_fit = P4pro_fit.Theta()*180./TMath::Pi();
      P_Phi_fit = P4pro_fit.Phi()*180./TMath::Pi();
      
      Ep_px_fit = P4pip_fit.Px();
      Ep_py_fit = P4pip_fit.Py();
      Ep_pz_fit = P4pip_fit.Pz();
      Ep_Theta_fit = P4pip_fit.Theta()*180./TMath::Pi();
      Ep_Phi_fit = P4pip_fit.Phi()*180./TMath::Pi();
      
      Em_px_fit = P4pim_fit.Px();
      Em_py_fit = P4pim_fit.Py();
      Em_pz_fit = P4pim_fit.Pz();
      Em_Theta_fit = P4pim_fit.Theta()*180./TMath::Pi();
      Em_Phi_fit = P4pim_fit.Phi()*180./TMath::Pi();
      
      P_Ptot_fit = sqrt(pow(P_px_fit,2) + pow(P_py_fit,2) + pow(P_pz_fit,2));
      Ep_Ptot_fit = sqrt(pow(Ep_px_fit,2) + pow(Ep_py_fit,2) + pow(Ep_pz_fit,2));
      Em_Ptot_fit = sqrt(pow(Em_px_fit,2) + pow(Em_py_fit,2) + pow(Em_pz_fit,2));
      
      P_dTOF_fit = ((prot_SC_time - prot_ST_time) - (prot_STSCpl/P_Ptot_fit)*sqrt(pow(M_P,2) + pow(P_Ptot_fit,2))/c);
      Ep_dTOF_fit = ((pip_SC_time - pip_ST_time) - (pip_STSCpl/Ep_Ptot_fit)*sqrt(pow(Melectron,2) + pow(Ep_Ptot_fit,2))/c);
      Em_dTOF_fit = ((pim_SC_time - pim_ST_time) - (pim_STSCpl/Em_Ptot_fit)*sqrt(pow(Melectron,2) + pow(Em_Ptot_fit,2))/c);
      
      
      Double_t P_beta_fit = P_Ptot_fit/sqrt(P_Ptot_fit*P_Ptot_fit + M_P * M_P);
      Double_t Ep_beta_fit = Ep_Ptot_fit/sqrt(Ep_Ptot_fit*Ep_Ptot_fit + Melectron * Melectron);
      Double_t Em_beta_fit = Em_Ptot_fit/sqrt(Em_Ptot_fit*Em_Ptot_fit + Melectron * Melectron);
      
      
      P_Timing_fit = tpho + prot_Tprop - (prot_SC_time -  prot_scLen/(c*P_beta_fit));
      Ep_Timing_fit = tpho + pip_Tprop - (pip_SC_time - pip_scLen/(c*Ep_beta_fit));
      Em_Timing_fit = tpho + pim_Tprop - (pim_SC_time - pim_scLen/(c*Em_beta_fit));
      
      
      //PIONS
      std::vector<string> piparticles(num_parts);
      
      std::vector<TLorentzVector> pPiGam(num_parts);
      std::vector<TLorentzVector> pPi(num_parts);
      std::vector<TLorentzVector> pPiPi0(num_parts);
      
      pPiGam[0] = P4pro;   pPi[0] = P4pro;    pPiPi0[0] = P4pro;
      pPiGam[1] = P4Pip;    pPi[1] = P4Pip;     pPiPi0[1] = P4Pip;
      pPiGam[2] = P4Pim;    pPi[2] = P4Pim;     pPiPi0[2] = P4Pim;
      
      piparticles[0] = "p";
      piparticles[1] = "pi+";
      piparticles[2] = "pi-";
      
      
      
      //PiGam
      covPiGamMatrix = CorrectCLAS_V(covTrack,piparticles,pPiGam,vert,multi,is_mc,experiment);
      
      Kstream PiGamfit;
      PiGamfit.StringNames(piparticles);
      PiGamfit.FitInput(e_gamma,pPiGam,covPiGamMatrix,m_targ);
      PiGamfit.Fit("gamma");
      
      PiGam_Pull_Chi = PiGamfit.Chi2();
      PiGam_Pull_Prob = PiGamfit.Prob();
      
      Eg_PiGam_fit = PiGamfit.FitPhotonEnergy();
      
      for(int i = 0; i < 3; i++) pPiGam[i] = PiGamfit.FitP4(i);
      
      P4pho_PiGam_fit.SetPxPyPzE(0.0,0.0,Eg_PiGam_fit,Eg_PiGam_fit);
      P4pro_PiGam_fit = pPiGam[0];
      P4ep_PiGam_fit = pPiGam[1];
      P4em_PiGam_fit = pPiGam[2];
      
      TLorentzVector MM_PVec_PiGam = (P4pho_PiGam_fit + P4target) - P4pro_PiGam_fit;
      TLorentzVector MM_PEpEmVec_PiGam = (P4pho_PiGam_fit + P4target) - (P4pro_PiGam_fit + P4ep_PiGam_fit + P4em_PiGam_fit);
      TLorentzVector IV_EpEmVec_PiGam = P4ep_PiGam_fit + P4em_PiGam_fit;
      
      
      mE_PEpEm_PiGam_fit = MM_PEpEmVec_PiGam.E();
      mm2_PEpEm_PiGam_fit = MM_PEpEmVec_PiGam.M2();
      mm_PEpEm_PiGam_fit = MM_PEpEmVec_PiGam.M();
      
      mm_P_PiGam_fit = MM_PVec_PiGam.M();
      mm2_P_PiGam_fit = MM_PVec_PiGam.M2();
      IV_EpEm_PiGam_fit = IV_EpEmVec_PiGam.M();
      
      //Pi
      covPiMatrix = CorrectCLAS_V(covTrack,piparticles,pPi,vert,multi,is_mc,experiment);
      
      Kstream Pifit;
      Pifit.StringNames(piparticles);
      Pifit.FitInput(e_gamma,pPi,covPiMatrix,m_targ);
      Pifit.Fit();
      
      PiPull_Chi = Pifit.Chi2();
      PiPull_Prob = Pifit.Prob();
      
      Eg_Pi_fit = Pifit.FitPhotonEnergy();
      
      for(int i = 0; i < 3; i++) pPi[i] = Pifit.FitP4(i);
      
      P4pho_Pi_fit.SetPxPyPzE(0.0,0.0,Eg_Pi_fit,Eg_Pi_fit);
      P4pro_Pi_fit = pPi[0];
      P4ep_Pi_fit = pPi[1];
      P4em_Pi_fit = pPi[2];
      
      TLorentzVector MM_PVec_Pi = (P4pho_Pi_fit + P4target) - P4pro_Pi_fit;
      TLorentzVector MM_PEpEmVec_Pi = (P4pho_Pi_fit + P4target) - (P4pro_Pi_fit + P4ep_Pi_fit + P4em_Pi_fit);
      TLorentzVector IV_EpEmVec_Pi = P4ep_Pi_fit + P4em_Pi_fit;
      
      
      mE_PEpEm_Pi_fit = MM_PEpEmVec_Pi.E();
      mm2_PEpEm_Pi_fit = MM_PEpEmVec_Pi.M2();
      mm_PEpEm_Pi_fit = MM_PEpEmVec_Pi.M();
      
      mm_P_Pi_fit = MM_PVec_Pi.M();
      mm2_P_Pi_fit = MM_PVec_Pi.M2();
      IV_EpEm_Pi_fit = IV_EpEmVec_Pi.M();
      
      //PiPi0
      covPiPi0Matrix = CorrectCLAS_V(covTrack,piparticles,pPiPi0,vert,multi,is_mc,experiment);
      
      Kstream PiPi0fit;
      PiPi0fit.StringNames(piparticles);
      PiPi0fit.FitInput(e_gamma,pPiPi0,covPiPi0Matrix,m_targ);
      PiPi0fit.Fit("pi0");
      
      PiPi0Pull_Chi = PiPi0fit.Chi2();
      PiPi0Pull_Prob = PiPi0fit.Prob();
      
      Eg_PiPi0_fit = PiPi0fit.FitPhotonEnergy();
      
      for(int i = 0; i < 3; i++) pPiPi0[i] = PiPi0fit.FitP4(i);
      
      P4pho_PiPi0_fit.SetPxPyPzE(0.0,0.0,Eg_PiPi0_fit,Eg_PiPi0_fit);
      P4pro_PiPi0_fit = pPiPi0[0];
      P4ep_PiPi0_fit = pPiPi0[1];
      P4em_PiPi0_fit = pPiPi0[2];
      
      TLorentzVector MM_PVec_PiPi0 = (P4pho_PiPi0_fit + P4target) - P4pro_PiPi0_fit;
      TLorentzVector MM_PEpEmVec_PiPi0 = (P4pho_PiPi0_fit + P4target) - (P4pro_PiPi0_fit + P4ep_PiPi0_fit + P4em_PiPi0_fit);
      TLorentzVector IV_EpEmVec_PiPi0 = P4ep_PiPi0_fit + P4em_PiPi0_fit;
      
      
      mE_PEpEm_PiPi0_fit = MM_PEpEmVec_PiPi0.E();
      mm2_PEpEm_PiPi0_fit = MM_PEpEmVec_PiPi0.M2();
      mm_PEpEm_PiPi0_fit = MM_PEpEmVec_PiPi0.M();
      
      mm_P_PiPi0_fit = MM_PVec_PiPi0.M();
      mm2_P_PiPi0_fit = MM_PVec_PiPi0.M2();
      IV_EpEm_PiPi0_fit= IV_EpEmVec_PiPi0.M();
      

      
      
      
      
      
      
      
      ///////
      std::vector<TLorentzVector> Etatest(num_parts);
      std::vector<bool> Eta_settest(num_parts);
      
      Eta_settest[0] = false;
      Eta_settest[1] = true;
      Eta_settest[2] = true;
      
      Etatest[0] = P4pro;
      Etatest[1] = P4Pip;
      Etatest[2] = P4Pim;
      
      
      covEta_testMatrix = CorrectCLAS_V(covTrack,piparticles,Etatest,vert,multi,is_mc,experiment);
      Kstream Eta_testfit;
      Eta_testfit.StringNames(particles);
      Eta_testfit.FitInput(e_gamma,Etatest,covEta_testMatrix,m_targ);
      bool Eta_include = true; //Include missing particle in constraint no (false).
      
      Eta_testfit.Fit(0,Eta_settest,Eta_include,0.54786);
      
      Eta_test_Pull_Chi = Eta_testfit.Chi2();
      Eta_test_Pull_Prob = Eta_testfit.Prob();
      
      
      //////
      
      P_px_PiPi0_fit = P4pro_PiPi0_fit.Px();
      P_py_PiPi0_fit = P4pro_PiPi0_fit.Py();
      P_pz_PiPi0_fit = P4pro_PiPi0_fit.Pz();
      P_Theta_PiPi0_fit = P4pro_PiPi0_fit.Theta()*180./TMath::Pi();
      P_Phi_PiPi0_fit = P4pro_PiPi0_fit.Phi()*180./TMath::Pi();
      
      Pip_px_PiPi0_fit = P4ep_PiPi0_fit.Px();
      Pip_py_PiPi0_fit = P4ep_PiPi0_fit.Py();
      Pip_pz_PiPi0_fit = P4ep_PiPi0_fit.Pz();
      Pip_Theta_PiPi0_fit = P4ep_PiPi0_fit.Theta()*180./TMath::Pi();
      Pip_Phi_PiPi0_fit = P4ep_PiPi0_fit.Phi()*180./TMath::Pi();
      
      Pim_px_PiPi0_fit = P4em_PiPi0_fit.Px();
      Pim_py_PiPi0_fit = P4em_PiPi0_fit.Py();
      Pim_pz_PiPi0_fit = P4em_PiPi0_fit.Pz();
      Pim_Theta_PiPi0_fit = P4em_PiPi0_fit.Theta()*180./TMath::Pi();
      Pim_Phi_PiPi0_fit = P4em_PiPi0_fit.Phi()*180./TMath::Pi();
      
      P_Ptot_PiPi0_fit = sqrt(pow(P_px_PiPi0_fit,2) + pow(P_py_PiPi0_fit,2) + pow(P_pz_PiPi0_fit,2));
      Pip_Ptot_PiPi0_fit = sqrt(pow(Pip_px_PiPi0_fit,2) + pow(Pip_py_PiPi0_fit,2) + pow(Pip_pz_PiPi0_fit,2));
      Pim_Ptot_PiPi0_fit = sqrt(pow(Pim_px_PiPi0_fit,2) + pow(Pim_py_PiPi0_fit,2) + pow(Pim_pz_PiPi0_fit,2));
      
      P_dTOF_PiPi0_fit = ((prot_SC_time - prot_ST_time) - (prot_STSCpl/P_Ptot_PiPi0_fit)*sqrt(pow(M_P,2) + pow(P_Ptot_PiPi0_fit,2))/c);
      Pip_dTOF_PiPi0_fit = ((pip_SC_time - pip_ST_time) - (pip_STSCpl/Pip_Ptot_PiPi0_fit)*sqrt(pow(Melectron,2) + pow(Pip_Ptot_PiPi0_fit,2))/c);
      Pim_dTOF_PiPi0_fit = ((pim_SC_time - pim_ST_time) - (pim_STSCpl/Pim_Ptot_PiPi0_fit)*sqrt(pow(Melectron,2) + pow(Pim_Ptot_PiPi0_fit,2))/c);
      
      
      Double_t P_beta_PiPi0_fit = P_Ptot_PiPi0_fit/sqrt(P_Ptot_PiPi0_fit*P_Ptot_PiPi0_fit + M_P * M_P);
      Double_t Pip_beta_PiPi0_fit = Pip_Ptot_PiPi0_fit/sqrt(Pip_Ptot_PiPi0_fit*Pip_Ptot_PiPi0_fit + Melectron * Melectron);
      Double_t Pim_beta_PiPi0_fit = Pim_Ptot_PiPi0_fit/sqrt(Pim_Ptot_PiPi0_fit*Pim_Ptot_PiPi0_fit + Melectron * Melectron);
      
      
      P_Timing_PiPi0_fit = tpho + prot_Tprop - (prot_SC_time -  prot_scLen/(c*P_beta_PiPi0_fit));
      Pip_Timing_PiPi0_fit = tpho + pip_Tprop - (pip_SC_time - pip_scLen/(c*Pip_beta_PiPi0_fit));
      Pim_Timing_PiPi0_fit = tpho + pim_Tprop - (pim_SC_time - pim_scLen/(c*Pim_beta_PiPi0_fit));
      
      Double_t Pip_beta = pip_P/sqrt(pip_P*pip_P + M_Pi * M_Pi);
      Double_t Pim_beta = pim_P/sqrt(pim_P*pim_P + M_Pi * M_Pi);
      
      Pip_Timing = tpho + pip_Tprop - (pip_SC_time - pip_scLen/(c*Pip_beta));
      Pim_Timing = tpho + pim_Tprop - (pim_SC_time - pim_scLen/(c*Pim_beta));
      
      Pip_dTOF = ((pip_SC_time - pip_ST_time) - (pip_STSCpl/pip_P)*sqrt(pow(M_Pi,2) + pow(pip_P,2))/c);
      Pim_dTOF = ((pim_SC_time - pim_ST_time) - (pim_STSCpl/pim_P)*sqrt(pow(M_Pi,2) + pow(pim_P,2))/c);

      
      t4->Fill();
    }
    //t4->Fill();
    
  }
  t4->Write();
  cout<<GOOD<<"  "<<BAD<<endl;
  
  cout<<"######-------------------------------------#####"<<endl;
  cout<<Good_events<<" are GOOD "<<Bad_events<<" are BAD "<<endl;
  cout<<"######-------------------------------------#####"<<endl;
  
  outFile.Write(); // write to the output file
  outFile.Close(); // close the output file
  
}//end of main



void Check_File(TString File_in){ //, TString *File_out
  
  
  
  TString in_dir = "/Volumes/DATA/CLAS_LEPTON_G12/RAW/";
  TString out_dir = "/Volumes/DATA/CLAS_LEPTON_G12/FITTED/";
  TString f_out = "SkimmedFiducial_";
  TString IN = in_dir + File_in;
  TString OUT = out_dir + f_out + File_in;
  
  //cout<<IN<<endl;
  //cout<<OUT<<endl;
  
  //  (*File_out) = OUT;
  //  (*File_in) = IN;
  
  TFile *fill = TFile::Open(OUT);
  if (!fill) {
    //cout<<"NOT HERE"<<endl;
    Fiducial_Clone_Tree(IN, OUT);
  }
  delete fill;
  //else{cout<<"ALREADY PROCESSED"<<endl;}
  
  
}

//void Copy_PPipPim_tree(){
int main(){
  int k=0;
  //ifstream in1("/Volumes/Seagate_Storage_2/WORK_DATA/ALL_POSSIBLE_LEPTONS/HADDED/totfile1");
  ifstream in1("/Volumes/DATA/CLAS_LEPTON_G12/RAW/totfile");

  string amount;
  while(in1 >> amount ){
    k++;
  }
  in1.close();
  cout <<"You have " << k <<" files"<<endl;
  
  //ifstream in("/Volumes/Seagate_Storage_2/WORK_DATA/ALL_POSSIBLE_LEPTONS/HADDED/totfile1");
  ifstream in("//Volumes/DATA/CLAS_LEPTON_G12/RAW/totfile");

  for(int i = 1; i<=k; i++){
    TString processed_file_;
    TString original_file;
    in >> original_file;
    Check_File(original_file); //, &processed_file_
  }
}

