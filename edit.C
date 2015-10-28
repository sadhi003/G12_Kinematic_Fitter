//For pi+pi-(gamma) fitting
TMatrixD covMatrix(10,10);

covMatrix = CorrectCLAS_V(covTrack,particles,p4,vert,multi,is_mc,experiment);
Kstream gamfit_diPion;
gamfit_diPion.StringNames(particles);
gamfit_diPion.FitInput(e_gamma,p4,covMatrix,m_targ);
gamfit_diPion.Fit("gamma");

Pull_Chi_gamfit_diPion = gamfit_diPion.Chi2();
Pull_Prob_gamfit_diPion = gamfit_diPion.Prob();

Eg_gamfit_diPion = gamfit_diPion.FitPhotonEnergy();

for(int i = 0; i < 3; i++) p4[i] = gamfit_diPion.FitP4(i);

P4pho_gamfit_diPion.SetPxPyPzE(0.0,0.0,Eg_gamfit_diPion,Eg_gamfit_diPion);
P4pro_gamfit_diPion = p4[0];
P4ep_gamfit_diPion = p4[1];
P4em_gamfit_diPion = p4[2];

TLorentzVector MM_PVec = (P4pho_gamfit_diPion + P4target) - P4pro_gamfit_diPion;
TLorentzVector MM_PEpEmVec = (P4pho_gamfit_diPion + P4target) - (P4pro_gamfit_diPion + P4ep_gamfit_diPion + P4em_gamfit_diPion);
TLorentzVector IV_EpEmVec = P4ep_gamfit_diPion + P4em_gamfit_diPion;
open_angle_gamfit_diPion = P4ep_gamfit_diPion.Vect().Angle(P4em_gamfit_diPion.Vect());
IV_EpEm_gamfit_diPion_P = sqrt(pow(IV_EpEmVec.Px(),2) + pow(IV_EpEmVec.Py(),2) + pow(IV_EpEmVec.Pz(),2));

CM_Theta_gamfit_diPion = OneBoostComparison(MM_PVec, (P4pho_gamfit_diPion + P4target));
mandelstam_t_gamfit_diPion = (P4target - P4pro_gamfit_diPion).M2();


mE_PEpEm_gamfit_diPion = MM_PEpEmVec.E();
mm2_PEpEm_gamfit_diPion = MM_PEpEmVec.M2();
mm_PEpEm_gamfit_diPion = MM_PEpEmVec.M();

mm_P_gamfit_diPion = MM_PVec.M();
mm2_P_gamfit_diPion = MM_PVec.M2();
IV_EpEm_gamfit_diPion = IV_EpEmVec.M();

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
TMatrixD covMatrixpi0(10,10);

covMatrixpi0 = CorrectCLAS_V(covTrack,particles,p4pi0,vert,multi,is_mc,experiment);
Kstream pi0fit_diPion;
pi0fit_diPion.StringNames(particles);
pi0fit_diPion.FitInput(e_gamma,p4pi0,covMatrix,m_targ);
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
TMatrixD covMatrixnothing(10,10);

covMatrixnothing = CorrectCLAS_V(covTrack,particles,p4nothing,vert,multi,is_mc,experiment);
Kstream nothingfit_diPion;
nothingfit_diPion.StringNames(particles);
nothingfit_diPion.FitInput(e_gamma,p4nothing,covMatrix,m_targ);
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