{
  gStyle->SetOptFit(11);
  //gStyle->SetOptFit(0);
  gStyle->SetStatW(0.21); gStyle->SetStatH(0.21);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleX(.25);
  //gStyle->SetLabelSize(0.095, "Y");
  //gStyle->SetLabelSize(0.095, "X");
  
  TChain *chain1 = new TChain("LepTree");
  chain1->Add("Lep_KinFit.root");
  
  //TFile* myfile = new TFile("Lep_KinFit.root");
  TCanvas *c1=new TCanvas("Pulls","Pulls",1400,1000);
  c1->SetFillColor(19);c1->Divide(3,4,0,0); c1->SetFillStyle(4000);
  //gStyle->SetStatFontSize(0.07);
  TH1F *h1 = new TH1F("h1","",10,-5,5);
  TH1F *h2 = new TH1F("h2","",10,-5,5);
  
    TF1 *FitFunc = new TF1("FitFunc","gaus",-3.5,3.5);
    FitFunc->SetFillColor(2);
    FitFunc->SetFillStyle(0);
    FitFunc->SetMarkerStyle(21);
    FitFunc->SetMarkerSize(0.3);
    FitFunc->SetLineWidth(3);
  
  c1->cd(1);chain1->Draw("Pull_One>>Pull_1(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(2);chain1->Draw("Pull_Two>>Pull_2(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(3);chain1->Draw("Pull_Three>>Pull_3(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(4);chain1->Draw("Pull_Four>>Pull_4(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(5);chain1->Draw("Pull_Five>>Pull_5(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(6);chain1->Draw("Pull_Six>>Pull_6(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(7);chain1->Draw("Pull_Seven>>Pull_7(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(8);chain1->Draw("Pull_Eight>>Pull_8(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(9);chain1->Draw("Pull_Nine>>Pull_9(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(10);h1->Draw();
  c1->cd(11);chain1->Draw("Pull_Zero>>Pull_0(100,-5,5)","Pull_Prob_nothingfit_diLep>0.01");
  c1->cd(12);h2->Draw();

  Pull_1->SetTitle("Proton p-pull");
  Pull_2->SetTitle("Proton #lambda-pull");
  Pull_3->SetTitle("Proton #phi-pull");
  
  Pull_4->SetTitle("e^{+} p-pull");
  Pull_5->SetTitle("e^{+} #lambda-pull");
  Pull_6->SetTitle("e^{+} #phi-pull");
  
  Pull_7->SetTitle("e^{-} p-pull");
  Pull_8->SetTitle("e^{-} #lambda-pull");
  Pull_9->SetTitle("e^{-} #phi-pull");
  
  Pull_0->SetTitle("#gamma E-pull");
  

  c1->cd(1);Pull_1->Fit("FitFunc","R");
  c1->cd(2);Pull_2->Fit("FitFunc","R");
  c1->cd(3);Pull_3->Fit("FitFunc","R");
  c1->cd(4);Pull_4->Fit("FitFunc","R");
  c1->cd(5);Pull_5->Fit("FitFunc","R");
  c1->cd(6);Pull_6->Fit("FitFunc","R");
  c1->cd(7);Pull_7->Fit("FitFunc","R");
  c1->cd(8);Pull_8->Fit("FitFunc","R");
  c1->cd(9);Pull_9->Fit("FitFunc","R");
  c1->cd(11);Pull_0->Fit("FitFunc","R");
//
//  //c1->cd(12);Prob->Draw("");
  c1->Print("Wrong_Lep_Pulls_fix_03_11_2015.pdf");
}
