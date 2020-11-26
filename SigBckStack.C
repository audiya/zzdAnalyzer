//==============================//
//     Au Diya 24 nov 2020      //
// calc cross section and stack //
//==============================//
// for Root histogram
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "THStack.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif

void SigBckStack(){
    
  gStyle->SetOptStat(0);


  TFile *f_sig = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_h125_4mu_zd60.root");
  TFile *f_bkg_zz4lep = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_Background_zz4lep.root");
  TFile *f_bkg_pph4l = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_Background_pph4l.root");
  
  TH1F *h_sig = (TH1F*)f_sig->Get("mZb_4mu")->Clone("h_sig");
  TH1F *h_bkg_zz4lep = (TH1F*)f_bkg_zz4lep->Get("mZb_4mu")->Clone("h_bkg_zz4lep");
  TH1F *h_bkg_pph4l = (TH1F*)f_bkg_pph4l->Get("mZb_4mu")->Clone("h_bkg_pph4l");

  float lumi = 250;

  float ntot_sig = 25000;
  float ntot_bkg_zz4lep = 25000;
  float ntot_bkg_pph4l = 25000;
    
  h_sig->Scale(ntot_sig/lumi);
    cout << "Signal cross section " << ntot_sig/lumi << endl;
  h_bkg_zz4lep->Scale(ntot_bkg_zz4lep/lumi);
    cout << "Bkg_zz4lep cross section " << ntot_bkg_zz4lep/lumi << endl;
  h_bkg_pph4l->Scale(ntot_bkg_pph4l/lumi);
  cout << "Bkg_pph4l cross section " << ntot_bkg_pph4l/lumi << endl;
    
  TCanvas *c = new TCanvas("c","c");
  
  h_sig->SetFillColor(kRed);
  h_bkg_zz4lep->SetFillColor(kGreen);
  h_bkg_pph4l->SetFillColor(kBlue);
    
  THStack *hs = new THStack("hs", "zd mass");
  hs->Add(h_bkg_zz4lep);
  hs->Add(h_bkg_pph4l);
  hs->Add(h_sig);
  hs->Draw("HIST");
  hs->GetXaxis()->SetTitle("mass of zd (4mu channel)");
  hs->GetYaxis()->SetTitle("Number of Events");
  
  printf("GlobalMaximum: %f\n",hs->GetMaximum());
  
  TH1 *h = ((TH1*)(hs->GetStack()->Last()));
  
  h->GetXaxis()->SetRangeUser(0,100000);
  printf("LocalMaximum: %f\n",h->GetMaximum());
    
  c->SaveAs("mZb_4mu_zd60.png");
}
