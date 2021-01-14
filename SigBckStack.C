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
    
// root file
  TFile *f_sig20 = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_h125_4e_zd20.root");
  TFile *f_sig40 = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_h125_4e_zd40.root");
  TFile *f_sig60 = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_h125_4e_zd60.root");
  TFile *f_bkg_zz4lep = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_Background_ZZ4LEP500k.root");
  TFile *f_bkg_pph4l = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_Background_pph4l500k.root");
    
 // get branch in tree
  TH1F *h_sig20 = (TH1F*)f_sig20->Get("mZb_4e")->Clone("h_sig20");
  TH1F *h_sig40 = (TH1F*)f_sig40->Get("mZb_4e")->Clone("h_sig40");
  TH1F *h_sig60 = (TH1F*)f_sig60->Get("mZb_4e")->Clone("h_sig60");
  TH1F *h_bkg_zz4lep = (TH1F*)f_bkg_zz4lep->Get("mZb_4e")->Clone("h_bkg_zz4lep");
  TH1F *h_bkg_pph4l = (TH1F*)f_bkg_pph4l->Get("mZb_4e")->Clone("h_bkg_pph4l");
    
 // getting num of event in bin of interest  
    double nsig20_Before = h_sig20->Integral(h_sig20->FindFixBin(15.), h_sig20->FindFixBin(120.));
    cout << "signal_20_4e nevent Before Normalization is = " << nsig20_Before << ".\n";
    
    double nsig40_Before = h_sig40->Integral(h_sig40->FindFixBin(15.), h_sig40->FindFixBin(120.));
    cout << "signal_40_4e nevent Before Normalization is = " << nsig40_Before << ".\n";
    
    double nsig60_Before = h_sig60->Integral(h_sig60->FindFixBin(15.), h_sig60->FindFixBin(120.));
    cout << "signal_60_4e nevent Before Normalization is = " << nsig60_Before << ".\n";
    
    double nbkg_zz4lepBefore = h_bkg_zz4lep->Integral(h_bkg_zz4lep->FindFixBin(15.), h_bkg_zz4lep->FindFixBin(120.));
    cout << "zz4lep nevent Before Normalization is = " << nbkg_zz4lepBefore << ".\n";
    
    double nbkg_pph4lBefore = h_bkg_pph4l->Integral(h_bkg_pph4l->FindFixBin(15.), h_bkg_pph4l->FindFixBin(120.));
    cout << "pph4lep nevent Before Normalization is = " << nbkg_pph4lBefore << ".\n";
   

  //250fb-1 = 250000pb
  float lumi = 250000;
  float lumi = 250000;
  double nevents = 25000;
  double backevents = 500000;
  double BR_SMHiggs = 1.25e-4;
  double xs = 48.58;
  double ATLASUP = 0.01;
  //based on ATLAS CL upper limit = 10e-3
    
  // normalize using scale to SM Higgs crossection & BR (total events 25k)
    double h = xs*ATLASUP*BR_SMHiggs*lumi/nevents;
    cout << "h =" << h << ".\n";
    
    double back = xs*ATLASUP*BR_SMHiggs*lumi/backevents;
    cout << "back =" << back << ".\n";
  
    h_sig20->Scale(h);
    h_sig40->Scale(h);
    h_sig60->Scale(h);
    h_bkg_pph4l->Scale(back);
    h_bkg_zz4lep->Scale(back);
    
    /*h_sig20->Scale(h/h_sig20->Integral(h_sig20->FindFixBin(15.), h_sig20->FindFixBin(120.)));
    h_sig40->Scale(h/h_sig40->Integral(h_sig40->FindFixBin(15.), h_sig40->FindFixBin(120.)));
    h_sig60->Scale(h/h_sig60->Integral(h_sig60->FindFixBin(15.), h_sig60->FindFixBin(120.)));
    h_bkg_zz4lep->Scale(back/h_bkg_zz4lep->Integral(h_bkg_zz4lep->FindFixBin(15.), h_bkg_zz4lep->FindFixBin(120.)));
    h_bkg_pph4l->Scale(back/h_bkg_pph4l->Integral(h_bkg_pph4l->FindFixBin(15.), h_bkg_pph4l->FindFixBin(120.)));*/
    
     double nsig20_AFTER = h_sig20->Integral(h_sig20->FindFixBin(15.), h_sig20->FindFixBin(120.));
     cout << "signal_20_4e nevent AFTER Normalization is = " << nsig20_AFTER << ".\n";
     
     double nsig40_AFTER = h_sig40->Integral(h_sig40->FindFixBin(15.), h_sig40->FindFixBin(120.));
     cout << "signal_40_4e nevent AFTER Normalization is = " << nsig40_AFTER << ".\n";
     
     double nsig60_AFTER = h_sig60->Integral(h_sig60->FindFixBin(15.), h_sig60->FindFixBin(120.));
     cout << "signal_60_4e nevent AFTER Normalization is = " << nsig60_AFTER << ".\n";
     
     double nbkg_zz4lepAFTER = h_bkg_zz4lep->Integral(h_bkg_zz4lep->FindFixBin(15.), h_bkg_zz4lep->FindFixBin(120.));
     cout << "zz4lep nevent AFTER Normalization is = " << nbkg_zz4lepAFTER << ".\n";
     
     double nbkg_pph4lAFTER = h_bkg_pph4l->Integral(h_bkg_pph4l->FindFixBin(15.), h_bkg_pph4l->FindFixBin(120.));
     cout << "pph4lep nevent AFTER Normalization is = " << nbkg_pph4lAFTER << ".\n";
    
  //stack normalized samples
  TCanvas *c = new TCanvas("c","c");
  h_sig20->SetFillColor(kMagenta-9);
  h_sig40->SetFillColor(kCyan-9);
  h_sig60->SetFillColor(kYellow-9);
  h_bkg_zz4lep->SetFillColor(kBlue+2);
  h_bkg_pph4l->SetFillColor(kRed+1);

    
 THStack *hs = new THStack("hs", "z_{b} mass");
  hs->Add(h_sig20);
  hs->Add(h_sig40);
  hs->Add(h_sig60);
  hs->Add(h_bkg_zz4lep);
  hs->Add(h_bkg_pph4l);
  hs->Draw("HIST");
  hs->GetXaxis()->SetTitle("mass of z_{b} (4e channel)");
  hs->GetYaxis()->SetTitle("Number of Events");
  hs->GetXaxis()->SetRangeUser(12,80);
  
  TLegend *leg = new TLegend(0.9,0.7,0.7,0.9);
  //(plg hujung x dekat 0,lebar dekat xaxis dr 0,lebar dekat yaxis dr 0,atasbawah on y axis)
  // dia count dr 0-1
  leg->SetHeader("Sample");
  leg->AddEntry("h_sig20","mZ_{b} = 20GeV","f");
  leg->AddEntry("h_sig40","mZ_{b} = 40GeV","f");
  leg->AddEntry("h_sig60","mZ_{b} = 60GeV","f");
  leg->AddEntry("h_bkg_zz4lep","zz4l","f");
  leg->AddEntry("h_bkg_pph4l","pph4l","f");
  leg->Draw();
    
  c->SaveAs("mZb_4mu_zd60.png");
}
