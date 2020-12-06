//==============================//
//     Au Diya 24 nov 2020      //
// calc cross section and stack //
//==============================//
// for Root histogram
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TAttMarker.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFile.h"
#include "math.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif
using namespace std;

void significance(){
    
  gStyle->SetOptStat(1111);// untk print no of events

 // root file
  TFile *f_sig = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_h125_4mu_zd50.root");
  TFile *f_bkg_zz4lep = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_Background_zz4lep.root");
  TFile *f_bkg_pph4l = new TFile("/Users/applestudio/Desktop/rootFILES/4LAnalyzer_Background_pph4l.root");
    
 // get branch in tree
  TH1F *h_sig = (TH1F*)f_sig->Get("mZb_4mu")->Clone("h_sig");
  TH1F *h_bkg_zz4lep = (TH1F*)f_bkg_zz4lep->Get("mZb_4mu")->Clone("h_bkg_zz4lep");
  TH1F *h_bkg_pph4l = (TH1F*)f_bkg_pph4l->Get("mZb_4mu")->Clone("h_bkg_pph4l");

 // book histogram for significance
 TH1D *S1= new TH1D("s1", "significane", 75, 0., 150.);
 S1->GetXaxis()->SetTitle("significane");
 S1->GetYaxis()->SetTitle("Number of Events");

 TH1D *S2= new TH1D("s1", "significane2", 75, 0., 150.);
 S2->GetXaxis()->SetTitle("significane2");
 S2->GetYaxis()->SetTitle("Number of Events");
    
    double nsigB = h_sig->Integral(h_sig->FindFixBin(15.), h_sig->FindFixBin(120.));
    cout << "signal nevent Before N is = " << nsigB << ".\n";
    
    double nbkg_zz4lepB = h_bkg_zz4lep->Integral(h_bkg_zz4lep->FindFixBin(15.), h_bkg_zz4lep->FindFixBin(120.));
    cout << "zz4lep nevent Before N is = " << nbkg_zz4lepB << ".\n";
    
    double nbkg_pph4lB = h_bkg_pph4l->Integral(h_bkg_pph4l->FindFixBin(15.), h_bkg_pph4l->FindFixBin(120.));
    cout << "pph4lep nevent Before N is = " << nbkg_pph4lB << ".\n";
    
 
 // normalize using scale to SM Higgs crossection & BR (total events 25k)
    h_sig->Scale(48.58*1.25e-4*250/25000);
    h_bkg_zz4lep->Scale(48.58*1.25e-4*250/25000);
    h_bkg_pph4l->Scale(48.58*1.25e-4*250/25000);
  
  //number of signal and background events determined by integrating over region of interest
  double nsig = h_sig->Integral(h_sig->FindFixBin(15.), h_sig->FindFixBin(120.));
  cout << "signal nevent is = " << nsig << ".\n";
  
  double nbkg_zz4lep = h_bkg_zz4lep->Integral(h_bkg_zz4lep->FindFixBin(15.), h_bkg_zz4lep->FindFixBin(120.));
  cout << "zz4lep nevent is = " << nbkg_zz4lep << ".\n";
  
  double nbkg_pph4l = h_bkg_pph4l->Integral(h_bkg_pph4l->FindFixBin(15.), h_bkg_pph4l->FindFixBin(120.));
  cout << "pph4lep nevent is = " << nbkg_pph4l << ".\n";
                     
  double S = nsig;
  double B = nbkg_zz4lep + nbkg_pph4l;

 // calculate significance
  double significance1 = 2*((sqrt(S+B))-(sqrt(B)));
    cout << "significance1 = " << significance1 << ".\n";
    S1 -> Fill(significance1);
  double significance2 = sqrt((2*(S+B)*log(1.0+S/B)-2*S));
    cout << "significance2 = " << significance2 << ".\n";
    S2 -> Fill(significance2);

    
                 
  //CREATE CANVAS
   TCanvas *c = new TCanvas("c","c");
    
  // CREATE AND OPEN NEW ROOT file to store info
    TFile fout("significance_4Mu_zd50.root", "RECREATE");
    S1->Write();
    S2->Write();
    
 // Clost the ROOT File
    fout.Close();
 
 // DRAW THE HISTOGRAM
    S1->Draw();
    S2->Draw();
     
}
