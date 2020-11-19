/*========================
  Au Diya 27 October 2020
 *========================
 //   fighting corona   //
****************************
Analysis strategy is to get the hist for 4mu 4e 2mu2e
****************************/

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <utility>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "TSystem.h"
#include "math.h"

// for Root histogram
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"

// for Delphes
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

template<typename T>
void CollectionFilter(const TClonesArray& inColl ,vector<T*>& outColl, Double_t ptMin=30, Double_t etaMax=2.5)
{

  const TObject *object;

  for (Int_t i = 0; i < inColl.GetEntriesFast(); i++)
  {

   object = inColl.At(i);
   const T *t = static_cast<const T*>(object);

   if(t->P4().Pt() < ptMin) continue;
   if(TMath::Abs(t->P4().Eta()) > etaMax) continue;

   outColl.push_back(t);

  }
}
//------------------------------------------------------------------------

void Audiya4lAnalyzer(const char *fileName="/Users/applestudio/Desktop/sample28sept/ALLROOTFILES"){
  
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  //4e
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e10.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e15.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e20.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e25.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e30.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e35.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e40.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e45.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e50.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e55.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e60.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e65.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/e70.root");
  
  //4mu
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu10.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu15.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu20.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu25.root");
  chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu30.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu35.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu40.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu45.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu50.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu55.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu60.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu65.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/mu70.root");
 
  //2mu2e
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e10.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e15.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e20.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e25.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e30.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e35.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e40.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e45.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e50.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e55.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e60.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e65.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2mu2e70.root");
 
  //2e2mu
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu10.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu15.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu20.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu25.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu30.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu35.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu40.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu45.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu50.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu55.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu60.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu65.root");
  //chain.Add("/Users/applestudio/Desktop/sample28sept/ALLROOTFILES/2e2mu70.root");
  
    
    
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    
 // Declare variable
    
    const Electron *elec1, *elec2, *elec3, *elec4;
    const Muon     *muon1, *muon2, *muon3, *muon4;
 
    
    vector<const Electron*>   *electrons = new vector<const Electron*>();;
    vector<const Muon*>           *muons = new vector<const Muon*>();;
    
    Double_t mass, massEE, massMM;
    
    
    double mZ12, mZ34, mZ13, mZ24, mZ14, mZ23;
    double mass4mu, mass4e, mass2mu2e;
    double mZa, mZb;
    double ZaPT, ZaPhi, ZaEta;
    double ZbPT, ZbPhi, ZbEta;
    

    double mZ;

    double dZ12, dZ34, dZ13, dZ24, dZ14, dZ23;
    double dZc1, dZc2, dZc3;
    
    
    
    
    TLorentzVector p4Za, p4Zb, p4H;
    TLorentzVector Z12_V, Z34_V, Z23_V, Z14_V, Z24_V, Z13_V;
    
  // Book histograms & set axis lable
    
  //**Dimuon
  // ZTo2mu mass
   TH1D *h_mZ_2mu = new TH1D("massZto2muon", "mass of Z to 2 muon",120, 40., 120.);
   h_mZ_2mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
   h_mZ_2mu->GetYaxis()->SetTitle("Number of Events");

   // ZTo2e mass
   TH1D *h_mZ_2e = new TH1D("massZto2e", "mass of Z to 2e", 120, 40., 120.);
   h_mZ_2e->GetXaxis()->SetTitle("Invariant Mass for Nelectron=2 (in GeV/c^2)");
   h_mZ_2e->GetYaxis()->SetTitle("Number of Events");
    
   // These histograms are for 4 MUON reconstruction with different combinations
   // First combination: 1234
    // Mass Z12
    TH1D *h_mZ12_4mu = new TH1D("mZ12_4mu", "mass of Z12", 75, 0., 150.);
    h_mZ12_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
    h_mZ12_4mu->GetYaxis()->SetTitle("Number of Events");

    // Mass Z34
    TH1D *h_mZ34_4mu = new TH1D("mZ34_4mu", "mass of Z34", 75, 0., 150.);
    h_mZ34_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
    h_mZ34_4mu->GetYaxis()->SetTitle("Number of Events");

    // Second combination: 1324
    // Mass Z13
    TH1D *h_mZ13_4mu = new TH1D("mZ13_4mu", "mass of Z13", 75, 0., 150.);
    h_mZ13_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
    h_mZ13_4mu->GetYaxis()->SetTitle("Number of Events");

    // Mass Z24
    TH1D *h_mZ24_4mu = new TH1D("mZ24_4mu", "mass of Z24", 75, 0., 150.);
    h_mZ24_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon");
    h_mZ24_4mu->GetYaxis()->SetTitle("Number of Events");

    // Third combination: 1423
    // Mass Z14
    TH1D *h_mZ14_4mu = new TH1D("mZ14_4mu", "mass of Z14", 75, 0., 150.);
    h_mZ14_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
    h_mZ14_4mu->GetYaxis()->SetTitle("Number of Events");

    // Mass Z23
    TH1D *h_mZ23_4mu = new TH1D("mZ23_4mu", "mass of Z23", 75, 0., 150.);
    h_mZ23_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
    h_mZ23_4mu->GetYaxis()->SetTitle("Number of Events");

    // Mass Za: mass of ZTo2mu closest to Z mass
    TH1D *h_mZa_4mu = new TH1D("mZa_4mu", "mass Za", 120, 0., 120.);
    h_mZa_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
    h_mZa_4mu->GetYaxis()->SetTitle("Number of Events");

    // Mass Zb: mass of ZTo2mu not closest to Z mass
    TH1D *h_mZb_4mu = new TH1D("mZb_4mu", "mass Zb", 120, 0., 120.);
    h_mZb_4mu->GetXaxis()->SetTitle("Invariant mass of dimuonZto2Z (in GeV/c^2)");
    h_mZb_4mu->GetYaxis()->SetTitle("Number of Events");
    
    // 4muon mass o muon (full mass range)
    TH1D *h_m4_m4mu = new TH1D("mass4mu_full", "mass of 4 muon", 300, 0., 900.);
    h_m4_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
    h_m4_m4mu->GetYaxis()->SetTitle("Number of Events");
    
    // These histograms are for 4 ELECTORON reconstruction with different combinations
    
    // First combination: 1234
    // Mass Z12
    TH1D *h_mZ12_4e = new TH1D("mZ12_4e", "mass of Z12", 75, 0., 150.);
    h_mZ12_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
    h_mZ12_4e->GetYaxis()->SetTitle("Number of Events");

    // Mass Z34
    TH1D *h_mZ34_4e = new TH1D("mZ34_4e", "mass of Z34", 75, 0., 150.);
    h_mZ34_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
    h_mZ34_4e->GetYaxis()->SetTitle("Number of Events");

    // Second combination: 1324
    // Mass Z13
    TH1D *h_mZ13_4e = new TH1D("mZ13_4e", "mass of Z13", 75, 0., 150.);
    h_mZ13_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
    h_mZ13_4e->GetYaxis()->SetTitle("Number of Events");

    // Mass Z24
    TH1D *h_mZ24_4e = new TH1D("mZ24_4e", "mass of Z24", 75, 0., 150.);
    h_mZ24_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
    h_mZ24_4e->GetYaxis()->SetTitle("Number of Events");

    // Third combination: 1423
    // Mass Z14
    TH1D *h_mZ14_4e = new TH1D("mZ14_4e", "mass of Z14", 75, 0., 150.);
    h_mZ14_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
    h_mZ14_4e->GetYaxis()->SetTitle("Number of Events");

    // Mass Z23
    TH1D *h_mZ23_4e = new TH1D("mZ23_4e", "mass of Z23", 75, 0., 150.);
    h_mZ23_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
    h_mZ23_4e->GetYaxis()->SetTitle("Number of Events");

    // Mass Za: mass of Z closest to Z mass
    TH1D *h_mZa_4e = new TH1D("mZa_4e", "mass Za", 120, 0., 120.);
    h_mZa_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
    h_mZa_4e->GetYaxis()->SetTitle("Number of Events");

    // Mass Zb: mass of Z not closest to Z mass
    TH1D *h_mZb_4e = new TH1D("mZb_4e", "mass Zb", 120, 0., 120.);
    h_mZb_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
    h_mZb_4e->GetYaxis()->SetTitle("Number of Events");

    // 4electron mass (full mass range)
    TH1D *h_m4_m4e = new TH1D("mass4e_full", "mass of 4 electron", 300, 0., 900.);
    h_m4_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
    h_m4_m4e->GetYaxis()->SetTitle("Number of Events");
    
    // These histograms are for 2mu2e reconstruction with different combinations
    
    // Mass of Z to 2mu from 2mu2e
    TH1D *h_mZmu_2mu2e = new TH1D("massZmu_2mu2e", "mass Z2mu:2mu2e", 75, 0., 150.);
    h_mZmu_2mu2e->GetXaxis()->SetTitle("Invariant mass of Z1 (in GeV/c^2)");
    h_mZmu_2mu2e->GetYaxis()->SetTitle("Number of Events");

    // Mass of Z to 2e from 2mu2e
    TH1D *h_mZe_2mu2e = new TH1D("massZe_2mu2e", "mass Z2e:2mu2e", 75, 0., 150.);
    h_mZe_2mu2e->GetXaxis()->SetTitle("Invariant Mass of Z2 (in GeV/c^2)");
    h_mZe_2mu2e->GetYaxis()->SetTitle("Number of Events");

    // Mass Za: mass of Z1 closest to Z mass
    TH1D *h_mZa_2mu2e = new TH1D("mZa_2mu2e", "mass Z higher", 120, 0., 120.);
    h_mZa_2mu2e->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
    h_mZa_2mu2e->GetYaxis()->SetTitle("Number of Events");

    // Mass Zb: mass of Z2 not closest to Z mass
    TH1D *h_mZb_2mu2e = new TH1D("mZb_2mu2e", "mass Z lower", 120, 0., 120.);
    h_mZb_2mu2e->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
    h_mZb_2mu2e->GetYaxis()->SetTitle("Number of Events");

    // 2muons 2electrons mass spectrum full mass range
    TH1D *h_m4_m2mu2e = new TH1D("mass2mu2e_full", "mass of 2mu2e", 300, 0., 900.);
    h_m4_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
    h_m4_m2mu2e->GetYaxis()->SetTitle("Number of Events");
    
    //PT for Za & Zb
    /*TH1D *h_4mu_ZaPT = new TH1D("Za_4mu channel", "Pt of Za", 200, 0., 200.);
    h_4mu_ZaPT->GetXaxis()->SetTitle("Za_Pt");
    h_4mu_ZaPT->GetYaxis()->SetTitle("Number of Events");
    
    TH1D *h_4mu_ZbPT = new TH1D("Zb_4mu channel", "Pt of Zb", 200, 0., 200.);
    h_4mu_ZbPT->GetXaxis()->SetTitle("Zb_Pt");
    h_4mu_ZbPT->GetYaxis()->SetTitle("Number of Events");
    
    TH1D *h_4e_ZaPT = new TH1D("Za_4e channel", "Pt of Za", 200, 0., 200.);
    h_4e_ZaPT->GetXaxis()->SetTitle("Za_Pt");
    h_4e_ZaPT->GetYaxis()->SetTitle("Number of Events");
    
    TH1D *h_4e_ZbPT = new TH1D("Zb_4e channel", "Pt of Zb", 200, 0., 200.);
    h_4e_ZbPT->GetXaxis()->SetTitle("Zb_Pt");
    h_4e_ZbPT->GetYaxis()->SetTitle("Number of Events");
    
    TH1D *h_2mu2e_ZaPT = new TH1D("Za_2mu2e channel", "Pt of Za", 200, 0., 200.);
    h_2mu2e_ZaPT->GetXaxis()->SetTitle("Za_Pt");
    h_2mu2e_ZaPT->GetYaxis()->SetTitle("Number of Events");
    
    TH1D *h_2mu2e_ZbPT = new TH1D("Zb_2mu2e channel", "Pt of Zb", 200, 0., 200.);
    h_2mu2e_ZbPT->GetXaxis()->SetTitle("Zb_Pt");
    h_2mu2e_ZbPT->GetYaxis()->SetTitle("Number of Events");*/
    
    // Initialize variables
    // select largest, init -

    mZ12 = -9999.; mZ34 = -9999.; mZ13 = -9999.; mZ24 = -9999.; mZ14 = -9999.; mZ23 = -9999.;
  
   
    dZ12 = 9999.; dZ34 = 9999.; dZ13 = 9999.; dZ24 = 9999.; dZ14 = 9999.; dZ23 = 9999.; // select smallest, init +
    dZc1 = 9999.; dZc2 = 9999.; dZc3 = 9999.;

    mZa = -9999.; ZaPT = -9999.; ZaEta = -9999.; ZaPhi = -9999.;
    mZb = -9999.; ZbPT = -9999.; ZbEta = -9999.; ZbPhi = -9999.;
    mass4mu = -9999.;
    mass4e = -9999.;
    mass2mu2e = -9999.;
    
    p4Zb.SetPtEtaPhiM(0.,0.,0.,0.);
    p4Za.SetPtEtaPhiM(0.,0.,0.,0.);
    p4H.SetPtEtaPhiM(0., 0., 0., 0.);
    Z12_V.SetPtEtaPhiM(0.,0.,0.,0.);
    Z34_V.SetPtEtaPhiM(0.,0.,0.,0.);
    Z23_V.SetPtEtaPhiM(0.,0.,0.,0.);
    Z14_V.SetPtEtaPhiM(0.,0.,0.,0.);
    Z13_V.SetPtEtaPhiM(0.,0.,0.,0.);
    Z24_V.SetPtEtaPhiM(0.,0.,0.,0.);

    mZ = 91.1876;
    
    // Loop over all events begins
    for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      // print event number in the file
      if(entry%1000==0) cout << "Reading Event " << entry << endl;
      
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
    
      electrons -> clear();
      muons     -> clear();

       // Select leptons e||mu with pT > 5 || 7  GeV and |eta| < 2.5 || 2.4
       CollectionFilter(*branchElectron, *electrons , 5.0 , 2.5);
       CollectionFilter(*branchMuon    , *muons     , 7.0 , 2.4);
   
        //std::cout << "debug" << muons->size() << " " << electrons->size() << std::endl;
        
        //if ( muons->size)
    //sort the PT from highest to lowest
       // std::sort(elec.begin(), elec.end(), [] (const  double &electron1, const double &electron2) {return electron1.Pt() < electron2.Pt()});

   // sort Muon PT, High to low
      //  std::sort(muons.begin(), muons.end(), [] (const double  &mu1 , const  double &mu2) { return mu1.Pt() < mu2.Pt()});
    

   
    //=======================Zto2Muon Start=========================//
    
    if (muons->size() >= 2)
    {
        //define leading and subleading leptons
        muon1 = muons->at(0);
        muon2 = muons->at(1);
        //std::cout << "muon1 = " << muon1 << "muon2 = " << muon2 << std::endl;
        //select event with opposite charge leptons
        if(muon1->Charge + muon2->Charge == 0)
        {
        //define muon pair invariant mass
        massMM = ((muon1->P4()) + (muon2->P4())).M();
        h_mZ_2mu -> Fill(massMM);
        }
    }
    
    //=========================Zto2e Start===========================//
    
    
    if(electrons->size() >= 2)
    {
        //define leading and subleading leptons
        elec1 = electrons->at(0);
        elec2 = electrons->at(1);
        
        //select event with opposite charge leptons
        if(elec1->Charge + elec2->Charge == 0)
        {
        // Define electron pair invariant mass
        massEE = ((elec1->P4()) + (elec2->P4())).M();
        h_mZ_2e -> Fill(massEE);
        }
    }
    
    //==========================START ZZDto4Mu===============================//
    
    if (muons->size() >= 4)
    {
        muon1 = muons->at(0);
        muon2 = muons->at(1);
        muon3 = muons->at(2);
        muon4 = muons->at(3);
        
     // if (muon1->Charge + muon2->Charge + muon3->Charge + muon4->Charge == 0)
        if (muon1->Charge + muon2->Charge + muon3->Charge + muon4->Charge == 0)
        {
            //first combination: Combine muon 1234
            if (muon1->Charge + muon2->Charge == 0)
            {
                
                mZ12 = ((muon1->P4()) + (muon2->P4())).M();
                Z12_V = muon1->P4() + muon2->P4();
        

                  if (muon3->Charge + muon4->Charge == 0)
                {
                  
                mZ34 = ((muon3->P4()) + (muon4->P4())).M();
                Z34_V = muon3->P4() + muon4->P4();
                                             
                  if (mZ12 > 0.) h_mZ12_4mu->Fill(mZ12);
                  if (mZ34 > 0.) h_mZ34_4mu->Fill(mZ34);
                }
                } //close first combination mu1234

              dZ12 = std::abs( mZ12 - mZ );
              dZ34 = std::abs( mZ34 - mZ );

              // take the smallest difference between mass
              // to use for 4muon combination
              dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34;
            
            // Second combination: Combine muon 1324
                if (muon1->Charge + muon3->Charge == 0)
                  {
                      mZ13 = ((muon1->P4()) + (muon3->P4())).M();
                      Z13_V = muon1->P4() + muon3->P4();
                 
                    if (muon2->Charge + muon4->Charge == 0)
                  {
                  
                    mZ24 = ((muon2->P4()) + (muon4->P4())).M();
                    Z24_V = muon2->P4() + muon4->P4();
                    

                    if (mZ13 > 0.) h_mZ13_4mu->Fill(mZ13);
                    if (mZ24 > 0.) h_mZ24_4mu->Fill(mZ24);
                  }
                  } //close second combination mu1324

                dZ13 = std::abs( mZ13 - mZ );
                dZ24 = std::abs( mZ24 - mZ );

                dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24;
            
            // Third combination: Combine muon 1423
            if (muon1->Charge + muon4->Charge == 0)
              {
                 mZ14 = ((muon1->P4()) + (muon4->P4())).M();
                 Z14_V = muon1->P4() + muon4->P4();

                if (muon2->Charge + muon3->Charge == 0)
              {
                 
                mZ23 = ((muon2->P4()) + (muon3->P4())).M();
                Z23_V = muon2->P4() + muon3->P4();
                  
                if (mZ14 > 0.) h_mZ14_4mu->Fill(mZ14);
                if (mZ23 > 0.) h_mZ23_4mu->Fill(mZ23);
              }
              } // close thrid combination mu1423

            dZ14 = std::abs( mZ14 - mZ );
            dZ23 = std::abs( mZ23 - mZ );

            dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23;
            
            bool ptZadaug = false;

            if (dZc1 < dZc2 && dZc1 < dZc3)
              {
                if (dZ12 < dZ34)
              {
                
                mZa  = mZ12;
                ZaPT = Z12_V.Pt();
                ZaEta = Z12_V.Eta();
                ZaPhi = Z12_V.Phi();
                
                if (muon1->PT > 20. and muon2->PT > 10.)
                ptZadaug = true;
                  
                mZb  = mZ34;
                  ZbPT = Z34_V.Pt();
                  ZbEta = Z34_V.Eta();
                  ZbPhi = Z34_V.Phi();
                  
              }
                else
              {
              
                mZa  = mZ34;
                  ZaPT = Z34_V.Pt();
                  ZaEta = Z34_V.Eta();
                  ZaPhi = Z34_V.Phi();

                if (muon3->PT > 20. and muon4->PT > 10.)
                ptZadaug = true;

                mZb  = mZ12;
                  ZbPT = Z12_V.Pt();
                  ZbEta = Z12_V.Eta();
                  ZbPhi = Z12_V.Phi();
                  
              }
              } // close (dZc1 < dZc2 && dZc1 < dZc3)

            else if (dZc2 < dZc1 && dZc2 < dZc3)
              {
                if (dZ13 < dZ24)
              {
                mZa  = mZ13;
                  ZaPT = Z13_V.Pt();
                  ZaEta = Z13_V.Eta();
                  ZaPhi = Z13_V.Phi();

                if (muon1->PT > 20. and muon3->PT > 10.)
                ptZadaug = true;

                mZb  = mZ24;
                  ZbPT = Z24_V.Pt();
                  ZbEta = Z24_V.Eta();
                  ZbPhi = Z24_V.Phi();
                  
              }
                else
              {
                
                mZa  = mZ24;
                  ZaPT = Z24_V.Pt();
                  ZaEta = Z24_V.Eta();
                  ZaPhi = Z24_V.Phi();

                if (muon2->PT > 20. and muon4->PT > 10.)
                ptZadaug = true;

                mZb  = mZ13;
                  ZbPT = Z13_V.Pt();
                  ZbEta = Z13_V.Eta();
                  ZbPhi = Z13_V.Phi();
              }
              } // close (dZc2 < dZc1 && dZc2 < dZc3)

            else if (dZc3 < dZc1 && dZc3 < dZc2)
              {
                if (dZ14 < dZ23)
              {
          
                mZa  = mZ14;
                ZaPT = Z14_V.Pt();
                ZaEta = Z14_V.Eta();
                ZaPhi = Z14_V.Phi();

                if (muon1->PT > 20. and muon4->PT > 10.)
                ptZadaug = true;

                mZb  = mZ23;
                  ZbPT = Z23_V.Pt();
                  ZbEta = Z23_V.Eta();
                  ZbPhi = Z23_V.Phi();
                  
              }
                else
              {
                
                mZa  = mZ23;
                  ZaPT = Z23_V.Pt();
                  ZaEta = Z23_V.Eta();
                  ZaPhi = Z23_V.Phi();

                if (muon2->PT > 20. and muon3->PT > 10.)
                ptZadaug = true;

                mZb  = mZ14;
                  ZbPT = Z14_V.Pt();
                  ZbEta = Z14_V.Eta();
                  ZbPhi = Z14_V.Phi();
                  
              }
              } // close (dZc3 < dZc1 && dZc3 < dZc2)
            
            if (ptZadaug) {
                   if (mZa > 40. && mZa < 120.)
                   {
                       if (mZb > 12. && mZb < 120.){

                   h_mZa_4mu->Fill(mZa);
                   h_mZb_4mu->Fill(mZb);

                   // 4 vector
                   p4Za.SetPtEtaPhiM(ZaPT, ZaEta, ZaPhi, mZa);
                   p4Zb.SetPtEtaPhiM(ZbPT, ZbEta, ZbPhi, mZb);
                   p4H = p4Za + p4Zb;

                   mass4mu = p4H.M();
                         
                if (mass4mu > 70.)
                   {
                     
                     h_m4_m4mu->Fill(mass4mu);
                    
                   }
         } // close zb
        } //close za
               // h_4mu_ZaPT->Fill(ZaPT);
              //  h_4mu_ZbPT->Fill(ZbPT);
        } //end ptZadaug
        }//end 4mu charge
    } // close All 4mu
    //****************************************END 4MU************************************//
    
    //===================================START ZZDTo4e==================================//
    
    if (electrons->size() >= 4 )
    {
        elec1 = electrons->at(0);
        elec2 = electrons->at(1);
        elec3 = electrons->at(2);
        elec4 = electrons->at(3);
        
        if (elec1->Charge + elec2->Charge + elec3->Charge + elec4->Charge == 0)
        {
        //FIRST combination : combine elec 1234
        if (elec1->Charge + elec2->Charge == 0)
        {
           mZ12 = ((elec1->P4()) + (elec2->P4())).M();
           Z12_V = elec1->P4() + elec2->P4();
 
         if (elec3->Charge + elec4->Charge == 0)
        {
           mZ34 = ((elec3->P4()) + (elec4->P4())).M();
           Z34_V = elec3->P4() + elec4->P4();
            
         if (mZ12 > 0.) h_mZ12_4e->Fill(mZ12);
         if (mZ34 > 0.) h_mZ34_4e->Fill(mZ34);
            
        }
        } // close combination 1, elec 1234
            
        dZ12 = std::abs( mZ12 - mZ );
        dZ34 = std::abs( mZ34 - mZ );

        // take the smallest diff between mass to use for 4electron combination
        dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34;
            
    // Second combination: Combine elec 1324
        if (elec1->Charge + elec3->Charge == 0)
          {
            mZ13 = ((elec1->P4()) + (elec3->P4())).M();
            Z13_V = elec1->P4() + elec3->P4();
         
            if (elec2->Charge + elec4->Charge == 0)
          {
            mZ24 = ((elec2->P4()) + (elec4->P4())).M();
            Z24_V = elec2->P4() + elec4->P4();
            
            if (mZ13 > 0.) h_mZ13_4e->Fill(mZ13);
            if (mZ24 > 0.) h_mZ24_4e->Fill(mZ24);
          }
          } //close second combination 1324

        dZ13 = std::abs( mZ13 - mZ );
        dZ24 = std::abs( mZ24 - mZ );

        dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24;
        
        // Third combination: Combine elec 1423
        if (elec1->Charge + elec4->Charge == 0)
            
          {
             mZ14 = ((elec1->P4()) + (elec4->P4())).M();
             Z14_V = elec1->P4() + elec4->P4();
              
            if (elec2->Charge + elec3->Charge == 0)
          {
             mZ23 = ((elec2->P4()) + (elec3->P4())).M();
             Z23_V = elec2->P4() + elec3->P4();
            
            if (mZ14 > 0.) h_mZ14_4e->Fill(mZ14);
            if (mZ23 > 0.) h_mZ23_4e->Fill(mZ23);
          }
          } // close third combination 1423

        dZ14 = std::abs( mZ14 - mZ );
        dZ23 = std::abs( mZ23 - mZ );

        dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23;
        
            bool ptZadaug = false;

             // Now whichever have the smallest diff is considered the best comb.
             if (dZc1 < dZc2 && dZc1 < dZc3)
               {
                 if (dZ12 < dZ34)
               {
                 
                 mZa  = mZ12;
                   ZaPT = Z12_V.Pt();
                   ZaEta = Z12_V.Eta();
                   ZaPhi = Z12_V.Phi();
                 

                 if (elec1->PT > 20. and elec2->PT > 10.)
                 ptZadaug = true;

                 mZb  = mZ34;
                   ZbPT = Z34_V.Pt();
                   ZbEta = Z34_V.Eta();
                   ZbPhi = Z34_V.Phi();
                
               }
                 else
               {
        
                 mZa  = mZ34;
                   ZaPT = Z34_V.Pt();
                   ZaEta = Z34_V.Eta();
                   ZaPhi = Z34_V.Phi();

                 if (elec3->PT > 20. and elec4->PT > 10.)
                 ptZadaug = true;

                 mZb  = mZ12;
                   ZbPT = Z12_V.Pt();
                   ZbEta = Z12_V.Eta();
                   ZbPhi = Z12_V.Phi();
                   
               }
               } // close (dZc1 < dZc2 && dZc1 < dZc3)

             else if (dZc2 < dZc1 && dZc2 < dZc3)
               {
                 if (dZ13 < dZ24)
               {
                
                 mZa  = mZ13;
                   ZaPT = Z13_V.Pt();
                   ZaEta = Z13_V.Eta();
                   ZaPhi = Z13_V.Phi();
                   

                 if (elec1->PT > 20. and elec3->PT > 10.)
                 ptZadaug = true;

                 mZb  = mZ24;
                   ZbPT = Z24_V.Pt();
                   ZbEta = Z24_V.Eta();
                   ZbPhi = Z24_V.Phi();
               }
                 else
               {
                
                 mZa  = mZ24;
                   ZaPT = Z24_V.Pt();
                   ZaEta = Z24_V.Eta();
                   ZaPhi = Z24_V.Phi();

                 if (elec2->PT > 20. and elec4->PT > 10.)
                 ptZadaug = true;
                   
                 mZb  = mZ13;
                   ZbPT = Z13_V.Pt();
                   ZbEta = Z13_V.Eta();
                   ZbPhi = Z13_V.Phi();
               }
               } //close (dZc2 < dZc1 && dZc2 < dZc3)

             else if (dZc3 < dZc1 && dZc3 < dZc2)
               {
                 if (dZ14 < dZ23)
               {
            
                 mZa  = mZ14;
                   ZaPT = Z14_V.Pt();
                   ZaEta = Z14_V.Eta();
                   ZaPhi = Z14_V.Phi();

                 if (elec1->PT > 20. and elec4->PT > 10.)
                 ptZadaug = true;
        
                 mZb  = mZ23;
                 ZbPT = Z23_V.Pt();
                   ZbEta = Z23_V.Eta();
                   ZbPhi = Z23_V.Phi();
               }
                 else
               {
                 mZa  = mZ23;
                   ZaPT = Z23_V.Pt();
                   ZaEta = Z23_V.Eta();
                   ZaPhi = Z23_V.Phi();

                 if (elec2->PT > 20. and elec3->PT > 10.)
                 ptZadaug = true;

                 mZb  = mZ14;
                   ZbPT = Z14_V.Pt();
                   ZbEta = Z14_V.Eta();
                   ZbPhi = Z14_V.Phi();
               }
               } //close (dZc3 < dZc1 && dZc3 < dZc2)

             if (ptZadaug) {
               if (mZa > 40. && mZa < 120.) {
                 if (mZb > 12. && mZb < 120.) {
               h_mZa_4e->Fill(mZa);
               h_mZb_4e->Fill(mZb);
                     
               // Calculate 4 elec
               p4Za.SetPtEtaPhiM(ZaPT, ZaEta, ZaPhi, mZa);
               p4Zb.SetPtEtaPhiM(ZbPT, ZbEta, ZbPhi, mZb);
                     
               p4H = p4Za + p4Zb;

               mass4e = p4H.M();

               if (mass4e > 70.)
                 {
                   
                   h_m4_m4e->Fill(mass4e);
               
                 }
        } //end zb
        } //end za
               //  h_4e_ZaPT->Fill(ZaPT);
               //  h_4e_ZbPT->Fill(ZbPT);
        } // end ptZadaug
    
        }//close 4e charge
    } //close 4e
    
    //************************************** END ZZDTo4e ************************************//
    
    //=================================== START ZZD To2mu2e =================================//
    
    if (electrons->size() >= 2 && muons->size() >= 2 )
    {
        elec1 = electrons->at(0);
        elec2 = electrons->at(1);
        muon1 = muons->at(0);
        muon2 = muons->at(1);
        
        if (muon1->Charge + muon2->Charge + elec1->Charge + elec2->Charge == 0)
        {
        // For case 2mu2e, there is only 1 combination
          if (muon1->Charge + muon2->Charge == 0)
            {
              mZ12 = ((muon1->P4()) + (muon2->P4())).M();
              Z12_V = muon1->P4() + muon2->P4();
               
              if (elec1->Charge + elec2->Charge == 0)
            {
               mZ34 = ((elec1->P4()) + (elec2->P4())).M();
               Z34_V = elec1->P4() + elec2->P4();
              
              if (mZ12 > 0.) h_mZmu_2mu2e->Fill(mZ12);
              if (mZ34 > 0.) h_mZe_2mu2e->Fill(mZ34);

            }
            }
            dZ12 = std::abs(mZ12 - mZ); // mu
            dZ34 = std::abs(mZ34 - mZ); // e
            
            bool ptZadaug = false;

                 if (dZ12 < dZ34)
                   {
                     
                     mZa  = mZ12;
                       ZaPT = Z12_V.Pt();
                       ZaEta = Z12_V.Eta();
                       ZaPhi = Z12_V.Phi();

                     if (muon1->PT > 20. and muon2->PT > 10.)
                   ptZadaug = true;

                     mZb  = mZ34;
                       ZbPT = Z34_V.Pt();
                       ZbEta = Z34_V.Eta();
                       ZbPhi = Z34_V.Phi();
                   }
                 else
                   {
             
                        mZa = mZ34;
                       ZaPT = Z34_V.Pt();
                       ZaEta = Z34_V.Eta();
                       ZaPhi = Z34_V.Phi();

                     if (elec1->PT > 20. and elec2->PT > 10.)
                   ptZadaug = true;

                     mZb  = mZ12;
                     ZbPT = Z12_V.Pt();
                     ZbEta = Z12_V.Eta();
                     ZbPhi = Z12_V.Phi();
                 }

                 if (ptZadaug) {
                   if (mZa > 40. && mZa < 120.) {
                     if (mZb > 12. && mZb < 120.) {
                   h_mZa_2mu2e->Fill(mZa);
                   h_mZb_2mu2e->Fill(mZb);

                   // Now combine these 2 muons and 2 electrons
            
                   // Calculate 4 lepton: 2muons 2electrons

                   p4Za.SetPtEtaPhiM(ZaPT, ZaEta, ZaPhi, mZa);
                   p4Zb.SetPtEtaPhiM(ZbPT, ZbEta, ZbPhi, mZb);
                   p4H = p4Za + p4Zb;
                         
                   mass2mu2e = p4H.M();
                         
                   if (mass2mu2e > 70.)
                     {
                     
                       h_m4_m2mu2e->Fill(mass2mu2e);
                      
                     }
                     } //close zb
                   }//close za
                  //   h_2mu2e_ZaPT->Fill(ZaPT);
                  //   h_2mu2e_ZbPT->Fill(ZbPT);
                 }//close ptZadaug
            
        }//close charge
       
    } // close 2mu2e
    
    //*************************************ZZdTo2e2mu end********************************//
} // END EVENT LOOP

  //CREATE CANVAS
    TCanvas *c = new TCanvas("C1","4LANALYZER",200,200,600,600);
    
  // CREATE AND OPEN NEW ROOT file to store info
    TFile fout("4LAnalyzer_4mu_zd30.root", "RECREATE");
    
 //Write the histograms
 //dimuon pair
 h_mZ_2mu->Write();
 h_mZ_2e->Write();
 //4mu
 h_mZ12_4mu->Write();
 h_mZ34_4mu->Write();
 h_mZ13_4mu->Write();
 h_mZ24_4mu->Write();
 h_mZ14_4mu->Write();
 h_mZ23_4mu->Write();
 h_mZa_4mu->Write();
 h_mZb_4mu->Write();
 h_m4_m4mu->Write();
 //4e
 h_mZ12_4e->Write();
 h_mZ34_4e->Write();
 h_mZ13_4e->Write();
 h_mZ24_4e->Write();
 h_mZ14_4e->Write();
 h_mZ23_4e->Write();
 h_mZa_4e->Write();
 h_mZb_4e->Write();
 h_m4_m4e->Write();
 //2e2mu
 h_mZmu_2mu2e->Write();
 h_mZe_2mu2e->Write();
 h_mZa_2mu2e->Write();
 h_mZb_2mu2e->Write();
 h_m4_m2mu2e->Write();
 // Pt for Za & Zb
 /*h_4mu_ZaPT->Write();
 h_4mu_ZbPT->Write();
 h_4e_ZaPT->Write();
 h_4e_ZbPT->Write();
 h_2mu2e_ZaPT->Write();
 h_2mu2e_ZbPT->Write();*/

    
 // Clost the ROOT File
    fout.Close();
 
//Draw the histogram
    //dimuon pair 2
    h_mZ_2mu->Draw();
    h_mZ_2e->Draw();
    //4mu 9
    h_mZ12_4mu->Draw();
    h_mZ34_4mu->Draw();
    h_mZ13_4mu->Draw();
    h_mZ24_4mu->Draw();
    h_mZ14_4mu->Draw();
    h_mZ23_4mu->Draw();
    h_mZa_4mu->Draw();
    h_mZb_4mu->Draw();
    h_m4_m4mu->Draw();
    //4e 9
    h_mZ12_4e->Draw();
    h_mZ34_4e->Draw();
    h_mZ13_4e->Draw();
    h_mZ24_4e->Draw();
    h_mZ14_4e->Draw();
    h_mZ23_4e->Draw();
    h_mZa_4e->Draw();
    h_mZb_4e->Draw();
    h_m4_m4e->Draw();
    //2e2mu 5
    h_mZmu_2mu2e->Draw();
    h_mZe_2mu2e->Draw();
    h_mZa_2mu2e->Draw();
    h_mZb_2mu2e->Draw();
    h_m4_m2mu2e->Draw();
    //Pt for Za & Zb
   /* h_4mu_ZaPT->Draw();
    h_4mu_ZbPT->Draw();
    h_4e_ZaPT->Draw();
    h_4e_ZbPT->Draw();
    h_2mu2e_ZaPT->Draw();
    h_2mu2e_ZbPT->Draw();*/
    
} // END VOID
 
