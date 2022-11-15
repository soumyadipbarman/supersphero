//HardQCD.cc
//histogram for basic distribution and tree for 4 vector
//CM energy 13Tev
//use the fastjet for jet finding

//Author : S.k.kundu
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/**

g++ HardQCD.cc -o HardQCD -I/home/soumyadip/Package/fastjet/fastjet-install/include -I/home/soumyadip/Package/Root/Root62606build/include/ -I../include -O2  -pedantic -W -Wall -Wshadow -fPIC -L../lib -Wl,-rpath,../lib -lpythia8 -ldl  `root-config --cflags` -Wl,-rpath, `/home/soumyadip/Package/Root/Root62606build/bin/root-config --glibs` -L/home/soumyadip/Package/fastjet/fastjet-install/lib -Wl,-rpath,/home/soumyadip/Package/fastjet/fastjet-install/lib -lfastjet

**/
/////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"                   // This is the minimal interface needed to access FastJet.
#include "fastjet/ClusterSequenceArea.hh"         //clusterSequence with the area support
#include "fastjet/ClusterSequence.hh"
#include "TH1.h"                                 // for histrograming
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"                     // ROOT, for saving file.
#include "TTree.h"                     //for Tree file 
#include "TROOT.h"
#include <vector>
using namespace std;
using namespace fastjet;
using namespace Pythia8;
///////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
  
  // Create Pythia instance and set it up to generate hard QCD processes
  Pythia pythia;                          // Generator
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");

  pythia.readString("HardQCD:all = on");  // 
  pythia.readString("PhaseSpace:pTHatMin = 15");
  
  pythia.readString("PhaseSpace:pTHatMax = 7000.0");
  pythia.readString("PhaseSpace:bias2Selection = on");
  pythia.readString("PhaseSpace:bias2SelectionPow = 4.0");

  pythia.readString("Tune:pp = 14");
  pythia.readString("Tune:ee = 7");
  pythia.readString("Beams:eCM = 13000.");
  
  pythia.readFile(argv[1]);
  pythia.init();
  
  TFile* outFile = new TFile(argv[2], "RECREATE");
  int nEvent = pythia.mode("Main:numberOfEvents");
//  int nEvent = 100;
//  int nEvent = argi[1];
////////////////////////////////////////////////////////////////////////
  
//Define and initialaize an Tree file 
  TTree *Event_tree = new TTree("Event_tree","Event_Information");
  
//Define the Variable for branch
  Double_t PX, PY, PZ, P_T, E_Jet, eta, phi, weight;
//   Double_t weight;  
 vector<Double_t> Px, Py, Pt, Pz, Eta, Phi, Ejet;
 
 Int_t  Njets , sjtevn=0 ,djtevn =0,mjtevn=0;
  
  /*
  Event_tree->Branch("Njets", &Njets, "Njets/I");
  Event_tree->Branch("Px", &Px, "Px[30]/D");
  Event_tree->Branch("Py", &Py, "Py[30]/D");
  Event_tree->Branch("Pz", &Pz, "Pz[30]/D");
  Event_tree->Branch("Pt", &Pt, "Pt[30]/D");
  Event_tree->Branch("Ejet", &Ejet, "Ejet[30]/D");
  Event_tree->Branch("Eta", &Eta, "Eta[30]/D");
  Event_tree->Branch("Phi", &Phi, "Phi[30]/D");
  Event_tree->Branch("weight", &weight, "weight/D");
  */
  Event_tree->Branch("Px", "vector<Double_t>", &Px);
  Event_tree->Branch("Py", "vector<Double_t>", &Py);
  Event_tree->Branch("Pz", "vector<Double_t>", &Pz);
  Event_tree->Branch("Pt", "vector<Double_t>", &Pt);
  Event_tree->Branch("Eta", "vector<Double_t>", &Eta);
  Event_tree->Branch("Phi", "vector<Double_t>", &Phi);
  Event_tree->Branch("Ejet", "vector<Double_t>", &Ejet);
  Event_tree->Branch("weight", &weight, "weight/D");
  Event_tree->Branch("Njets", &Njets, "Njets/I");



  
  TH1D *PT_inclusive = new TH1D("PT_inclusive","PT_of_all_jets",400,20,2020.0);
  PT_inclusive->Sumw2();
  TH1D *Eta_inclusive = new TH1D("Eta_inclusive","Eta_of all_Jets",100,-2.5, 2.5);
  Eta_inclusive->Sumw2();
  TH1D *Phi_inclusive = new TH1D("phi_inclusive","phi_of_all_Jets",100,-3.14,3.14);
  Phi_inclusive->Sumw2();
  TH1D *HT2=new TH1D("HT2","HT2", 400,20.0 ,2020.0);
  HT2->Sumw2();
  TH1D *PT_1st_jet = new TH1D("PT_1st_jet","PT_of_1st_Jet",400,20.0,2020);
  PT_1st_jet->Sumw2();
  TH1D *Eta_1st_jet = new TH1D("Eta_1st_jet","Eta of 1st Jet",100,-2.5, 2.5);
  Eta_1st_jet->Sumw2();
  TH1D *PT_2nd_jet = new TH1D("PT_2nd_jet","PT_of_second_leading_jet",400,20.0,2020.0);
  PT_2nd_jet->Sumw2();
  TH1D *Eta_2nd_jet =  new TH1D("second_jet_eta","second_jet_eta",100,-2.5, 2.5);
  Eta_2nd_jet->Sumw2();
  TH1D *Delta_PT12 = new TH1D("Delta_PT12","Delta_PT_of_ist_and Second_jets",100,20.0,500.0);
  Delta_PT12->Sumw2();
  TH1D *Delta_phi12 = new TH1D("Delta_phi12","Delta_phi_of_ist_and_2nd_jets",100,-3.14,3.14);
  Delta_phi12->Sumw2();
  TH1D *comppt1pt2= new TH1D("comppt1pt2","Component_of_2nd_jet_on_leading_jet",60,0.0,1.0);
  comppt1pt2->Sumw2();
  TH1D *Numjet= new TH1D("numjet","Number_of_jets",60,0.0,30);
  Numjet->Sumw2();
  TH1D *Numjet1= new TH1D("numjet1","Number_of_all_jet_without_cuts",60,0.0,30);
  Numjet1->Sumw2();
  
  // Select common parameters for FastJet analyses.
  double R = 0.4;                        
  
  //start the fastjet analysis 
  
  JetDefinition jet_def(antikt_algorithm, R);              //Choose the jet difinition in :jet_def
  
  vector<PseudoJet> fjInputs;   //Define Fastjet input :fjInput
  
  //*********************************************************************************
  // Begin event loop. Generate event; skip if generation aborted.
  
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {                          /*start event loop*/
  if (!pythia.next()) continue;
   
   fjInputs.resize(0);   //Reset Fastjet input 

   //cout <<pythia.event.size()<< endl;                     This give the total number of particle created in the events
   
   for (int i = 0 ; i < pythia.event.size() ; ++i){
  // cout << "eta = " << pythia.event[i].eta() <<endl;
     if (!pythia.event[i].isFinal())        continue;  //allow final state particle only
    fjInputs.push_back(PseudoJet(pythia.event[i].px(), pythia.event[i].py(),pythia.event[i].pz(), pythia.event[i].e()));            
      
   } 
   
   //check the event  
   if (fjInputs.size() == 0){
     cout << "Error:event with no final state particles"<<endl;
     continue; 
   }
    
    //Run Fastjet algorithm//**********************************************************************************
   ClusterSequence clust_seq(fjInputs, jet_def);  // run the jet clustering with the above jet 
   
   //sorted the jet with the minimum pT
   double ptmin =30.0;
   Selector jet_selector = SelectorAbsEtaMax(4.7);
   vector<PseudoJet> inclusive_jets = sorted_by_pt(jet_selector(clust_seq.inclusive_jets(ptmin)));
   
   //vector <PseudoJet> inclusive_jets = sorted_by_pt(clustSeq.inclusive_jets(ptmin));   // Define inclusive_jets as fastjet object
   //vector <PseudoJet> inclusive_jets = clustSeq.inclusive_jets();
   //**********************************************************************************************************
   
   //cout << "Clustering with " << jet_def.description() << endl;

   Px.clear(); 
    Py.clear();
     Pz.clear();
      Pt.clear(); 
       Ejet.clear(); 
        Eta.clear(); 
         Phi.clear();
            

 
  Double_t pt1= 0.0, pt2 = 0.0;
  Double_t phi1 = 0.0, phi2 = 0.0;
  Double_t Deltapt = 0.0, Deltaphi = 0.0, ht2 = 0.0;
  Double_t pt2_pt1 = 0.0;

   Njets =0;
   int numjet = inclusive_jets.size();                    //number of jets in any events 
   weight=0.0;
   P_T = 0.0; PX = 0.0; PY = 0.0; PZ = 0.0; eta =0.0; E_Jet =0.0; phi = 0.0;
   Numjet1->Fill(numjet,weight);

////////////////////////////////////////////
     if(inclusive_jets.size()<2) continue ;
    
     if((inclusive_jets[0].perp()) < 30 || fabs(inclusive_jets[0].eta()) >2.4) continue;
     if((inclusive_jets[1].perp()) < 30 || fabs(inclusive_jets[1].eta()) >2.4) continue;
     if(((inclusive_jets[0].perp() + inclusive_jets[1].perp())/2) < 73.0 ) continue;

///////////////////////////////////////////////
 if(numjet ==1){sjtevn++;}
 if(numjet ==2){djtevn++;}
 if(numjet >2){mjtevn++;}

///////////////////////////////////////////////
  if(pythia.info.weight() < 0) continue; 
   weight = pythia.info.weight();
  
//for ISR FSR issue check
//if (((inclusive_jets[0].perp() + inclusive_jets[1].perp())/2) > 165.0 && ((inclusive_jets[0].perp() + inclusive_jets[1].perp())/2) <225.0 ){
//	if(weight>1) cout << "Evt No. = " << iEvent << ": jets count = " << numjet << ": HT2 = "<<(inclusive_jets[0].perp() + inclusive_jets[1].perp())/2 <<setprecision(10) << "  :weight = " << weight << endl;
//}

  
//    cout << "Evt No. = " << iEvent << ": jets count = " << numjet << setprecision(10) << "  :weight = " << weight << endl; 
    int jjet = 0;
  
     for (int i = 0; i < numjet; ++i) {              //start jets loop
     if((inclusive_jets[i].perp()) < 30 || fabs(inclusive_jets[i].eta()) >2.4) continue;

   jjet++; 

     P_T = inclusive_jets[i].perp();              //cout << P_T << endl;
     eta = inclusive_jets[i].eta();
     phi = inclusive_jets[i].phi_std();
     PX = inclusive_jets[i].px();
     PY = inclusive_jets[i].py();
     PZ = inclusive_jets[i].pz();
     E_Jet = inclusive_jets[i].e();
  
   
  Pt.push_back(P_T);
  Px.push_back(PX);
  Py.push_back(PY);
  Pz.push_back(PZ);
  Eta.push_back(eta);
  Phi.push_back(phi);
  Ejet.push_back(E_Jet);


     PT_inclusive->Fill(P_T,weight);
     Eta_inclusive->Fill(eta,weight);
     Phi_inclusive->Fill(phi,weight);
   
     if(jjet==1){
	PT_1st_jet->Fill(P_T,weight);
        Eta_1st_jet->Fill(eta,weight);
        pt1 = P_T;
        phi1 = phi;	      
	}

    if(jjet==2){       
        PT_2nd_jet->Fill(P_T,weight);
        Eta_2nd_jet->Fill(eta,weight);
        pt2 = P_T;
        phi2 = phi;
        }
 

   } //end jets loop
   
  Njets = jjet;  
  


   Event_tree->Fill();

    ht2 = (pt1+pt2)/2;
    HT2->Fill(ht2,weight);
    Deltapt = (pt1 - pt2);
    Deltaphi = (phi1 - phi2);
    pt2_pt1 = ((pt2 * sin(Deltaphi))/pt1);
   
   Delta_PT12->Fill(Deltapt,weight);
   Delta_phi12->Fill(Deltaphi,weight);
   comppt1pt2->Fill(pt2_pt1,weight);   
   Numjet->Fill(Njets,weight);
  }                                                  /*end of the event loop*/
 /*****************************************************************************/
 pythia.stat();
 
 Event_tree->Write();

 PT_inclusive->Write();
 Eta_inclusive->Write();
 Phi_inclusive->Write();
 PT_1st_jet->Write();
 Eta_1st_jet->Write();
 Eta_2nd_jet->Write();
 PT_2nd_jet->Write();
 HT2->Write();
 Delta_PT12->Write();
 Delta_phi12->Write();
 comppt1pt2->Write();
 Numjet->Write();
 Numjet1->Write();

 cout << "sinjet = " << sjtevn << "  doublejt = " << djtevn << " mjtevn = " <<mjtevn  <<endl;
 delete outFile;
 //cout << "Clustering with " << jet_def.description() << endl;
 return 0;

}
