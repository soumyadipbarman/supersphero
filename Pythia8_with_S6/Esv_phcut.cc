//Esv_phcut.cc
//here vector array is used to save events

//Read the Pyhtia sample data and calculate the Eventshape variable with the phase space cut 
//Author :S.k.kundu   cerated on Fed 2018      update on 24 sep 2018
//Modified: S. Barman  dated: 15 Nov 2022

//Run this code with CLHEP along with root
//g++ Esv_phcut.cc EventShape_vector.cc -o Esv_phcut -I/home/soumyadip/Package/clhep/clhep2452/CLHEP_build -L/home/soumyadip/Package/clhep/clhep2452/CLHEP_install/lib -O -ansi -pedantic -Wall -D_GNU_SOURCE -std=c++11 -D_GNU_SOURCE -pthread -O2 -g -DNDEBUG `root-config --cflags`-Wl,-rpath,./ `root-config --glibs` -lCLHEP-Vector-2.4.5.2 -L/home/soumyadip/Package/Root/Root62606build/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -I/home/soumyadip/Package/Root/Root62606build/include/

//Run: ./Esv_phcut input.root output.root
//---------------------------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <TArrayC.h>
#include <string>
#include "TH1.h"                       // for histrograming
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TFile.h"                     // ROOT, for saving file.
#include  <string>
#include  <map>
#include "TTree.h"                     //for Tree file 
#include "TROOT.h"
#include "EventShape_vector.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
using namespace CLHEP;


int main(int argc, char* argv[] ){
  
  TFile *f = new TFile(argv[1]);  // load the root file 
  //TFile *f = new TFile("/home/suman/Fpythia/Final_Pythia_MCdata/PP_4vec_PThat80-7000_array.root");  // load the root file 
 
  TFile* outfile = new TFile(argv[2],"RECREATE"); // output file


//********************************************************************************************************

  TTree *Event_tree = (TTree*)f->Get("Event_tree");                       // define another new tree for Events
//  f.GetObject("Event_tree",Event_tree);


// define the variables
   Double_t weight;

   vector<Double_t> *Px = new vector<Double_t> ;
   vector<Double_t> *Py = new vector<Double_t>;
   vector<Double_t> *Pz = new vector<Double_t>;
   vector<Double_t> *Pt = new vector<Double_t>;
   vector<Double_t> *Ejet = new vector<Double_t>;

  Event_tree->SetBranchAddress("Px", &Px);
  Event_tree->SetBranchAddress("Py", &Py);
  Event_tree->SetBranchAddress("Pz", &Pz);
  Event_tree->SetBranchAddress("Pt", &Pt);
  Event_tree->SetBranchAddress("Ejet", &Ejet);


 Int_t Njets ;

  Event_tree->SetBranchAddress("Njets", &Njets);
  Event_tree->SetBranchAddress("weight", &weight);



//*****************************************************************************************
// Define the input jet vector
  vector<HepLorentzVector> Jet_object;

// Define EventShape vector
  vector<double> recovar;


//***********************************************************************************************

  
  static const int nvariable = 4;                   // number of variable studies
  int usedvars[nvariable] = {3, 9, 18, 33};     // four eventshape varaiable(Thurst,Transverse jet mass, jet boarddening , total jet mass)


  int nHTcut = 8;                              //Number of phase space cut used
  double phasecut[9] = {73, 93, 165, 225, 298, 365, 452, 557, 3000.0}; //The trigger value of the cut
  
  
  int ipt;
  //**************************************************************************************************
  //define the bins width 
  const int nmxbins = 38;
  
  int Nbins[nvariable]={33, 38, 33, 34};
  double binEdges[nvariable][nmxbins+1]={
    {-11.8, -10.8, -9.7, -8.9, -8.19, -7.56, -7.05, -6.62, -6.24, -5.9, -5.59, -5.31, -5.05, -4.8, -4.56, -4.33, -4.11, -3.9, -3.69, -3.49, -3.29, -3.09, -2.9, -2.71, -2.53, -2.35, -2.17, -2, -1.83, -1.67, -1.51, -1.36, -1.21, -1.07, 0, 0, 0, 0, 0},//thrust 3
    {-8.55, -8.1, -7.65, -7.3, -6.99, -6.71, -6.44, -6.18, -5.92, -5.67, -5.42, -5.18, -4.94, -4.7, -4.46, -4.23, -4, -3.77, -3.54, -3.31, -3.09, -2.87, -2.65, -2.44, -2.23, -2.03, -1.83, -1.64, -1.46, -1.28, -1.11, -0.95, -0.8, -0.66, -0.52, -0.39, -0.27, -0.16, -0.06 },//rho total 9
    /* {-12.2, -11.4, -10.7, -10.1, -9.49, -8.84, -8.24, -7.67, -7.13, -6.6, -6.09, -5.58, -5.08, -4.59, -4.12, -3.66, -3.23, -2.82, -2.44, -2.09, -1.78, -1.51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},    //y23 15*/
    {-6.30, -5.80, -5.40, -5.04, -4.74, -4.46, -4.19, -3.93, -3.68, -3.44, -3.21, -3, -2.8, -2.61, -2.43, -2.25, -2.08, -1.92, -1.77, -1.63, -1.49, -1.36, -1.23, -1.11, -0.99, -0.87, -0.76, -0.65, -0.54, -0.44, -0.34, -0.24, -0.14, -0.04, 0, 0, 0, 0, 0}, //jet broad 18
    {-9.8, -9.1, -8.5, -8.0, -7.59, -7.12, -6.71, -6.34, -6, -5.68, -5.38, -5.1, -4.83, -4.57, -4.31, -4.06, -3.82, -3.58, -3.35, -3.12, -2.9, -2.69, -2.48, -2.28, -2.09, -1.91, -1.73, -1.56, -1.4, -1.25, -1.11, -0.98, -0.86, -0.75, -0.65, 0, 1, 2, 3}  //Total trans jet 24
  };
  
  //****************************************************************************************************
  
  TH1F *Hist[10][10];      
  int Htcut[9]={73, 93, 165, 225, 298, 365, 452, 557, 3000};
  const char* typname[nvariable] = {"trans_thrust", "total_jet_mass", "jet_boardening", "trans_jet_mass"};
  char Hist_name[200];
  char Hist_title[200];
  for (int icut = 0; icut < nHTcut ; icut++){                  // books the histogram
    for ( int iivar = 0 ; iivar < nvariable ; iivar++){
      sprintf(Hist_name, "Esv_var%i_pt%i", usedvars[iivar], icut);
      sprintf(Hist_title, "%s_ht_%i-%i", typname[iivar], Htcut[icut], Htcut[icut+1]);
      Hist[iivar][icut] = new TH1F(Hist_name, Hist_title, Nbins[iivar], binEdges[iivar]);
      Hist[iivar][icut]->Sumw2();
    }
  }

  TH1F *Supersphero[10];
  for (int icut = 0; icut < nHTcut ; icut++){
  Supersphero[icut] = new TH1F(Form("SuperSphero_s6_pt%i", icut), Form("SuperSphero_ht_%i-%i", Htcut[icut], Htcut[icut+1]), 100,-8.0,0.5);
   
  }

  TH2F *S6_ncharge[10];
  for (int icut = 0; icut < nHTcut ; icut++){
  S6_ncharge[icut] = new TH2F(Form("s6_ncharge%i", icut), Form("S6_vs_Ncharge%i-%i", Htcut[icut], Htcut[icut+1]),10,.5,10.5, 100,0,1);

  }

  //************************************************************************************************************* 
//  Int_t jj=0; 
  Double_t  HT2;
  Int_t numjet;
  Double_t eventweight;
  Double_t  Pt1=0.0;
  Double_t Pt2=0.0;
  HepLorentzVector tmp3jet;
  Double_t px =0.0;
  Double_t py = 0.0;
  Double_t ptxy = 0.0;
  Double_t recterm = 0.0;
  Double_t PX=0.0, PY =0.0, PZ = 0.0 , E_Jet =0.0; 
  double tmppx ;
  double tmppy ;
  double tmppt ;
  double tmprec;
  
  
  //Start the Event loop to get the four vector of corresponding Jets

  Int_t nentries1 = Int_t(Event_tree->GetEntries());
  cout << "Total Event count ="<< nentries1 << endl;  



 for (Int_t i=0; i<nentries1 ; i++){                     //start event loop
    ipt =-1;
    numjet = 0;
    eventweight = 0.0;
    HT2=0.0;
    Pt1=0.0;
    Pt2=0.0;
    
    tmppx = 0.0;
    tmppy = 0.0;
    tmppt = 0.0;
    tmprec = 0.0;
 
    Event_tree->GetEntry(i);
    Jet_object.clear();
//    cout << "Even Loop = " << i << "     ";
//   cout << "jet size before entry = " << Jet_object.size() << endl;

    numjet = Njets;
    eventweight = weight;
//    cout << "number of jets = " << numjet<< endl;
//    cout << "weight of event= " << eventweight<< endl;
    
      //************************************************************************
//      cout << "HT2 before calcutation = " << HT2 << endl;
//      cout << "No of Jets after condition = " << numjet << endl;
      
      px =0.0;
      py = 0.0;
      ptxy = 0.0;
      recterm = 0.0;
      //tmp3jet.clear();

      for (Int_t ii =0 ; ii <numjet ; ii++){
	PX = (Px->at(ii));            //get the Px ,py ,pz form root file 
        PY = (Py->at(ii));
        PZ = (Pz->at(ii));
        E_Jet = (Ejet->at(ii));     


	HepLorentzVector tmp4v(PX,PY,PZ,E_Jet);                 // create a temporary fourvector
	
	  Jet_object.push_back(tmp4v);
     
	
	//check the phase space cut range for this event 
	if(ii==0){Pt1 = (Pt->at(ii));
  //         cout << "pt1 = " << Pt1  << endl;
          }
	if(ii==1){Pt2 = (Pt->at(ii));
    //             cout << "pt2 = " << Pt2  << endl;  
     }

//	cout << "HT2 = " << HT2 << endl;
	
      }    // end of jets loop
	HT2=((Pt1+Pt2)/2);



  //    cout << "after jet loop no of jet = " << Jet_object.size() <<endl; 
      
      if( HT2>=phasecut[0]){                                                      // Condition that events with HT2> 73 is only allowable          
	for(Int_t iipt =0 ; iipt <nHTcut+1 ; iipt++){                                                     // Loop to check the HT range
	  if(HT2 > phasecut[iipt] && HT2 < phasecut[iipt+1]){ipt = iipt;
//              cout << "save at   "<< ipt << "th phasesape"  <<endl;  
              };
	  
	}
	
	//////////////////////////////////////////////////////////////////////////////
	//Check the order of jet pt
	/*
	  vector<HepLorentzVector> Jet_object1;
	  Jet_object1.clear();
	  for(unsigned int  isize = 0; isize<Jet_object.size(); isize++){
	  for(unsigned  int  jsize =0 ; jsize < Jet_object.size() ; jsize++){
	  if (Jet_object[isize].perp() < Jet_object[jsize].perp()){Jet_object1.push_back(Jet_object[isize]);
	  Jet_object[isize] = Jet_object[jsize];
	  Jet_object[jsize] = Jet_object1[0];
	  Jet_object1.clear(); };
	  }
	  }                    //end of the loop for ordering of jet in order of pt 
	*/
	///////////////////////////////////////////////////////////////////////////////
	
       //Start EventShape variable calculation
	recovar.clear();
	EventShape_vector esv_pythia(Jet_object , 2.4, 1, 2, -1);
	
	recovar = esv_pythia.getEventShapes();
	//////////////////////////////////////////////////////////////////////
	if(ipt >=0){
	for( Int_t ifill = 0; ifill < nvariable ; ifill ++){
	  Hist[ifill][ipt]->Fill(recovar[usedvars[ifill]],eventweight);
	  cout << " | Var: " << usedvars[ifill] << " : " << recovar[usedvars[ifill]];
          	}
           Supersphero[ipt]->Fill(recovar[33], eventweight);
           S6_ncharge[ipt]->Fill(numjet,recovar[33], eventweight);
	 
	     
            } // condition for the phasespace check
          }// end of the condition of event with  HT2 >73  is allowable
cout << endl;
    
    }//end of event loop

////////////////////////////////////////////////////////////  
//calculation for rebin
int arrayind[4] ={3,9,18,33};

TH1F* rebin_hist(TH1F* thin, int ijetpt, int ivar);


TH1F *MC_hist_rebined[10][10];

for(Int_t iivar =0 ; iivar < 4 ; iivar++){
            for(Int_t icut =0; icut < 8 ; icut++){

//MC_hist_rebined[iivar][icut] = (TH1F*)(rebin_hist(Hist[iivar][icut] , icut , arrayind[iivar]));
   }
}

////////////////////////////////////////////////////////////////
  for (int icut = 0; icut < nHTcut ; icut++){                  // Write the histogram in output file
    for ( int iivar = 0 ; iivar < nvariable ; iivar++){
     // MC_hist_rebined[iivar][icut]->Write();
      Hist[iivar][icut]->Write();
    }
    Supersphero[icut]->Write();
    S6_ncharge[icut]->Write();
  }
                
  delete outfile;
  return 0;
  
} // end of main program

//////////////////////////////////////////////////////////////
//Function for rebin hist

TH1F* rebin_hist(TH1F* thin, int ijetpt, int ivar){
TH1F* thout;

cout << "Eventshape histogram rebined" <<endl;
int nxmod2 = 0;
double xmod2[200];

//int numHT =8;
//int numvar =4;

const char* namex;
char namey[100];
const char* titlex;
char titley[100];



// define the number of bin excluded form lower end

int arrayvar_first[4][8] = {{4,6,5,4,2,3,2,4},  // thrustc
                               {5,3,3,4,4,5,6,5}, //t3mass
                               {10,9,8,8,8,6,5,4}, //broadt
                               {8,7,7,7,7,4,5,5}}; //ttmass

// define the number of bin excluded from higher end

int arrayvar_last[4][8] =  {{4,3,3,3,3,4,5,4},  // thrustc
                               {4,4,5,5,8,7,7,7}, //t3mass
                               {7,6,7,9,5,12,11,10}, //broadt
                               {7,6,7,5,7,7,9,8}}; //ttmass

int arrayind[4] ={3,9,18,33};

  double yvl[100]={0.0};
  double erryvl[100]={0.0};
  int nbinx = thin->GetNbinsX();
  int indx=-1;
  for (int ij=0; ij<4; ij++) {
    if (ivar==arrayind[ij]) { indx=ij; break;}
  }
  if (indx<0) return thout;

  int ifirst = 0;
  nxmod2= -1;
  if (nxmod2<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_first[indx][ijetpt]) {
      } else if ( ij > (nbinx+2 - arrayvar_last[indx][ijetpt])) {
      } else {
        xmod2[ifirst] = thin->GetBinLowEdge(ij);
        xmod2[ifirst+1] = thin->GetBinLowEdge(ij+1);
        yvl[ifirst+1] =thin->GetBinContent(ij);
        erryvl[ifirst+1] +=thin->GetBinError(ij)*thin->GetBinError(ij);
        ifirst++;
      }
    }
    nxmod2 = ifirst;
  }

  namex = thin->GetName();
  sprintf(namey, "Roo_%s", namex);   // Roo means rebinded 
  titlex = thin->GetTitle();
  sprintf(titley, "Roo_%s", titlex);
  thout = new TH1F(namey, titley, nxmod2, xmod2);

  for (int ix=0; ix<nxmod2; ix++) {
    thout->SetBinContent(ix+1, yvl[ix+1]);
    thout->SetBinError(ix+1, sqrt(erryvl[ix+1]));
     }
/*  thout->GetXaxis()->SetTitleFont(42);
  thout->GetXaxis()->SetLabelFont(42);
  thout->GetXaxis()->SetLabelSize(0.07);
  thout->GetXaxis()->SetLabelOffset(.01);

  thout->GetYaxis()->SetTitleFont(42);
  thout->GetYaxis()->SetLabelFont(42);
  thout->GetYaxis()->SetLabelSize(0.07);
  thout->GetYaxis()->SetLabelOffset(.01);

  thout->GetXaxis()->SetTitleOffset(0.9);
  thout->GetXaxis()->SetTitleSize(0.07);

  thout->SetLineWidth(1);
  thout->GetXaxis()->SetTitle(thout->GetTitle());
*/
  return thout;

}   //end of rebin_hist functiom 
