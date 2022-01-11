#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"



void plot_various() {

  bool idebug = 1;


  const char* inputfilename="test.root";
  const char* outputfilename="test_hist.root";

  TH1F *hTotalE = new TH1F("hTotalE","energy depisited /true",600,0.,1.1);
  TH1F *hWorldE = new TH1F("hWorldE","energy total world /true",600,0.,1.1);
  TH2F *h_nonel = new TH2F("h_nonel","number inelastic vs proton E deposit",
			   100,0.,5.,600,0.,600);
  TH2F *h_nonel2 = new TH2F("h_nonel2","number inelastic vs  deposit after 4ns",
			   100,0.,5.,600,0.,600);

  TFile *f = new TFile(inputfilename);
  TTree *t1 = (TTree*)f->Get("tree");

  vector<float> *inputMomentum = new vector<float>;
  float Ninelastic;
  float depositedEnergyECAL_f;
  float kineticEnergyEscapeWorld;
  float depositedEnergyTotal,depositedEnergyWorld;
  float depositedIonEnergyTotal;
  float depositedEnergyECAL_absorb_f_particleID[9];
  float depositedIonEnergyECAL_pidtime[9][80];

  t1->SetBranchAddress("inputMomentum",&inputMomentum);
  t1->SetBranchAddress("kineticEnergyEscapeWorld",&kineticEnergyEscapeWorld);
  t1->SetBranchAddress("depositedEnergyTotal",&depositedEnergyTotal);
  t1->SetBranchAddress("depositedIonEnergyTotal",&depositedIonEnergyTotal);
  t1->SetBranchAddress("depositedEnergyECAL_f",&depositedEnergyECAL_f);
  t1->SetBranchAddress("depositedEnergyECAL_absorb_f_particleID",&depositedEnergyECAL_absorb_f_particleID);
  t1->SetBranchAddress("depositedIonEnergyECAL_pidtime",&depositedIonEnergyECAL_pidtime);
  t1->SetBranchAddress("Ninelastic",&Ninelastic);


  float acheck[9];

  Int_t nentries = (Int_t)t1->GetEntries();
  for(Int_t i=0;i<nentries; i++) {
    t1->GetEntry(i);
    float trueE=9999999.;
    if((*inputMomentum)[3]>0) trueE=(*inputMomentum)[3];
    float Eabs=depositedEnergyECAL_f;

    
    if(idebug) std::cout<<endl<<"event number "<<i<<std::endl;
    if(idebug) std::cout<<(*inputMomentum)[0]<<","<<(*inputMomentum)[1]<<","<<(*inputMomentum)[2]<<","<<(*inputMomentum)[3]<<std::endl;

    
    if(idebug) std::cout<<"total energy deposited is "<<depositedEnergyTotal<<std::endl;    
    if(idebug) std::cout<<"total ionizing energy deposited is "<<depositedIonEnergyTotal<<std::endl;
    if(idebug) std::cout<<"kinetic energy escape deposited is "<<kineticEnergyEscapeWorld<<std::endl;


    for (int j=0;j<9;j++) {
      acheck[j]=0.;
    }
    float alate=0.;
    for (int j=0;j<9;j++) {
      for (int k=0;k<80;k++) {
	acheck[j]+=depositedIonEnergyECAL_pidtime[j][k];
      }
      for (int k=15;k<80;k++) {
	alate+=depositedIonEnergyECAL_pidtime[j][k];
      }
    }


    if(idebug) std::cout<<"energy pi- "<<depositedEnergyECAL_absorb_f_particleID[0]<<" "<<acheck[0]<<std::endl;
    if(idebug) std::cout<<"energy e+ "<<depositedEnergyECAL_absorb_f_particleID[1]<<" "<<acheck[1]<<std::endl;
    if(idebug) std::cout<<"energy e- "<<depositedEnergyECAL_absorb_f_particleID[2]<<" "<<acheck[2]<<std::endl;
    if(idebug) std::cout<<"energy gamma "<<depositedEnergyECAL_absorb_f_particleID[3]<<" "<<acheck[3]<<std::endl;
    if(idebug) std::cout<<"energy pi+ "<<depositedEnergyECAL_absorb_f_particleID[4]<<" "<<acheck[4]<<std::endl;
    if(idebug) std::cout<<"energy K "<<depositedEnergyECAL_absorb_f_particleID[5]<<" "<<acheck[5]<<std::endl;
    if(idebug) std::cout<<"energy n "<<depositedEnergyECAL_absorb_f_particleID[6]<<" "<<acheck[6]<<std::endl;
    if(idebug) std::cout<<"energy p "<<depositedEnergyECAL_absorb_f_particleID[7]<<" "<<acheck[7]<<std::endl;
    if(idebug) std::cout<<"energy other "<<depositedEnergyECAL_absorb_f_particleID[8]<<" "<<acheck[8]<<std::endl;
    float asum=0.;
    for (int j=0;j<9;j++ ) {
      asum+=depositedEnergyECAL_absorb_f_particleID[j];
    }
    if(idebug) std::cout<<" sum over particles "<<asum<<std::endl;
    if(fabs(asum-depositedEnergyTotal)>0.01) std::cout<<" energy conservation problem "<<depositedEnergyTotal<<" "<<asum<<std::endl;


    hTotalE->Fill(depositedEnergyTotal/trueE);
    hWorldE->Fill((depositedEnergyTotal+kineticEnergyEscapeWorld)/trueE);
    h_nonel->Fill(depositedEnergyECAL_absorb_f_particleID[7],Ninelastic);
    h_nonel2->Fill(alate,Ninelastic);
		    
  }

  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hTotalE->Write();
  hWorldE->Write();
  h_nonel->Write();
  h_nonel2->Write();
  out->Close();

}


