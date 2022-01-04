#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"



void plot_various() {

  const char* inputfilename="test.root";
  const char* outputfilename="test_hist.root";

  TH1F *hTotalE = new TH1F("hTotalE","energy depisited /true",600,0.,1.1);
  TH1F *hWorldE = new TH1F("hWorldE","energy total world /true",600,0.,1.1);
  TH2F *h_nonel = new TH2F("h_nonel","number inelastic vs proton E deposit",
			   100,0.,5.,600,0.,600);


  TFile *f = new TFile(inputfilename);
  TTree *t1 = (TTree*)f->Get("tree");

  vector<float> *inputMomentum = new vector<float>;
  float Ninelastic;
  float depositedEnergyECAL_f;
  float kineticEnergyEscapeWorld;
  float depositedEnergyTotal,depositedEnergyWorld;
  float depositedEnergyECAL_absorb_f_particleID[8];

  t1->SetBranchAddress("inputMomentum",&inputMomentum);
  t1->SetBranchAddress("kineticEnergyEscapeWorld",&kineticEnergyEscapeWorld);
  t1->SetBranchAddress("depositedEnergyTotal",&depositedEnergyTotal);
  t1->SetBranchAddress("depositedEnergyECAL_f",&depositedEnergyECAL_f);
  t1->SetBranchAddress("depositedEnergyECAL_absorb_f_particleID",&depositedEnergyECAL_absorb_f_particleID);
  t1->SetBranchAddress("Ninelastic",&Ninelastic);


  Int_t nentries = (Int_t)t1->GetEntries();
  for(Int_t i=0;i<nentries; i++) {
    t1->GetEntry(i);
    float trueE=9999999.;
    if((*inputMomentum)[3]>0) trueE=(*inputMomentum)[3];
    float Eabs=depositedEnergyECAL_f;

    
    std::cout<<endl<<"event number "<<i<<std::endl;
    std::cout<<(*inputMomentum)[0]<<","<<(*inputMomentum)[1]<<","<<(*inputMomentum)[2]<<","<<(*inputMomentum)[3]<<std::endl;

    
    std::cout<<"total energy deposited is "<<depositedEnergyTotal<<std::endl;
    std::cout<<"kinetic energy escape deposited is "<<kineticEnergyEscapeWorld<<std::endl;



    hTotalE->Fill(depositedEnergyTotal/trueE);
    hWorldE->Fill((depositedEnergyTotal+kineticEnergyEscapeWorld)/trueE);
    h_nonel->Fill(depositedEnergyECAL_absorb_f_particleID[7],Ninelastic);
		    
  }

  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hTotalE->Write();
  hWorldE->Write();
  h_nonel->Write();
  out->Close();

}


