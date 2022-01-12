#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"



void plot_various(char* input) {

  bool idebug = 1;

  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  myStyle->SetOptStat(0);



  const char* inputfilename=input;
  const char* outputfilename="hist.root";

  TH1F *hTotalE = new TH1F("hTotalE","energy depisited /true",600,0.,1.1);
  TH1F *hWorldE = new TH1F("hWorldE","energy total world /true",600,0.,1.1);
  TH2F *h_nonel = new TH2F("h_nonel","number inelastic vs proton E deposit",
			   100,0.,5.,600,0.,600);
  TH2F *h_nonel2 = new TH2F("h_nonel2","number inelastic vs  deposit after 1.25ns",
			   100,0.,10.,600,0.,600);

  TH2F *h_nonel3 = new TH2F("h_nonel3","number inelastic vs  total E deposition",
			   100,0.,15.,600,0.,600);



  TH1F *h_timepim = new TH1F("h_timepim","",80,0.,20.);
  TH1F *h_timeep = new TH1F("h_timeep","",80,0.,20.);
  TH1F *h_timeem = new TH1F("h_timeem","",80,0.,20.);
  TH1F *h_timegam = new TH1F("h_timegam","",80,0.,20.);
  TH1F *h_timepip = new TH1F("h_timepip","",80,0.,20.);
  TH1F *h_timek = new TH1F("h_timek","",80,0.,20.);
  TH1F *h_timen = new TH1F("h_timen","",80,0.,20.);
  TH1F *h_timep = new TH1F("h_timep","",80,0.,20.);
  TH1F *h_timeo = new TH1F("h_timeo","",80,0.,20.);



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
  float acheck2[9];
  for (int j=0;j<9;j++) {
    acheck2[j]=0.;
  }

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
	acheck2[j]+=depositedIonEnergyECAL_pidtime[j][k];
      }
      //      for (int k=15;k<80;k++) {
      for (int k=5;k<80;k++) {
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
    h_nonel3->Fill(depositedEnergyTotal,Ninelastic);
		    
    for (int j=0;j<80;j++) {
      float atime = j*0.25+(0.25/2);

      
      h_timepim->Fill(atime,depositedIonEnergyECAL_pidtime[0][j]);
      h_timeep->Fill(atime,depositedIonEnergyECAL_pidtime[1][j]);
      h_timeem->Fill(atime,depositedIonEnergyECAL_pidtime[2][j]);
      h_timegam->Fill(atime,depositedIonEnergyECAL_pidtime[3][j]);
      h_timepip->Fill(atime,depositedIonEnergyECAL_pidtime[4][j]);
      h_timek->Fill(atime,depositedIonEnergyECAL_pidtime[5][j]);
      h_timen->Fill(atime,depositedIonEnergyECAL_pidtime[6][j]);
      h_timep->Fill(atime,depositedIonEnergyECAL_pidtime[7][j]);
      h_timeo->Fill(atime,depositedIonEnergyECAL_pidtime[8][j]);
    }


  }

  h_timepim->Scale(1./acheck2[0]);
  h_timeep->Scale(1./acheck2[1]);
  h_timeem->Scale(1./acheck2[2]);
  h_timegam->Scale(1./acheck2[3]);
  h_timepip->Scale(1./acheck2[4]);
  h_timek->Scale(1./acheck2[5]);
  h_timen->Scale(1./acheck2[6]);
  h_timep->Scale(1./acheck2[7]);
  h_timeo->Scale(1./acheck2[8]);


  f->Close();
  /*
  for(int ii=0;ii<80;ii++) {
    float aaaa=h_timeep->GetBinContent(ii);
    std::cout<<" bin "<<ii<<" = "<<aaaa<<std::endl;
  }
  */

  TFile * out = new TFile(outputfilename,"RECREATE");
  hTotalE->Write();
  hWorldE->Write();
  h_nonel->Write();
  h_nonel2->Write();
  h_nonel3->Write();
  h_timepim->Write();
  h_timeep->Write();
  h_timeem->Write();
  h_timegam->Write();
  h_timepip->Write();
  h_timek->Write();
  h_timen->Write();
  h_timep->Write();
  h_timeo->Write();
  out->Close();







  h_timepim->SetMarkerColor(kBlack);
  h_timeep->SetMarkerColor(kBlue);
  h_timeem->SetMarkerColor(kRed);
  h_timegam->SetMarkerColor(kGreen);
  h_timepip->SetMarkerColor(kCyan);
  h_timek->SetMarkerColor(kOrange);
  h_timen->SetMarkerColor(kGray);
  h_timep->SetMarkerColor(kPink);
  h_timeo->SetMarkerColor(kYellow);



  h_timepim->SetLineColor(kBlack);
  h_timeep->SetLineColor(kBlue);
  h_timeem->SetLineColor(kRed);
  h_timegam->SetLineColor(kGreen);
  h_timepip->SetLineColor(kCyan);
  h_timek->SetLineColor(kOrange);
  h_timen->SetLineColor(kGray);
  h_timep->SetLineColor(kPink);
  h_timeo->SetLineColor(kYellow);


  h_timepim->SetTitle("pi-");
  h_timeep->SetTitle("e+");
  h_timeem->SetTitle("e-;time (ns); fraction of ionizing energy");
  h_timegam->SetTitle("gamma");
  h_timepip->SetTitle("pi+");
  h_timek->SetTitle("k");
  h_timen->SetTitle("n");
  h_timep->SetTitle("p");
  h_timeo->SetTitle("other");


  TCanvas *c = new TCanvas();
  c->SetLogy(1);

  h_timeem->Draw("HIST");
  h_timepim->Draw("same HIST");
  h_timeep->Draw("same HIST");
  h_timegam->Draw("same HIST");
  h_timepip->Draw("same HIST");
  //h_timek->Draw("same HIST");
  //h_timen->Draw("same HIST");
  h_timep->Draw("same HIST");
  h_timeo->Draw("same HIST");




  c->BuildLegend(0.3,0.6,0.6,0.9);

}


