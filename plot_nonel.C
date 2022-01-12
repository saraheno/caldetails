#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"



void plot_nonel(char* input) {

  bool idebug = 1;

  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  myStyle->SetOptStat(0);

  TFile* ff=new TFile(input);
  TH2D* haha = ( (TH2D*) ff->Get("h_nonel2"));



  haha->SetMarkerColor(kBlack);


  haha->SetTitle("nonel; energy deposited after 1.25 ns; number of inelastic collisions");


  TCanvas *c = new TCanvas();


  haha->Draw("HIST");


}


