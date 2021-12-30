


void plot_beta() {

    TFile* ff=new TFile("data/pion_10GeV_N50_PWO_0degree_3000mm_A1000mm_sipm_1000mm_small_world_pdg_beta_ionenergy_step.root");
   Float_t         depositedIonEnergyTotal;
   Float_t         depositedIonEnergyECAL_f[3];
   TBranch        *b_depositedIonEnergyTotal;   //!
   TBranch        *b_depositedIonEnergyECAL_f;   //!
   TTree *tree;
   ff->GetObject("tree",tree);
   tree->SetBranchAddress("depositedIonEnergyTotal", &depositedIonEnergyTotal, &b_depositedIonEnergyTotal);
   tree->SetBranchAddress("depositedIonEnergyECAL_f", depositedIonEnergyECAL_f, &b_depositedIonEnergyECAL_f);
   double total_ion=0;
   for (Long64_t jentry=0; jentry<tree->GetEntries();jentry++) {
      tree->GetEntry(jentry);
      total_ion+=depositedIonEnergyTotal;
      //total_ion+=depositedIonEnergyECAL_f[0];
   }
   cout<<total_ion<<endl;

TH1D* pip = ((TH2D*) ff->Get("pdg_beta"))->ProjectionY("",5211,5212);
TH1D *hpip = (TH1D*) pip->Clone();
TH1D* pim = ((TH2D*) ff->Get("pdg_beta"))->ProjectionY("",4789,4791);
TH1D *hpim = (TH1D*) pim->Clone();
TH1D* pi0 = ((TH2D*) ff->Get("pdg_beta"))->ProjectionY("",5111,5113);
TH1D *hpi0 = (TH1D*) pi0->Clone();
TH1D* ep = ((TH2D*) ff->Get("pdg_beta"))->ProjectionY("",4990,4990);
TH1D *hep = (TH1D*) ep->Clone();

TH1D* em = ((TH2D*) ff->Get("pdg_beta"))->ProjectionY("",5012,5012);
TH1D *hem = (TH1D*) em->Clone();
TH1D* p = ((TH2D*) ff->Get("pdg_beta"))->ProjectionY("",7212,7214);
TH1D *hp = (TH1D*) p->Clone();
TH1D* n = ((TH2D*) ff->Get("pdg_beta"))->ProjectionY("",7112,7114);
TH1D *hn = (TH1D*) n->Clone();


TCanvas*c=new TCanvas();
hpip->Draw("HIST");
hpim->Draw("same HIST");
hpi0->Draw("same HIST");
hep->Draw("same HIST");
hem->Draw("same HIST");
hp->Draw("same HIST");
hn->Draw("same HIST");

hpip->SetTitle("pi+;#beta;fraction of ionizing energy");
hpim->SetTitle("pi-");
hpi0->SetTitle("pi0");
hep->SetTitle("e+");
hem->SetTitle("e-");
hp->SetTitle("proton");
hn->SetTitle("neutron");

   hpip->SetLineColor(kBlack);
   hpim->SetLineColor(kBlue);
   hpi0->SetLineColor(kRed);
   hep->SetLineColor(kGreen);
   hem->SetLineColor(kCyan);
   hp->SetLineColor(kOrange);
   hn->SetLineColor(kGray);

   hpip->SetLineWidth(3);
   hpim->SetLineWidth(3);
   hpi0->SetLineWidth(3);
   hep->SetLineWidth(3);
   hem->SetLineWidth(3);
   hp->SetLineWidth(3);
   hn->SetLineWidth(3);

   hpip->Scale(1.0/total_ion);
   hpim->Scale(1.0/total_ion);
   hpi0->Scale(1.0/total_ion);
   hep->Scale(1.0/total_ion);
   hem->Scale(1.0/total_ion);
   hp->Scale(1.0/total_ion);
   hn->Scale(1.0/total_ion);

cout<<hpip->Integral()+hpim->Integral()+hpi0->Integral()+hep->Integral()+hem->Integral()+hp->Integral()+hn->Integral()<<endl;
cout<<((TH2D*)ff->Get("pdg_beta"))->Integral(0,10000,0,101)<<endl;

/*
   hpip->Scale(1.0/hpip->Integral());
   hpim->Scale(1.0/hpim->Integral());
   hpi0->Scale(1.0/hpi0->Integral());
   hep->Scale(1.0/hep->Integral());
   hem->Scale(1.0/hem->Integral());
   hp->Scale(1.0/hp->Integral());
   hn->Scale(1.0/hn->Integral());
*/
c->BuildLegend();

}






