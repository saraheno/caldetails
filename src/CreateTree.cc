#include "CreateTree.hh"
#include <algorithm>

using namespace std;

CreateTree *CreateTree::fInstance = NULL;

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

CreateTree::CreateTree(TString name)
{
  if (fInstance)
  {
    return;
  }

  this->fInstance = this;
  this->fname = name;
  this->ftree = new TTree(name, name);

  this->GetTree()->Branch("Event", &this->Event, "Event/I");



  this->GetTree()->Branch("inputE1Thick", &this->inputE1Thick, "inputE1Thick/F");

  this->GetTree()->Branch("inputE1Width", &this->inputE1Width, "inputE1Width/F");


  inputInitialPosition = new vector<float>(3, 0.);
  inputMomentum = new vector<float>(4, 0.);
  primaryPosT1 = new vector<float>(3, 0.);
  primaryMomT1 = new vector<float>(4, 0.);
  primaryPosE1 = new vector<float>(3, 0.);
  primaryMomE1 = new vector<float>(4, 0.);

pdgid_escape  = new vector<int>;
KineticEnergy_escape  = new vector<float>;
positionx_escape = new vector<float>;
positiony_escape = new vector<float>;
positionz_escape = new vector<float>;

  this->GetTree()->Branch("inputInitialPosition", "vector<float>", &inputInitialPosition);
  this->GetTree()->Branch("inputMomentum", "vector<float>", &inputMomentum);
  this->GetTree()->Branch("primaryPosT1", "vector<float>", &primaryPosT1);
  this->GetTree()->Branch("primaryMomT1", "vector<float>", &primaryMomT1);
  this->GetTree()->Branch("primaryPosE1", "vector<float>", &primaryPosE1);
  this->GetTree()->Branch("primaryMomE1", "vector<float>", &primaryMomE1);


  this->GetTree()->Branch("nTracksE1", &this->nTracksE1, "nTracksE1/I");



  //integrated per longitudinal layer
  this->GetTree()->Branch("depositedEnergyTotal", &this->depositedEnergyTotal, "depositedEnergyTotal/F");
  this->GetTree()->Branch("kineticEnergyEscapeWorld", &this->kineticEnergyEscapeWorld, "kineticEnergyEscapeWorld/F");
  this->GetTree()->Branch("depositedEnergyECAL_f", &this->depositedEnergyECAL_f, "depositedEnergyECAL_f/F");
  this->GetTree()->Branch("depositedEnergyWrap", &this->depositedEnergyWrap, "depositedEnergyWrap/F");

  this->GetTree()->Branch("depositedEnergyWorld", &this->depositedEnergyWorld, "depositedEnergyWorld/F");

  this->GetTree()->Branch("depositedIonEnergyTotal", &this->depositedIonEnergyTotal, "depositedIonEnergyTotal/F");


  this->GetTree()->Branch("depositedIonEnergyECAL_f", &this->depositedIonEnergyECAL_f, "depositedIonEnergyECAL_f/F");
  this->GetTree()->Branch("depositedIonEnergyWorld", &this->depositedIonEnergyWorld, "depositedIonEnergyWorld/F");
  this->GetTree()->Branch("depositedIonEnergyWrap", &this->depositedIonEnergyWrap, "depositedIonEnergyWrap/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_absorb_f_particleID", &this->depositedIonEnergyECAL_absorb_f_particleID, "depositedIonEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("betaparticleID", &this->betaparticleID, "betaparticleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_absorb_f_particleID", &this->depositedEnergyECAL_absorb_f_particleID, "depositedEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("Ninelastic", &this->Ninelastic, "Ninelastic/F");
  this->GetTree()->Branch("nNeutrons", &this->nNeutrons, "nNeutrons/F");


  this->GetTree()->Branch("depositedElecEnergyTotal", &this->depositedElecEnergyTotal, "depositedElecEnergyTotal/F");
  this->GetTree()->Branch("depositedHadronIonEnergyTotal", &this->depositedHadronIonEnergyTotal, "depositedHadronIonEnergyTotal/F");
  this->GetTree()->Branch("depositedElecEnergyECAL_f", &this->depositedElecEnergyECAL_f, "depositedElecEnergyECAL_f[3]/F");

  this->GetTree()->Branch("depositedElecEnergyWorld", &this->depositedElecEnergyWorld, "depositedElecEnergyWorld/F");
  this->GetTree()->Branch("depositedElecEnergyWrap", &this->depositedElecEnergyWrap, "depositedElecEnergyWrap/F");


  //Cerenkov photons





  //S and C in crystal bar
  //total generated in front and rear part
  this->GetTree()->Branch("ECAL_f_total_S", &this->ECAL_f_total_S, "ECAL_f_total_S/I");
  this->GetTree()->Branch("ECAL_f_total_C", &this->ECAL_f_total_C, "ECAL_f_total_C/I");
  this->GetTree()->Branch("ECAL_r_total_S", &this->ECAL_r_total_S, "ECAL_r_total_S/I");
  this->GetTree()->Branch("ECAL_r_total_C", &this->ECAL_r_total_C, "ECAL_r_total_C/I");  
  //detected in two detectors
  this->GetTree()->Branch("SDdetected_ff_S", &this->SDdetected_ff_S, "SDdetected_ff_S/I");
  this->GetTree()->Branch("SDdetected_ff_C", &this->SDdetected_ff_C, "SDdetected_ff_C/I");
  this->GetTree()->Branch("SDdetected_rr_S", &this->SDdetected_rr_S, "SDdetected_rr_S/I");
  this->GetTree()->Branch("SDdetected_rr_C", &this->SDdetected_rr_C, "SDdetected_rr_C/I");


  // basic plots
  h_totaliondepositedE = new TH1F("h_totaliondepositedE","",150,0.,1.5);
  h_totaldepositedE = new TH1F("h_totaldepositedE","",150,0.,1.5);
  h_totalpkineticenergyescape = new TH1F("h_totalpkineticenergyescape","",150,0.,1.5);
  pdg_ke = new TH2F("pdg_ke","",10000,-5000,5000,100,0,10.);
  h_nneutrons = new TH1F("h_nneutrons","",2000,0.,10000);
  h_keneutrons = new TH1F("h_keneutrons","",150,0.,0.01);
  h_nonelvpe = new TH2F("h_nonelvpe","",100,0.,5.,600,0.,600.);
  h_nonelvlst = new TH2F("h_nonelvlst","",100,0.5,1.1,600,0.,600.);
  h_nnvlst = new TH2F("h_nnvlst","",100,0.5,1.1,1000,0.,10000.);
  h_nnvnonel = new TH2F("h_nnvnonel","",2000,0,10000,600,0.,600.);
  h_pevlst = new TH2F("h_pevlst","",100,0.5,1.1,100,0.,5.);


  //detected photons 
  h_phot_lambda_ECAL_f_Scin = new TH1F("h_phot_lambda_ECAL_f_Scin", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_r_Scin = new TH1F("h_phot_lambda_ECAL_r_Scin", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_f_Ceren = new TH1F("h_phot_lambda_ECAL_f_Ceren", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_r_Ceren = new TH1F("h_phot_lambda_ECAL_r_Ceren", "", 1250, 0., 1250.);
  //generated photons
  h_phot_lambda_ECAL_f_produce_Scin = new TH1F("h_phot_lambda_ECAL_f_produce_Scin", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_r_produce_Scin = new TH1F("h_phot_lambda_ECAL_r_produce_Scin", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_f_produce_Ceren = new TH1F("h_phot_lambda_ECAL_f_produce_Ceren", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_r_produce_Ceren = new TH1F("h_phot_lambda_ECAL_r_produce_Ceren", "", 1250, 0., 1250.);

  //photons detected time
  ion_z = new TH1F("ion_z","",1500,-1500,1500);
  ion_rz = new TH2F("ion_rz","",1500,0,1000, 15,-1500,1500);
  ion_r1 = new TH1F("ion_r1","",1500,0,1000);
  ion_r2 = new TH1F("ion_r2","",1500,0,1000);
  ion_r3 = new TH1F("ion_r3","",1500,0,1000);
  ion_r4 = new TH1F("ion_r4","",1500,0,1000);
  ion_r5 = new TH1F("ion_r5","",1500,0,1000);
  ion_r6 = new TH1F("ion_r6","",1500,0,1000);
  ion_r7 = new TH1F("ion_r7","",1500,0,1000);

  h_time_z_egamma = new TH2F("h_time_z_egamma","",1250,0.,1500, 1250,0,0.1);
  h_time_z_other = new TH2F("h_time_z_other","",1250,0.,1500,1250,0,0.1);


  //photon position
  h_photon_2D_produce_Scin = new TH2F("h_photon_2D_produce_Scin", "", 500, -50, 50, 500, 0., 200);
  h_photon_2D_receive_Scin = new TH2F("h_photon_2D_receive_Scin", "", 500, -50, 50, 500, 0., 200);
  h_photon_2D_produce_Ceren = new TH2F("h_photon_2D_produce_Ceren", "", 500, -50, 50, 500, 0., 200);
  h_photon_2D_receive_Ceren = new TH2F("h_photon_2D_receive_Ceren", "", 500, -50, 50, 500, 0., 200);
  
  h_photon_hits_map = new TH2F("h_photon_hits_map", "", 100, -5, 5,100, -5, 5);
  pdg_beta = new TH2F("pdg_beta","",10000,-5000,5000,100,0,1);
  inelasticEK = new TH1F("inelasticEK","",1000,0,1);
  this->GetTree()->Branch("pdgid_escape", "vector<int>", &pdgid_escape);
  this->GetTree()->Branch("KineticEnergy_escape", "vector<float>", &KineticEnergy_escape);
  this->GetTree()->Branch("positionx_escape", "vector<float>", &positionx_escape);
  this->GetTree()->Branch("positiony_escape", "vector<float>", &positiony_escape);
  this->GetTree()->Branch("positionz_escape", "vector<float>", &positionz_escape);

  this->Clear();

}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

CreateTree::~CreateTree()
{
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

int CreateTree::Fill()
{
//  this->GetTree()->Write(NULL, TObject::kOverwrite );

  std::cout<<"depositedEnergyTotal is "<<depositedEnergyTotal<<std::endl;
  std::cout<<"kineticEnergyEscapeWorld is "<<kineticEnergyEscapeWorld<<std::endl;
  float sum = depositedEnergyTotal+kineticEnergyEscapeWorld;
  std::cout<<"sum is "<<sum<<std::endl;
  std::cout<<"depositedEnergyECAL_f is "<<depositedEnergyECAL_f<<std::endl;
  std::cout<<"depositedEnergyWrap is "<<depositedEnergyWrap<<std::endl;
  std::cout<<"depositedEnergyworld is "<<depositedEnergyWorld<<std::endl;
  float diff = depositedEnergyTotal - depositedEnergyECAL_f-depositedEnergyWrap-depositedEnergyWorld;
  std::cout<<"diff is "<<diff<<std::endl;


  h_totaliondepositedE->Fill( depositedIonEnergyTotal/(inputMomentum->at(3)) );
  h_totaldepositedE->Fill( depositedEnergyTotal/(inputMomentum->at(3)) );
  h_totalpkineticenergyescape->Fill((depositedEnergyTotal+kineticEnergyEscapeWorld)/(inputMomentum->at(3)));
  h_nneutrons->Fill(nNeutrons);
  h_nonelvpe->Fill(depositedEnergyECAL_absorb_f_particleID[7],Ninelastic);// proton energy deposit versus number of inelastic
  h_nonelvlst->Fill(depositedEnergyECAL_f/(inputMomentum->at(3)),Ninelastic);
  h_nnvlst->Fill(depositedEnergyECAL_f/(inputMomentum->at(3)),nNeutrons);
  h_nnvnonel->Fill(nNeutrons,Ninelastic);
  h_pevlst->Fill(depositedEnergyECAL_f/(inputMomentum->at(3)),depositedEnergyECAL_absorb_f_particleID[7]);


  return this->GetTree()->Fill();
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

bool CreateTree::Write(TFile *outfile)
{
  outfile->cd();
  ftree->Write();
  h_phot_lambda_ECAL_f_Scin->Write();
  h_phot_lambda_ECAL_r_Scin->Write();
  h_phot_lambda_ECAL_f_produce_Scin->Write();
  h_phot_lambda_ECAL_r_produce_Scin->Write();
  h_photon_2D_produce_Scin->Write();
  h_photon_2D_receive_Scin->Write();

  h_phot_lambda_ECAL_f_Ceren->Write();
  h_phot_lambda_ECAL_r_Ceren->Write();
  h_phot_lambda_ECAL_f_produce_Ceren->Write();
  h_phot_lambda_ECAL_r_produce_Ceren->Write();
  h_photon_2D_produce_Ceren->Write();
  h_photon_2D_receive_Ceren->Write();
  h_photon_hits_map->Write();
  pdg_beta->Write();
  inelasticEK->Write();
  ion_z->Write();
  ion_rz->Write();
  ion_r1->Write();
  ion_r2->Write();
  ion_r3->Write();
  ion_r4->Write();
  ion_r5->Write();
  ion_r6->Write();
  ion_r7->Write();
  h_time_z_egamma->Write();
  h_time_z_other->Write();


  pdg_ke->Write();
  h_totaldepositedE->Write();
  h_totaliondepositedE->Write();
  h_totalpkineticenergyescape->Write();
  h_nneutrons->Write();
  h_keneutrons->Write();
  h_nonelvpe->Write();
  h_nonelvlst->Write();
  h_nnvlst->Write();
  h_nnvnonel->Write();
  h_pevlst->Write();

  return true;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void CreateTree::Clear()
{
  Event = 0;


  nTracksE1 = 0;




  kineticEnergyEscapeWorld = 0.;
  TotalEnergyEscapeWorld = 0.;

  depositedEnergyTotal = 0.;

  depositedEnergyECAL_f = 0.;
  depositedEnergyWorld = 0.;
  depositedEnergyWrap = 0.;


  depositedIonEnergyTotal = 0.;
  depositedIonEnergyECAL_f = 0.;
  depositedIonEnergyWorld = 0.;
  depositedIonEnergyWrap = 0.;
  depositedHadronIonEnergyTotal = 0.;

  depositedElecEnergyECAL_f = 0.;
  depositedElecEnergyWorld = 0.;
  depositedElecEnergyWrap = 0.;


  ECAL_f_total_C = 0.;
  ECAL_r_total_S = 0.;
  ECAL_f_total_S = 0.;
  ECAL_r_total_C = 0.;
  SDdetected_ff_S = 0.;
  SDdetected_ff_C = 0.;
  SDdetected_rr_S = 0.;
  SDdetected_rr_C = 0.; 
  Ninelastic=0;
  nNeutrons=0;
  for (int iparticle = 0; iparticle < 8; iparticle++)
  {
    depositedEnergyECAL_absorb_f_particleID[iparticle] = 0.;






    
    depositedIonEnergyECAL_absorb_f_particleID[iparticle] = 0.;
    betaparticleID[iparticle] = 0.;






    tot_phot_cer_ECAL_cheren_r_particleID[iparticle] = 0.;
  }




  for (int iCh = 0; iCh < 6400; iCh++)
  {
    Edep_ECAL_f_ch[iCh] = 0.;
    Edep_ECAL_r_ch[iCh] = 0.;

    IonEdep_ECAL_f_ch[iCh] = 0.;


  }

  for (int iZ = 0; iZ < 2500; iZ++)
  {
    E_Zdep_0to5000mm_total[iZ]=0.;
    E_Zdep_0to5000mm_Pion_n[iZ]=0.;
    E_Zdep_0to5000mm_Positron[iZ]=0.;
    E_Zdep_0to5000mm_Electron[iZ]=0.;
    E_Zdep_0to5000mm_Photon[iZ]=0.;
    E_Zdep_0to5000mm_Pion_p[iZ]=0.;
    E_Zdep_0to5000mm_Kaon[iZ]=0.;
    E_Zdep_0to5000mm_Neutron[iZ]=0.;
    E_Zdep_0to5000mm_Proton[iZ]=0.;
  }

  for (int i = 0; i < 3; ++i)
  {
    inputInitialPosition->at(i) = 0.;
    primaryPosT1->at(i) = 0.;
    primaryPosE1->at(i) = 0.;
  }
  for (int i = 0; i < 4; ++i)
  {
    inputMomentum->at(i) = 0.;
    primaryMomT1->at(i) = 0.;
    primaryMomE1->at(i) = 0.;
  }

  for (int i = 0; i < 3; ++i)
  {
    inputInitialPosition->at(i) = 0.;
    primaryPosT1->at(i) = 0.;
    primaryPosE1->at(i) = 0.;
  }

KineticEnergy_escape->clear();
pdgid_escape->clear();
positionx_escape->clear();
positiony_escape->clear();
positionz_escape->clear();


}
