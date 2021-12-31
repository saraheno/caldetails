#ifndef CreateTree_H
#define CreateTree_H 1

#include <iostream>
#include <vector>
#include "TString.h"
#include <map>

#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

class CreateTree
{
private:
  TTree *ftree;
  TString fname;

public:
  CreateTree(TString name);
  ~CreateTree();

  TTree *GetTree() const { return ftree; };
  TString GetName() const { return fname; };
  void AddEnergyDeposit(int index, float deposit);
  void AddScintillationPhoton(int index);
  void AddCerenkovPhoton(int index);
  int Fill();
  bool Write(TFile *);
  void Clear();

  static CreateTree *Instance() { return fInstance; };
  static CreateTree *fInstance;

  int Event;




  int inputE1Thick;
  int inputE1Width;


  std::vector<float> *inputMomentum;        // Px Py Pz E
  std::vector<float> *inputInitialPosition; // x, y, z

  std::vector<float> *primaryMomT1; // Px Py Pz E
  std::vector<float> *primaryPosT1; // x, y, z

  std::vector<float> *primaryMomE1; // Px Py Pz E
  std::vector<float> *primaryPosE1; // x, y, z

  int nTracksE1;



  //integrated energy in each longitudinal layer
  float depositedEnergyEscapeWorld;
  float depositedTotalEnergyEscapeWorld;

  float depositedEnergyTotal;
  float depositedEnergyECAL_f;
  float depositedEnergyEcalGap;
  float depositedEnergyEcalDet;
  float depositedEnergyWorld;

  float depositedIonEnergyTotal;
  float depositedIonEnergyECAL_f;
  float depositedIonEnergyEcalGap;
  float depositedIonEnergyEcalDet;
  float depositedIonEnergyWorld;

  float depositedElecEnergyTotal;
  float depositedHadronIonEnergyTotal;
  float depositedElecEnergyECAL_f;
  float depositedElecEnergyEcalGap;
  float depositedElecEnergyEcalDet;
  float depositedElecEnergyWorld;

  //store the energy deposition by components

  float depositedEnergyECAL_absorb_f_particleID[8];
  float depositedIonEnergyECAL_absorb_f_particleID[8];
  float betaparticleID[8];
  float Ninelastic;







  int ECAL_f_total_S;
  int ECAL_r_total_S;
  int ECAL_f_total_C;
  int ECAL_r_total_C;

  int tot_phot_cer_ECAL_scinti_r_particleID[8];
  int tot_phot_cer_ECAL_cheren_f_total;
  int tot_phot_cer_ECAL_cheren_r_total;
  int tot_phot_cer_ECAL_cheren_f_particleID[8];
  int tot_phot_cer_ECAL_cheren_r_particleID[8];


  int SDdetected_ff_S;
  int SDdetected_ff_C;
  int SDdetected_rr_S;
  int SDdetected_rr_C;
  /***************** begin to seperate energy into different channels    ******************/


  //energy deposit in each trasnversally segmented channel
  float Edep_Timing_f_ch[18];
  float Edep_Timing_r_ch[18];

  float Edep_ECAL_f_ch[6400];
  float Edep_ECAL_r_ch[6400];

  float IonEdep_ECAL_f_ch[6400];
  float IonEdep_ECAL_r_ch[6400];

  float E_Zdep_0to5000mm_total[2500];
  float E_Tdep_0to5ns_total[2500];
  float E_Zdep_0to5000mm_Pion_n[2500];
  float E_Tdep_0to5ns_Pion_n[2500];
  float E_Zdep_0to5000mm_Positron[2500];
  float E_Tdep_0to5ns_Positron[2500];
  float E_Zdep_0to5000mm_Electron[2500];
  float E_Tdep_0to5ns_Electron[2500];
  float E_Zdep_0to5000mm_Photon[2500];
  float E_Tdep_0to5ns_Photon[2500];
  float E_Zdep_0to5000mm_Pion_p[2500];
  float E_Tdep_0to5ns_Pion_p[2500];
  float E_Zdep_0to5000mm_Kion[2500];
  float E_Tdep_0to5ns_Kion[2500];
  float E_Zdep_0to5000mm_Neutron[2500];
  float E_Tdep_0to5ns_Neutron[2500];
  float E_Zdep_0to5000mm_Proton[2500];
  float E_Tdep_0to5ns_Proton[2500];

  TH1F *h_phot_lambda_ECAL_f_Ceren;
  TH1F *h_phot_lambda_ECAL_r_Ceren;
  TH1F *h_phot_lambda_ECAL_f_produce_Ceren;
  TH1F *h_phot_lambda_ECAL_r_produce_Ceren;
  TH2F *h_photon_2D_produce_Ceren;
  TH2F *h_photon_2D_receive_Ceren;


  TH1F *h_phot_lambda_ECAL_f_Scin;
  TH1F *h_phot_lambda_ECAL_r_Scin;
  TH1F *h_phot_lambda_ECAL_f_produce_Scin;
  TH1F *h_phot_lambda_ECAL_r_produce_Scin;
  TH2F *h_photon_2D_produce_Scin;
  TH2F *h_photon_2D_receive_Scin;
  TH2F *h_photon_hits_map;
  TH2F *ion_rz;
  TH1F *ion_r1;
  TH1F *ion_r2;
  TH1F *ion_r3;
  TH1F *ion_r4;
  TH1F *ion_r5;
  TH1F *ion_r6;
  TH1F *ion_r7;

  TH1F *ion_z;
  TH1F *h_phot_detect_time_f_Scin;
  TH1F *h_phot_detect_time_r_Scin;
  TH1F *h_phot_detect_time_f_Ceren;
  TH1F *h_phot_detect_time_r_Ceren;
  TH2F *h_time_z_egamma;
  TH2F *h_time_z_other;
  TH2F *pdg_beta;
  TH2F *pdg_ke;
  TH1F *inelasticEK;

  TH1F *h_totaldepositedE;
  TH1F *h_totaldepositedEpescapeke;



  std::vector<int> *pdgid_escape;
  std::vector<float> *KineticEnergy_escape;
  std::vector<float> *positionx_escape;
  std::vector<float> *positiony_escape;
  std::vector<float> *positionz_escape;

};


#endif
