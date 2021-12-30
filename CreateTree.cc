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

  this->GetTree()->Branch("inputTrackerX0", &this->inputTrackerX0, "inputTrackerX0/F");
  this->GetTree()->Branch("inputServiceAlmm", &this->inputServiceAlmm, "inputServiceAlmm/F");
  this->GetTree()->Branch("inputTimingThick", &this->inputTimingThick, "inputTimingThick/F");
  this->GetTree()->Branch("inputE1Thick", &this->inputE1Thick, "inputE1Thick/F");
  this->GetTree()->Branch("inputE2Thick", &this->inputE2Thick, "inputE2Thick/F");
  this->GetTree()->Branch("inputE1Width", &this->inputE1Width, "inputE1Width/F");
  this->GetTree()->Branch("inputTimingECAL_dist", &this->inputTimingECAL_dist, "inputTimingECAL_dist/F");

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

  this->GetTree()->Branch("nTracksT1", &this->nTracksT1, "nTracksT1/I");
  this->GetTree()->Branch("nTracksT2", &this->nTracksT2, "nTracksT2/I");
  this->GetTree()->Branch("nTracksE1", &this->nTracksE1, "nTracksE1/I");
  this->GetTree()->Branch("nTracksE2", &this->nTracksE2, "nTracksE2/I");
  this->GetTree()->Branch("nTracksTRK", &this->nTracksTRK, "nTracksTRK[6]/F");

  //integrated per longitudinal layer
  this->GetTree()->Branch("depositedEnergyTotal", &this->depositedEnergyTotal, "depositedEnergyTotal/F");
  this->GetTree()->Branch("depositedEnergyEscapeWorld", &this->depositedEnergyEscapeWorld, "depositedEnergyEscapeWorld/F");
  this->GetTree()->Branch("depositedEnergyTiming_f", &this->depositedEnergyTiming_f, "depositedEnergyTiming_f/F");
  this->GetTree()->Branch("depositedEnergyTiming_r", &this->depositedEnergyTiming_r, "depositedEnergyTiming_r/F");
  this->GetTree()->Branch("depositedEnergyECAL_f", &this->depositedEnergyECAL_f, "depositedEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedEnergyECAL_r", &this->depositedEnergyECAL_r, "depositedEnergyECAL_r[3]/F");
  this->GetTree()->Branch("depositedEnergyHCALAct", &this->depositedEnergyHCALAct, "depositedEnergyHCALAct/F");
  this->GetTree()->Branch("depositedEnergyHCALPas", &this->depositedEnergyHCALPas, "depositedEnergyHCALPas/F");
  this->GetTree()->Branch("depositedEnergyWorld", &this->depositedEnergyWorld, "depositedEnergyWorld/F");
  this->GetTree()->Branch("depositedEnergyTimingGap", &this->depositedEnergyTimingGap, "depositedEnergyTimingGap/F");
  this->GetTree()->Branch("depositedEnergyServices", &this->depositedEnergyServices, "depositedEnergyServices/F");
  this->GetTree()->Branch("depositedEnergyEcalGap", &this->depositedEnergyEcalGap, "depositedEnergyEcalGap/F");
  this->GetTree()->Branch("depositedEnergyEcalDet", &this->depositedEnergyEcalDet, "depositedEnergyEcalDet/F");
  this->GetTree()->Branch("depositedEnergySolenoid", &this->depositedEnergySolenoid, "depositedEnergySolenoid/F");
/*
  this->GetTree()->Branch("depositedEnergyECAL_absorb_f_particleID", &this->depositedEnergyECAL_absorb_f_particleID, "depositedEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_absorb_r_particleID", &this->depositedEnergyECAL_absorb_r_particleID, "depositedEnergyECAL_absorb_r_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_scinti_f_particleID", &this->depositedEnergyECAL_scinti_f_particleID, "depositedEnergyECAL_scinti_f_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_scinti_r_particleID", &this->depositedEnergyECAL_scinti_r_particleID, "depositedEnergyECAL_scinti_r_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_cheren_f_particleID", &this->depositedEnergyECAL_cheren_f_particleID, "depositedEnergyECAL_cheren_f_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_cheren_r_particleID", &this->depositedEnergyECAL_cheren_r_particleID, "depositedEnergyECAL_cheren_r_particleID[8]/F");
*/
  this->GetTree()->Branch("depositedIonEnergyTotal", &this->depositedIonEnergyTotal, "depositedIonEnergyTotal/F");
  this->GetTree()->Branch("depositedIonEnergyTiming_f", &this->depositedIonEnergyTiming_f, "depositedIonEnergyTiming_f/F");
  this->GetTree()->Branch("depositedIonEnergyTiming_r", &this->depositedIonEnergyTiming_r, "depositedIonEnergyTiming_r/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_f", &this->depositedIonEnergyECAL_f, "depositedIonEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_r", &this->depositedIonEnergyECAL_r, "depositedIonEnergyECAL_r[3]/F");
  this->GetTree()->Branch("depositedIonEnergyHCALAct", &this->depositedIonEnergyHCALAct, "depositedIonEnergyHCALAct/F");
  this->GetTree()->Branch("depositedIonEnergyHCALPas", &this->depositedIonEnergyHCALPas, "depositedIonEnergyHCALPas/F");
  this->GetTree()->Branch("depositedIonEnergyWorld", &this->depositedIonEnergyWorld, "depositedIonEnergyWorld/F");
  this->GetTree()->Branch("depositedIonEnergyTimingGap", &this->depositedIonEnergyTimingGap, "depositedIonEnergyTimingGap/F");
  this->GetTree()->Branch("depositedIonEnergyServices", &this->depositedIonEnergyServices, "depositedIonEnergyServices/F");
  this->GetTree()->Branch("depositedIonEnergyEcalGap", &this->depositedIonEnergyEcalGap, "depositedIonEnergyEcalGap/F");
  this->GetTree()->Branch("depositedIonEnergyEcalDet", &this->depositedIonEnergyEcalDet, "depositedIonEnergyEcalDet/F");
  this->GetTree()->Branch("depositedIonEnergySolenoid", &this->depositedIonEnergySolenoid, "depositedIonEnergySolenoid/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_absorb_f_particleID", &this->depositedIonEnergyECAL_absorb_f_particleID, "depositedIonEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("betaparticleID", &this->betaparticleID, "betaparticleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_absorb_f_particleID", &this->depositedEnergyECAL_absorb_f_particleID, "depositedEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("Ninelastic", &this->Ninelastic, "Ninelastic/F");

/*
  this->GetTree()->Branch("depositedIonEnergyECAL_absorb_f_particleID", &this->depositedIonEnergyECAL_absorb_f_particleID, "depositedIonEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("betaparticleID", &this->betaparticleID, "betaparticleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_scinti_f_particleID", &this->depositedIonEnergyECAL_scinti_f_particleID, "depositedIonEnergyECAL_scinti_f_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_scinti_r_particleID", &this->depositedIonEnergyECAL_scinti_r_particleID, "depositedIonEnergyECAL_scinti_r_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_cheren_f_particleID", &this->depositedIonEnergyECAL_cheren_f_particleID, "depositedIonEnergyECAL_cheren_f_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_cheren_r_particleID", &this->depositedIonEnergyECAL_cheren_r_particleID, "depositedIonEnergyECAL_cheren_r_particleID[8]/F");
*/
  this->GetTree()->Branch("depositedElecEnergyTotal", &this->depositedElecEnergyTotal, "depositedElecEnergyTotal/F");
  this->GetTree()->Branch("depositedHadronIonEnergyTotal", &this->depositedHadronIonEnergyTotal, "depositedHadronIonEnergyTotal/F");
  this->GetTree()->Branch("depositedElecEnergyTiming_f", &this->depositedElecEnergyTiming_f, "depositedElecEnergyTiming_f/F");
  this->GetTree()->Branch("depositedElecEnergyTiming_r", &this->depositedElecEnergyTiming_r, "depositedElecEnergyTiming_r/F");
  this->GetTree()->Branch("depositedElecEnergyECAL_f", &this->depositedElecEnergyECAL_f, "depositedElecEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedElecEnergyECAL_r", &this->depositedElecEnergyECAL_r, "depositedElecEnergyECAL_r[3]/F");
  this->GetTree()->Branch("depositedElecEnergyHCALAct", &this->depositedElecEnergyHCALAct, "depositedElecEnergyHCALAct/F");
  this->GetTree()->Branch("depositedElecEnergyHCALPas", &this->depositedElecEnergyHCALPas, "depositedElecEnergyHCALPas/F");
  this->GetTree()->Branch("depositedElecEnergyWorld", &this->depositedElecEnergyWorld, "depositedElecEnergyWorld/F");
  this->GetTree()->Branch("depositedElecEnergyTimingGap", &this->depositedElecEnergyTimingGap, "depositedElecEnergyTimingGap/F");
  this->GetTree()->Branch("depositedElecEnergyServices", &this->depositedElecEnergyServices, "depositedElecEnergyServices/F");
  this->GetTree()->Branch("depositedElecEnergyEcalGap", &this->depositedElecEnergyEcalGap, "depositedElecEnergyEcalGap/F");
  this->GetTree()->Branch("depositedElecEnergyEcalDet", &this->depositedElecEnergyEcalDet, "depositedElecEnergyEcalDet/F");
  this->GetTree()->Branch("depositedElecEnergySolenoid", &this->depositedElecEnergySolenoid, "depositedElecEnergySolenoid/F");

  //single channels
/*
  this->GetTree()->Branch("Edep_Tracker_layer", &this->Edep_Tracker_layer, "Edep_Tracker_layer[6]/F");
  this->GetTree()->Branch("Edep_Timing_f_ch", &this->Edep_Timing_f_ch, "Edep_Timing_f_ch[18]/F");
  this->GetTree()->Branch("Edep_Timing_r_ch", &this->Edep_Timing_r_ch, "Edep_Timing_r_ch[18]/F");
  this->GetTree()->Branch("Edep_ECAL_f_ch", &this->Edep_ECAL_f_ch, "Edep_ECAL_f_ch[6400]/F");
  this->GetTree()->Branch("Edep_ECAL_r_ch", &this->Edep_ECAL_r_ch, "Edep_ECAL_r_ch[6400]/F");

  this->GetTree()->Branch("IonEdep_ECAL_f_ch", &this->IonEdep_ECAL_f_ch, "IonEdep_ECAL_f_ch[6400]/F");
  this->GetTree()->Branch("IonEdep_ECAL_r_ch", &this->IonEdep_ECAL_r_ch, "IonEdep_ECAL_r_ch[6400]/F");
*/
  //Cerenkov photons
  //  this -> GetTree() -> Branch("tot_phot_sci_Timing",        &this->tot_phot_sci_Timing,               "tot_phot_sci_Timing/I");
  this->GetTree()->Branch("tot_phot_cer_Timing_f_total", &this->tot_phot_cer_Timing_f_total, "tot_phot_cer_Timing_f_total/I");
  this->GetTree()->Branch("tot_phot_cer_Timing_r_total", &this->tot_phot_cer_Timing_r_total, "tot_phot_cer_Timing_r_total/I");
/*
  this->GetTree()->Branch("tot_phot_cer_ECAL_scinti_f_particleID", &this->tot_phot_cer_ECAL_scinti_f_particleID, "tot_phot_cer_ECAL_scinti_f_particleID[8]/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_scinti_r_particleID", &this->tot_phot_cer_ECAL_scinti_r_particleID, "tot_phot_cer_ECAL_scinti_r_particleID[8]/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_cheren_f_particleID", &this->tot_phot_cer_ECAL_cheren_f_particleID, "tot_phot_cer_ECAL_cheren_f_particleID[8]/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_cheren_r_particleID", &this->tot_phot_cer_ECAL_cheren_r_particleID, "tot_phot_cer_ECAL_cheren_r_particleID[8]/I");
  this->GetTree()->Branch("tot_phot_cer_HCAL", &this->tot_phot_cer_HCAL, "tot_phot_cer_HCAL/I");
  this->GetTree()->Branch("E_Zdep_0to5000mm_total", &this->E_Zdep_0to5000mm_total, "E_Zdep_0to5000mm_total[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Pion_n", &this->E_Zdep_0to5000mm_Pion_n, "E_Zdep_0to5000mm_Pion_n[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Positron", &this->E_Zdep_0to5000mm_Positron, "E_Zdep_0to5000mm_Positron[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Electron", &this->E_Zdep_0to5000mm_Electron, "E_Zdep_0to5000mm_Electron[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Photon", &this->E_Zdep_0to5000mm_Photon, "E_Zdep_0to5000mm_Photon[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Pion_p", &this->E_Zdep_0to5000mm_Pion_p, "E_Zdep_0to5000mm_Pion_p[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Kion", &this->E_Zdep_0to5000mm_Kion, "E_Zdep_0to5000mm_Kion[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Neutron", &this->E_Zdep_0to5000mm_Neutron, "E_Zdep_0to5000mm_Neutron[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Proton", &this->E_Zdep_0to5000mm_Proton, "E_Zdep_0to5000mm_Proton[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_total", &this->E_Tdep_0to5ns_total, "E_Tdep_0to5ns_total[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Pion_n", &this->E_Tdep_0to5ns_Pion_n, "E_Tdep_0to5ns_Pion_n[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Positron", &this->E_Tdep_0to5ns_Positron, "E_Tdep_0to5ns_Positron[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Electron", &this->E_Tdep_0to5ns_Electron, "E_Tdep_0to5ns_Electron[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Photon", &this->E_Tdep_0to5ns_Photon, "E_Tdep_0to5ns_Photon[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Pion_p", &this->E_Tdep_0to5ns_Pion_p, "E_Tdep_0to5ns_Pion_p[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Kion", &this->E_Tdep_0to5ns_Kion, "E_Tdep_0to5ns_Kion[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Neutron", &this->E_Tdep_0to5ns_Neutron, "E_Tdep_0to5ns_Neutron[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Proton", &this->E_Tdep_0to5ns_Proton, "E_Tdep_0to5ns_Proton[2500]/F");
*/



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

  h_time_z_egamma = new TH2F("h_time_z_egamma","",1250,-1500,1500, 1250,0,50);
  h_time_z_other = new TH2F("h_time_z_other","",1250,-1500,1500,1250,0,50);

  h_phot_detect_time_f_Scin = new TH1F("h_phot_detect_time_f_Scin","",1250,0,25);
  h_phot_detect_time_r_Scin = new TH1F("h_phot_detect_time_r_Scin","",1250,0,25);
  h_phot_detect_time_f_Ceren = new TH1F("h_phot_detect_time_f_Ceren","",1250,0,25);
  h_phot_detect_time_r_Ceren  = new TH1F("h_phot_detect_time_r_Ceren","",1250,0,25);
  h_totaldepositedE = new TH1F("h_totaldepositedE","",100,0.,10.);
  h_totaldepositedEpescapeke = new TH1F("h_totaldepositedEpescapeke","",100,0.,10.);

  //photon position
  h_photon_2D_produce_Scin = new TH2F("h_photon_2D_produce_Scin", "", 500, -50, 50, 500, 0., 200);
  h_photon_2D_receive_Scin = new TH2F("h_photon_2D_receive_Scin", "", 500, -50, 50, 500, 0., 200);
  h_photon_2D_produce_Ceren = new TH2F("h_photon_2D_produce_Ceren", "", 500, -50, 50, 500, 0., 200);
  h_photon_2D_receive_Ceren = new TH2F("h_photon_2D_receive_Ceren", "", 500, -50, 50, 500, 0., 200);
  
  h_photon_hits_map = new TH2F("h_photon_hits_map", "", 100, -5, 5,100, -5, 5);
  pdg_beta = new TH2F("pdg_beta","",10000,-5000,5000,100,0,1);
  pdg_ke = new TH2F("pdg_ke","",10000,-5000,5000,100,0,10.);
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
  std::cout<<"depositedEnergyEscapeWorld is "<<depositedEnergyEscapeWorld<<std::endl;
  float sum = depositedEnergyTotal+depositedEnergyEscapeWorld;
  std::cout<<"sum is "<<sum<<std::endl;
  std::cout<<"depositedEnergyEcalDet is "<<depositedEnergyECAL_f[0]<<std::endl;
  float diff = depositedEnergyTotal - depositedEnergyECAL_f[0];
  std::cout<<"diff is "<<diff<<std::endl;
 
  h_totaldepositedE->Fill(depositedEnergyTotal);
  h_totaldepositedEpescapeke->Fill(depositedEnergyTotal+depositedEnergyEscapeWorld);



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
  pdg_ke->Write();
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
  h_totaldepositedE->Write();
  h_totaldepositedEpescapeke->Write();
  h_phot_detect_time_f_Scin->Write();
  h_phot_detect_time_r_Scin->Write();
  h_phot_detect_time_f_Ceren->Write();
  h_phot_detect_time_r_Ceren->Write();
  h_time_z_egamma->Write();
  h_time_z_other->Write();
  return true;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void CreateTree::Clear()
{
  Event = 0;

  nTracksT1 = 0;
  nTracksT2 = 0;
  nTracksE1 = 0;
  nTracksE2 = 0;

  for (int iLayer = 0; iLayer < 6; iLayer++)
  {
    nTracksTRK[iLayer] = 0;
  }

  depositedEnergyEscapeWorld = 0.;
  depositedTotalEnergyEscapeWorld = 0.;

  depositedEnergyTotal = 0.;
  depositedEnergyTiming_f = 0.;
  depositedEnergyTiming_r = 0.;
  for (int i = 0; i < 3; i++){
    depositedEnergyECAL_f[i] = 0.;
    depositedEnergyECAL_r[i] = 0.;
  }
  depositedEnergyHCALAct = 0.;
  depositedEnergyHCALPas = 0.;
  depositedEnergyWorld = 0.;
  depositedEnergyServices = 0.;
  depositedEnergyTimingGap = 0.;
  depositedEnergyEcalGap = 0.;
  depositedEnergyEcalDet = 0.;
  depositedEnergySolenoid = 0.;

  depositedIonEnergyTotal = 0.;
  depositedIonEnergyTiming_f = 0.;
  depositedIonEnergyTiming_r = 0.;
  for (int i = 0; i < 3; i++){
    depositedIonEnergyECAL_f[i] = 0.;
    depositedIonEnergyECAL_r[i] = 0.;
  }
  depositedIonEnergyHCALAct = 0.;
  depositedIonEnergyHCALPas = 0.;
  depositedIonEnergyWorld = 0.;
  depositedIonEnergyServices = 0.;
  depositedIonEnergyTimingGap = 0.;
  depositedIonEnergyEcalGap = 0.;
  depositedIonEnergyEcalDet = 0.;
  depositedIonEnergySolenoid = 0.;

  depositedElecEnergyTotal = 0.;
  depositedHadronIonEnergyTotal = 0.;
  depositedElecEnergyTiming_f = 0.;
  depositedElecEnergyTiming_r = 0.;
  for (int i = 0; i < 3; i++){
    depositedElecEnergyECAL_f[i] = 0.;
    depositedElecEnergyECAL_r[i] = 0.;
  }
  depositedElecEnergyHCALAct = 0.;
  depositedElecEnergyHCALPas = 0.;
  depositedElecEnergyWorld = 0.;
  depositedElecEnergyServices = 0.;
  depositedElecEnergyTimingGap = 0.;
  depositedElecEnergyEcalGap = 0.;
  depositedElecEnergyEcalDet = 0.;
  depositedElecEnergySolenoid = 0.;

  tot_phot_cer_Timing_f_total = 0.;
  tot_phot_cer_Timing_r_total = 0.;
  ECAL_f_total_C = 0.;
  ECAL_r_total_S = 0.;
  ECAL_f_total_S = 0.;
  ECAL_r_total_C = 0.;
  tot_phot_cer_HCAL = 0.;
  SDdetected_ff_S = 0.;
  SDdetected_ff_C = 0.;
  SDdetected_rr_S = 0.;
  SDdetected_rr_C = 0.; 
  Ninelastic=0;
  for (int iparticle = 0; iparticle < 8; iparticle++)
  {
    depositedEnergyECAL_absorb_f_particleID[iparticle] = 0.;
    depositedEnergyECAL_absorb_r_particleID[iparticle] = 0.;
    depositedEnergyECAL_scinti_f_particleID[iparticle] = 0.;
    depositedEnergyECAL_scinti_r_particleID[iparticle] = 0.;
    depositedEnergyECAL_cheren_f_particleID[iparticle] = 0.;
    depositedEnergyECAL_cheren_r_particleID[iparticle] = 0.;

    
    depositedIonEnergyECAL_absorb_f_particleID[iparticle] = 0.;
    betaparticleID[iparticle] = 0.;
    depositedIonEnergyECAL_scinti_f_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_scinti_r_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_cheren_f_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_cheren_r_particleID[iparticle] = 0.;

    tot_phot_cer_ECAL_scinti_f_particleID[iparticle] = 0.;
    tot_phot_cer_ECAL_scinti_r_particleID[iparticle] = 0.;
    tot_phot_cer_ECAL_cheren_f_particleID[iparticle] = 0.;
    tot_phot_cer_ECAL_cheren_r_particleID[iparticle] = 0.;
  }

  for (int iLayer = 0; iLayer < 6; iLayer++)
  {
    Edep_Tracker_layer[iLayer] = 0.;
  }

  for (int iBar = 0; iBar < 18; iBar++)
  {
    Edep_Timing_f_ch[iBar] = 0.;
    Edep_Timing_r_ch[iBar] = 0.;
  }
  for (int iCh = 0; iCh < 6400; iCh++)
  {
    Edep_ECAL_f_ch[iCh] = 0.;
    Edep_ECAL_r_ch[iCh] = 0.;

    IonEdep_ECAL_f_ch[iCh] = 0.;
    IonEdep_ECAL_r_ch[iCh] = 0.;

  }

  for (int iZ = 0; iZ < 2500; iZ++)
  {
    E_Zdep_0to5000mm_total[iZ]=0.;
    E_Zdep_0to5000mm_Pion_n[iZ]=0.;
    E_Zdep_0to5000mm_Positron[iZ]=0.;
    E_Zdep_0to5000mm_Electron[iZ]=0.;
    E_Zdep_0to5000mm_Photon[iZ]=0.;
    E_Zdep_0to5000mm_Pion_p[iZ]=0.;
    E_Zdep_0to5000mm_Kion[iZ]=0.;
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