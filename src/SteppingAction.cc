#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "DR_PMTSD.hh"
#include "DetectorConstruction.hh"
#include "TString.h"
#include "TRandom3.h"
//#include "TCint.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include <time.h>

#include "G4EventManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"

//long int CreateSeed();

using namespace std;
using namespace CLHEP;

int to_int(string name)
{
  int Result; // int which will contain the result
  stringstream convert(name);
  string dummy;
  convert >> dummy;
  convert >> Result;
  return Result;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::SteppingAction(DetectorConstruction *detectorConstruction,
                               const G4int &scint, const G4int &cher) : fDetectorConstruction(detectorConstruction),
                                                                        propagateScintillation(scint),
                                                                        propagateCerenkov(cher)
{
  maxtracklength = 500000. * mm;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::~SteppingAction()
{
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void SteppingAction::UserSteppingAction(const G4Step *theStep)
{

  G4Track *theTrack = theStep->GetTrack();

  //  const G4ThreeVector& theTrackDirection = theTrack->GetMomentumDirection();
  //  const G4ThreeVector& theTrackVertexDirection = theTrack->GetVertexMomentumDirection();
  //  TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformaation());

  G4ParticleDefinition *particleType = theTrack->GetDefinition();
  //G4int trackID = theTrack->GetTrackID();

  G4StepPoint *thePrePoint = theStep->GetPreStepPoint();
  G4StepPoint *thePostPoint = theStep->GetPostStepPoint();
  const G4ThreeVector &thePrePosition = thePrePoint->GetPosition();
  const G4ThreeVector &thePostPosition = thePostPoint->GetPosition();
  G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
  G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
  G4String thePrePVName = "";
  if (thePrePV)
    thePrePVName = thePrePV->GetName();
  G4String thePostPVName = "";
  if (thePostPV)
    thePostPVName = thePostPV->GetName();

  //  G4VSolid* thePreS = thePrePV->GetLogicalVolume()->GetSolid();

  G4int nStep = theTrack->GetCurrentStepNumber();
  G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();

  //        cout << " step length = " << theStep->GetStepLength() << endl;
  //-------------

  // get position
  G4double global_x = thePrePosition.x() / mm;
  G4double global_y = thePrePosition.y() / mm;
  G4double global_z = thePrePosition.z() / mm;

  //if(thePrePVName.contains("ecalDetP_rr")) cout<<"z "<<global_z<<endl;

  G4double energy = theStep->GetTotalEnergyDeposit();
  G4double energyIon = 0;

  double weight = 1.;
  double birkCut = 0.1;
  double charge = thePrePoint->GetCharge();
  if (charge != 0. && theStep->GetStepLength() > 0)
  {
    const G4Material *matbirk = thePrePoint->GetMaterial();
    double density = matbirk->GetDensity();
    double dedx = theStep->GetTotalEnergyDeposit() / theStep->GetStepLength();
    double rkb = 0.03333 * g / (MeV * cm2) / density;
    if (dedx > 0)
    {
      weight = 1. - 0.253694 * log(rkb * dedx);
      if (weight < birkCut)
        weight = birkCut;
      else if (weight > 1.)
        weight = 1.;
    }
  }
  bool use_birk = false;
  bool time_cut = false;
  if (use_birk)
  {
    energyIon = (energy - theStep->GetNonIonizingEnergyDeposit()) * weight;
  }
  else
  {
    energyIon = (energy - theStep->GetNonIonizingEnergyDeposit());
  }
  if (time_cut)
  {
    if (((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3.) > 75)
      energyIon = 0;
  }
  //------------- categorize the particles generated -------------
  G4double energyElec = 0.;
  //total energy by particle types
  G4double energyPion_n = 0.;
  G4double energyPositron = 0.;
  G4double energyElectron = 0.;
  G4double energyPhoton = 0.;
  G4double energyPion_p = 0.;
  G4double energyKaon = 0.;
  G4double energyNeutron = 0.;
  G4double energyProton = 0.;
  //ion energy by particle types
  G4double energyIonPion_n = 0.;
  G4double energyIonPositron = 0.;
  G4double energyIonElectron = 0.;
  G4double energyIonPhoton = 0.;
  G4double energyIonPion_p = 0.;
  G4double energyIonKion = 0.;
  G4double energyIonNeutron = 0.;
  G4double energyIonProton = 0.;

  if (TrPDGid == (-211))
  {
    energyPion_n = energy;
    energyIonPion_n = energyIon;
  }
  else if (TrPDGid == (-11))
  {
    energyPositron = energy;
    energyIonPositron = energyIon;
  }
  else if (TrPDGid == (11))
  {
    energyElectron = energy;
    energyIonElectron = energyIon;
  }
  else if (TrPDGid == (22))
  {
    energyPhoton = energy;
    energyIonPhoton = energyIon;
  }
  else if (TrPDGid == (211))
  {
    energyPion_p = energy;
    energyIonPion_p = energyIon;
  }
  else if (TrPDGid == (321))
  {
    energyKaon = energy;
    energyIonKion = energyIon;
  }
  else if (TrPDGid == (2112))
  {
    energyNeutron = energy;
    energyIonNeutron = energyIon;
  }
  else if (TrPDGid == (2212))
  {
    energyProton = energy;
    energyIonProton = energyIon;
  }
  energyElec = energyIonPositron + energyIonElectron + energyIonPhoton;

  //std::cout<<"TrPDGid energy energyIon enegyElec are "<<TrPDGid<<" "<<energy<<" "<<energyIon<<" "<<energyElec<<std::endl;

  CreateTree::Instance()->depositedEnergyTotal += theStep->GetTotalEnergyDeposit() / GeV;
  CreateTree::Instance()->depositedIonEnergyTotal += energyIon / GeV;
  CreateTree::Instance()->depositedElecEnergyTotal += energyElec / GeV;

  if(theTrack->GetVelocity()/299.792 > 1.0/2.2 && abs(TrPDGid)!=11&&abs(TrPDGid)!=22){
      CreateTree::Instance()->depositedHadronIonEnergyTotal += energyIon / GeV;
      //std::cout<<theTrack->GetDefinition()->GetParticleName()<<" "<<theTrack->GetVelocity()/299.792<<" "<<energyIon / GeV<<std::endl;
  }

//if(nStep == 1) CreateTree::Instance()->pdg_beta->Fill(TrPDGid,theTrack->GetVelocity()/299.792,);
CreateTree::Instance()->pdg_beta->Fill(TrPDGid,thePrePoint->GetBeta(),energyIon / GeV);
 CreateTree::Instance()->pdg_ke->Fill(TrPDGid,thePrePoint->GetKineticEnergy(),energyIon / GeV);

//if(TrPDGid>=111) cout<<energyIon / GeV<<endl;

  if(thePostPoint->GetProcessDefinedStep()->GetProcessName().contains("Inelast")){
    CreateTree::Instance()->Ninelastic++;
    CreateTree::Instance()->inelasticEK->Fill((theStep->GetPreStepPoint()->GetKineticEnergy() - theStep->GetPostStepPoint()->GetKineticEnergy())/GeV);

/*
    cout<<"pdg:"<<TrPDGid<<" process "<<thePostPoint->GetProcessDefinedStep()->GetProcessName()<<"nonion enegry "<<theStep->GetNonIonizingEnergyDeposit()/GeV <<" tot "<<theStep->GetTotalEnergyDeposit() / GeV<<" EK change "<<(theStep->GetPreStepPoint()->GetKineticEnergy() - theStep->GetPostStepPoint()->GetKineticEnergy())/GeV<<endl;
        const std::vector<const G4Track*>* secondaries = theStep->GetSecondaryInCurrentStep();
        for (auto t = secondaries->begin(); t != secondaries->end(); t++) {
            std::cout << "**"<<(*t)->GetParticleDefinition()->GetParticleName()
                      << " ParentID: "<<(*t)->GetParentID()
                      << " IsBelowThreshold: "<<(*t)->IsBelowThreshold()
                      <<std::endl;
            std::cout << " Etot: "<<(*t)->GetTotalEnergy()/GeV<<" Ek: "<<(*t)->GetKineticEnergy()/GeV//<<" p: "<<(*t)->GetMomentum()
                      <<" v: "<<(*t)->GetVelocity()/299.792<<"c "<<" m: "<<  (*t)->GetParticleDefinition()->GetPDGMass()/GeV<<" GeV"
                      <<std::endl;
        }
*/
  }

//  if(thePostPoint->GetProcessDefinedStep()->GetProcessName().contains("Inelast")) cout<<"pdg:"<<TrPDGid<<" process "<<thePostPoint->GetProcessDefinedStep()->GetProcessName()<<"ion enegry"<<energyIon / GeV<<" tot "<<theStep->GetTotalEnergyDeposit() / GeV<<endl;

/*
if (theTrack->GetParentID() == 0 &&  ((thePostPoint->GetProcessDefinedStep()->GetProcessName() == "pi-Inelastic") || (thePostPoint->GetProcessDefinedStep()->GetProcessName() == "pi+Inelastic"))) {

        const std::vector<const G4Track*>* secondaries = theStep->GetSecondaryInCurrentStep();
        if(secondaries->size()==0) std::cout<<"No secondaries"<<std::endl;
        for (auto t = secondaries->begin(); t != secondaries->end(); t++) {
            std::cout << "**"<<(*t)->GetParticleDefinition()->GetParticleName()
                      << " ParentID: "<<(*t)->GetParentID()
                      << " IsBelowThreshold: "<<(*t)->IsBelowThreshold()
                      <<std::endl;
            std::cout << " Etot: "<<(*t)->GetTotalEnergy()/GeV<<" Ek: "<<(*t)->GetKineticEnergy()/GeV//<<" p: "<<(*t)->GetMomentum()
                      <<" v: "<<(*t)->GetVelocity()/299.792<<"c "<<" m: "<<  (*t)->GetParticleDefinition()->GetPDGMass()/GeV<<" GeV"
                      <<std::endl;
        }
}
*/

if(abs(TrPDGid)==22){
    CreateTree::Instance()->ion_z->Fill(global_z, energyIon);
}
if(abs(TrPDGid)==22 || abs(TrPDGid)==11){
    CreateTree::Instance()->h_time_z_egamma->Fill(global_z,thePostPoint->GetGlobalTime() / ns - (global_z+1500) / 300.,energy / GeV);
}else{
    CreateTree::Instance()->h_time_z_other->Fill(global_z,thePostPoint->GetGlobalTime() / ns - (global_z+1500) / 300.,energy / GeV);
}
    double inter_len = 223; //mm
    CreateTree::Instance()->ion_rz->Fill(sqrt(global_x*global_x+global_y*global_y), global_z, energyIon);
    if(global_z>-1500 && global_z<-1500 + inter_len) CreateTree::Instance()->ion_r1->Fill(sqrt(global_x*global_x+global_y*global_y),energyIon);
    if(global_z>-1500 + inter_len && global_z<-1500 + inter_len*2) CreateTree::Instance()->ion_r2->Fill(sqrt(global_x*global_x+global_y*global_y),energyIon);
    if(global_z>-1500 + inter_len*2 && global_z<-1500 + inter_len*3) CreateTree::Instance()->ion_r3->Fill(sqrt(global_x*global_x+global_y*global_y),energyIon);
    if(global_z>-1500 + inter_len*3 && global_z<-1500 + inter_len*4) CreateTree::Instance()->ion_r4->Fill(sqrt(global_x*global_x+global_y*global_y),energyIon);
    if(global_z>-1500 + inter_len*4 && global_z<-1500 + inter_len*5) CreateTree::Instance()->ion_r5->Fill(sqrt(global_x*global_x+global_y*global_y),energyIon);
    if(global_z>-1500 + inter_len*5 && global_z<-1500 + inter_len*6) CreateTree::Instance()->ion_r6->Fill(sqrt(global_x*global_x+global_y*global_y),energyIon);
    if(global_z>-1500 + inter_len*6 && global_z<-1500 + inter_len*7) CreateTree::Instance()->ion_r7->Fill(sqrt(global_x*global_x+global_y*global_y),energyIon);

/*
  if(theTrack->GetDefinition()->GetParticleName().contains("pi") || theTrack->GetDefinition()->GetParticleName().contains("eta") || theTrack->GetDefinition()->GetParticleName().contains("kaon")){
    if(theStep->GetPostStepPoint()->GetKineticEnergy()==0)
      CreateTree::Instance()->depositedEnergyTotal += theTrack->GetDefinition()->GetPDGMass()/ GeV;
  }
*/
/*
  if(theTrack->GetParentID() == 0 && thePostPoint->GetProcessDefinedStep()->GetProcessName().contains("pi-Inelastic")){
      cout<<"pi-Inelastic"<<endl;
      CreateTree::Instance()->depositedEnergyTotal += 100;
  }
*/
/*
  if(theTrack->GetDefinition()->GetParticleName().contains("anti")){
    if(theStep->GetPostStepPoint()->GetKineticEnergy()==0)
      CreateTree::Instance()->depositedEnergyTotal += 2*theTrack->GetDefinition()->GetPDGMass()/ GeV;
  }
*/
  bool outworld = ((theStep->GetPostStepPoint())->GetStepStatus()) == fWorldBoundary;
  if (outworld)
  {
/*
            if(theTrack->GetDefinition()->GetParticleName()=="pi0") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName()=="pi-") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName()=="pi+") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName()=="eta_prime") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName()=="kaon0S") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName()=="kaon0L") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName()=="kaon+") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName()=="kaon-") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName()=="eta") CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
            else if(theTrack->GetDefinition()->GetParticleName().contains("anti")) CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV + theTrack->GetDefinition()->GetPDGMass()/ GeV;
            else CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetKineticEnergy()/GeV;
*/
  CreateTree::Instance()->depositedEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetKineticEnergy()/GeV;
  CreateTree::Instance()->depositedTotalEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
  CreateTree::Instance()->pdgid_escape->push_back( TrPDGid);
  CreateTree::Instance()->KineticEnergy_escape->push_back((theStep->GetPostStepPoint())->GetKineticEnergy() / GeV);
  CreateTree::Instance()->positionx_escape->push_back(thePostPosition.x() / mm);
  CreateTree::Instance()->positiony_escape->push_back(thePostPosition.y() / mm);
  CreateTree::Instance()->positionz_escape->push_back(thePostPosition.z() / mm);
  }

  //if(nStep==1 && theTrack->GetDefinition()->GetParticleName().contains("anti")) std::cout<<"place:"<<thePostPVName<<" "<<theTrack->GetDefinition()->GetParticleName()<<std::endl;

//if (theTrack->GetParentID() == 0 &&  ((thePostPoint->GetProcessDefinedStep()->GetProcessName() == "pi-Inelastic") || (thePostPoint->GetProcessDefinedStep()->GetProcessName() == "pi+Inelastic"))) {
//if  (theTrack->GetParentID() == 0 && ((thePostPoint->GetProcessDefinedStep()->GetProcessName() == "pi-Inelastic") || (thePostPoint->GetProcessDefinedStep()->GetProcessName() == "pi+Inelastic"))) {
/*
if(theTrack->GetDefinition()->GetParticleName().contains("gamma")){

//if (thePostPoint->GetProcessDefinedStep()->GetProcessName() == "Cerenkov"){


        std::cout<<"place:"<<thePostPVName<<" "<<std::endl;//theTrack->GetDefinition()->GetParticleName()<<std::endl;
        std::cout << "**process name: " << thePostPoint->GetProcessDefinedStep()->GetProcessName() << std::endl;
        std::cout << "**TrackID: " << theStep->GetTrack()->GetTrackID() << std::endl;
        std::cout << "**total energy pre_step: " << thePrePoint->GetTotalEnergy() / GeV<<"GeV, post_step: "<<thePostPoint->GetTotalEnergy() / GeV <<" GeV"<< std::endl;
        std::cout << "**kinetic energy pre_step: " << thePrePoint->GetKineticEnergy() / GeV<<"GeV, post_step: "<<thePostPoint->GetKineticEnergy() / GeV <<" GeV"<< std::endl;
        std::cout << "PreEnergy - PostEnergy: "<<thePrePoint->GetTotalEnergy() / GeV - thePostPoint->GetTotalEnergy() / GeV<<std::endl;
        std::cout << "GetTotalEnergyDeposit: "<<theStep->GetTotalEnergyDeposit() /GeV<<" ion deposit"<<energyIon /GeV<<std::endl;
        std::cout << "Total Lost:  "<<thePrePoint->GetTotalEnergy() / GeV - thePostPoint->GetTotalEnergy() / GeV - theStep->GetTotalEnergyDeposit() /GeV << std::endl;
        std::cout << " Secondary: ------------------------------ " << std::endl;

        const std::vector<const G4Track*>* secondaries = theStep->GetSecondaryInCurrentStep();
        if(secondaries->size()==0) std::cout<<"No secondaries"<<std::endl;
        double sum_Etrans=0;
        for (auto t = secondaries->begin(); t != secondaries->end(); t++) {
        //if((*t)->GetTrackStatus()==fAlive){
            if((*t)->GetParticleDefinition()->GetParticleName()=="pi0") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName()=="pi-") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName()=="pi+") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName()=="eta_prime") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName()=="kaon0S") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName()=="kaon0L") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName()=="kaon+") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName()=="kaon-") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName()=="eta") sum_Etrans+=(*t)->GetTotalEnergy()/GeV;
            else if((*t)->GetParticleDefinition()->GetParticleName().contains("anti")) sum_Etrans+=(*t)->GetTotalEnergy()/GeV + (*t)->GetParticleDefinition()->GetPDGMass()/ GeV;
            else sum_Etrans+=(*t)->GetKineticEnergy()/GeV;

            std::cout << "**"<<(*t)->GetParticleDefinition()->GetParticleName()
                      << " ParentID: "<<(*t)->GetParentID()
                      << " IsBelowThreshold: "<<(*t)->IsBelowThreshold()
                      <<std::endl;
            std::cout << " Etot: "<<(*t)->GetTotalEnergy()/GeV<<" Ek: "<<(*t)->GetKineticEnergy()/GeV//<<" p: "<<(*t)->GetMomentum()
                      <<" v: "<<(*t)->GetVelocity()/299.792<<"c "<<" m: "<<  (*t)->GetParticleDefinition()->GetPDGMass()/GeV<<" GeV"
                      <<std::endl;

        }
        std::cout << " Finish: ------------------------------ " << std::endl;
            //std::cout<<"imbalance: "<<thePrePoint->GetTotalEnergy() / GeV - thePostPoint->GetTotalEnergy() / GeV - theStep->GetTotalEnergyDeposit() /GeV - sum_Etrans<<std::endl;
            //CreateTree::Instance()->primary_Eimbalance->push_back(thePrePoint->GetTotalEnergy() / GeV - thePostPoint->GetTotalEnergy() / GeV - theStep->GetTotalEnergyDeposit() /GeV - sum_Etrans);

}

*/

  //------------- optical photon -------------
  if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
  {
    //if optics
    G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
    float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy() / eV);
    //only consider 300 to 1000nm
    if (photWL > 1000 || photWL < 300)
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    else
    {
      //only consider Cerenkov and Scintillation
      if ((processName == "Cerenkov") || processName == "Scintillation")
      {

        TrackInformation *theTrackInfo = (TrackInformation *)(theTrack->GetUserInformation());
        G4int aapdgid = theTrackInfo->GetParentPDGid();
        //if (thePostPVName.contains("ecalDetWindowP_ff_"))
        if ((thePostPVName.contains("ecalGapP_ff_")))
        {
          if (processName == "Scintillation")
          {
            CreateTree::Instance()->h_photon_2D_receive_Scin->Fill(theTrack->GetVertexPosition()[1], theTrack->GetVertexPosition()[2]);
            CreateTree::Instance()->h_phot_lambda_ECAL_f_Scin->Fill(photWL);
            CreateTree::Instance()->SDdetected_ff_S++;
            CreateTree::Instance()->h_phot_detect_time_f_Scin->Fill(thePostPoint->GetGlobalTime() / ns);
          }
          if (processName == "Cerenkov")
          {
            CreateTree::Instance()->h_photon_2D_receive_Ceren->Fill(theTrack->GetVertexPosition()[1], theTrack->GetVertexPosition()[2]);
            CreateTree::Instance()->h_phot_lambda_ECAL_f_Ceren->Fill(photWL);
            CreateTree::Instance()->SDdetected_ff_C++;
            CreateTree::Instance()->h_phot_detect_time_f_Ceren->Fill(thePostPoint->GetGlobalTime() / ns);
          }
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
        //else if (thePostPVName.contains("ecalDetWindowP_fr_"))
        else if (thePostPVName.contains("ecalGapP_fr_"))
        {
          if (processName == "Scintillation")
          {
            CreateTree::Instance()->h_photon_2D_receive_Scin->Fill(theTrack->GetVertexPosition()[1], theTrack->GetVertexPosition()[2]);
            CreateTree::Instance()->h_phot_lambda_ECAL_r_Scin->Fill(photWL);
            CreateTree::Instance()->SDdetected_rr_S++;
            CreateTree::Instance()->h_phot_detect_time_r_Scin->Fill(thePostPoint->GetGlobalTime() / ns);
          }
          if (processName == "Cerenkov")
          {
            CreateTree::Instance()->h_photon_2D_receive_Ceren->Fill(theTrack->GetVertexPosition()[1], theTrack->GetVertexPosition()[2]);
            CreateTree::Instance()->h_phot_lambda_ECAL_r_Ceren->Fill(photWL);
            CreateTree::Instance()->SDdetected_rr_C++;
            CreateTree::Instance()->h_phot_detect_time_r_Ceren->Fill(thePostPoint->GetGlobalTime() / ns);
          }
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }


        if (thePrePVName.contains("world"))
        {
          //theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
        else if (thePrePVName.contains("ecalGapP"))
        {
          //std::cout << "in Gap " << thePrePVName << std::endl;
        }
        else if (thePrePVName.contains("wrap"))
        {
          //std::cout << "Crystal " << thePrePVName << std::endl;
        }
        else
        {
          //std::cout << "weird PrePVName " << thePrePVName << std::endl;
        }
        if (thePostPVName.contains("world"))
        {
          //theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
        //std::cout << "out of world " << thePrePVName << std::endl;
        G4double tracklength = theStep->GetTrack()->GetTrackLength();
        if (tracklength > maxtracklength)
        {
          theStep->GetTrack()->SetTrackStatus(fStopAndKill);
          std::cout << "maximum " << thePrePVName << std::endl;
        }
        if (!propagateCerenkov && (processName == "Cerenkov"))
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);

        if (!propagateScintillation && (processName == "Scintillation"))
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);

        if ((nStep == 1) && thePrePVName.contains("ecal"))
        {
          //theTrack->GetDefinition()->GetParticleName()
          //cout<<"position x:"<<global_x<<"  z: "<<global_z<<endl;
          if (processName == "Scintillation")
            CreateTree::Instance()->h_photon_2D_produce_Scin->Fill(global_y, global_z);
          if (processName == "Cerenkov")
            CreateTree::Instance()->h_photon_2D_produce_Ceren->Fill(global_y, global_z);
          if (thePrePVName.contains("ecalCrystalP_f"))
          {
            if (processName == "Scintillation")
            {
              CreateTree::Instance()->h_phot_lambda_ECAL_f_produce_Scin->Fill(photWL);
              CreateTree::Instance()->ECAL_f_total_S += 1;
            }
            if (processName == "Cerenkov")
            {
              CreateTree::Instance()->h_phot_lambda_ECAL_f_produce_Ceren->Fill(photWL);
              CreateTree::Instance()->ECAL_f_total_C += 1;
            }
          }
          if (thePrePVName.contains("ecalCrystalP_r"))
          {
            CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_total += 1;
            if (processName == "Scintillation")
            {
              CreateTree::Instance()->h_phot_lambda_ECAL_r_produce_Scin->Fill(photWL);
              CreateTree::Instance()->ECAL_r_total_S += 1;
            }
            if (processName == "Cerenkov")
            {
              CreateTree::Instance()->h_phot_lambda_ECAL_r_produce_Ceren->Fill(photWL);
              CreateTree::Instance()->ECAL_r_total_C += 1;
            }
          }
          TrackInformation *theTrackInfo = (TrackInformation *)(theTrack->GetUserInformation());
          G4int aapdgid = theTrackInfo->GetParentPDGid();

        }

      }
      else
      {

        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      }
    }
  }
  else
  {

    //count tracks before SCEPCAL at the tracker layers
    if (thePrePVName.contains("world") && thePostPVName.contains("trackerPV_Layer")) // interface between T1 and T2
    {
      for (int iLayer = 0; iLayer < 6; iLayer++)
      {
        if (thePostPVName == Form("trackerPV_Layer%d", iLayer)) // interface between T1 and T2
          CreateTree::Instance()->nTracksTRK[iLayer]++;         //counting tracks crossing the boundary
      }
    }

    //save primary particle position and momentum before timing layer T1 and before ECAL E1
    else if (thePrePVName.contains("world") && thePostPVName.contains("corePV_front")) // interface between world and T1
    {

      CreateTree::Instance()->nTracksT1++; //counting tracks crossing the boundary

      if (theTrack->GetParentID() == 0) // select only the primary particle
      {
        CreateTree::Instance()->primaryPosT1->at(0) = global_x;
        CreateTree::Instance()->primaryPosT1->at(1) = global_y;
        CreateTree::Instance()->primaryPosT1->at(2) = global_z;

        CreateTree::Instance()->primaryMomT1->at(0) = thePrePoint->GetMomentum().x() / GeV;
        CreateTree::Instance()->primaryMomT1->at(1) = thePrePoint->GetMomentum().y() / GeV;
        CreateTree::Instance()->primaryMomT1->at(2) = thePrePoint->GetMomentum().z() / GeV;
        CreateTree::Instance()->primaryMomT1->at(3) = thePrePoint->GetTotalEnergy() / GeV;
      }
    }

    else if (thePrePVName.contains("world") && thePostPVName.contains("corePV_rear")) // interface between T1 and T2
    {
      CreateTree::Instance()->nTracksT2++; //counting tracks crossing the boundary
    }

    else if ((thePrePVName.contains("world") || thePrePVName.contains("ecalGapP_f") || thePrePVName.contains("ecalDetP_f")) && thePostPVName.contains("ecalCrystalP_f") // interface between world and E1
    )
    {
      CreateTree::Instance()->nTracksE1++; //counting tracks crossing the boundary

      if (theTrack->GetParentID() == 0) // select only the primary particle
      {
        CreateTree::Instance()->primaryPosE1->at(0) = global_x;
        CreateTree::Instance()->primaryPosE1->at(1) = global_y;
        CreateTree::Instance()->primaryPosE1->at(2) = global_z;

        CreateTree::Instance()->primaryMomE1->at(0) = thePrePoint->GetMomentum().x() / GeV;
        CreateTree::Instance()->primaryMomE1->at(1) = thePrePoint->GetMomentum().y() / GeV;
        CreateTree::Instance()->primaryMomE1->at(2) = thePrePoint->GetMomentum().z() / GeV;
        CreateTree::Instance()->primaryMomE1->at(3) = thePrePoint->GetTotalEnergy() / GeV;
      }
    }

    else if ((thePrePVName.contains("ecalCrystalP_f") || thePrePVName.contains("world")) && thePostPVName.contains("ecalCrystalP_r")) // interface between E1 and E2
    {
      CreateTree::Instance()->nTracksE2++; //counting tracks crossing the boundary
    }

    //tracker
    if (thePrePVName.contains("trackerPV_Layer"))
    {
      for (int iLayer = 0; iLayer < 6; iLayer++)
      {
        if (thePrePVName == Form("trackerPV_Layer%d", iLayer))
          CreateTree::Instance()->Edep_Tracker_layer[iLayer] += energy / GeV;
      }
    }

    //timing
    if (thePrePVName.contains("corePV_front"))
    {
      CreateTree::Instance()->depositedEnergyTiming_f += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTiming_f += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTiming_f += energyElec / GeV;
      for (int iBar = 0; iBar < 18; iBar++)
      {
        if (thePrePVName == Form("corePV_front_%d", iBar))
          CreateTree::Instance()->Edep_Timing_f_ch[iBar] += energy / GeV;
      }
    }
    if (thePrePVName.contains("corePV_rear"))
    {
      CreateTree::Instance()->depositedEnergyTiming_r += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTiming_r += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTiming_r += energyElec / GeV;
      for (int iBar = 0; iBar < 18; iBar++)
      {
        if (thePrePVName == Form("corePV_rear_%d", iBar))
          CreateTree::Instance()->Edep_Timing_r_ch[iBar] += energy / GeV;
      }
    }

    //ecal
    int chn = 64;
    if (thePrePVName.contains("ecalCrystalP_f"))
    {
      for (int iCh = 0; iCh < chn; iCh++)
      {
        if (thePrePVName.contains(Form("_%d", iCh)))
        {
          CreateTree::Instance()->Edep_ECAL_f_ch[iCh] += energy / GeV;
          CreateTree::Instance()->IonEdep_ECAL_f_ch[iCh] += energyIon / GeV;
          //          CreateTree::Instance()->ElecEdep_ECAL_f_ch[iCh] += energyElec / GeV;
        }
      }

      if (thePrePVName.contains("ecalCrystalP_f"))
      {

        CreateTree::Instance()->depositedEnergyECAL_f[0] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_f[0] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_f[0] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[5] += energyKaon / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[7] += energyIonProton / GeV;
      }
    }
/*
      if (thePrePVName.contains("ecalCrystalP_f_fiber_scinti"))
      {
        CreateTree::Instance()->depositedEnergyECAL_f[1] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_f[1] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_f[1] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[5] += energyKaon / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_f_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_f_particleID[7] += energyIonProton / GeV;
       
        for (int iCh = 0; iCh < chn; iCh++)
        {
          if (thePrePVName == Form("ecalCrystalP_f_fiber_scinti_%d", iCh))
          {
            CreateTree::Instance()->Edep_ECAL_f_scinti_ch[iCh] += energy / GeV;
            CreateTree::Instance()->IonEdep_ECAL_f_scinti_ch[iCh] += energyIon / GeV;
            CreateTree::Instance()->ElecEdep_ECAL_f_scinti_ch[iCh] += energyElec / GeV;
          }
        }
      }
      if (thePrePVName.contains("ecalCrystalP_f_fiber_cheren"))
      {
        CreateTree::Instance()->depositedEnergyECAL_f[2] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_f[2] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_f[2] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[5] += energyKaon / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_f_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_f_particleID[7] += energyIonProton / GeV;
      }
    }

    if (thePrePVName.contains("ecalCrystalP_r"))
    {
      for (int iCh = 0; iCh < chn; iCh++)
      {
        if (thePrePVName.contains(Form("_%d", iCh)))
        {
          CreateTree::Instance()->Edep_ECAL_r_ch[iCh] += energy / GeV;
          CreateTree::Instance()->IonEdep_ECAL_r_ch[iCh] += energyIon / GeV;
          //          CreateTree::Instance()->ElecEdep_ECAL_r_ch[iCh] += energyElec / GeV;
        }
      }
      if (thePrePVName.contains("ecalCrystalP_r_absorb"))
      {
        CreateTree::Instance()->depositedEnergyECAL_r[0] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_r[0] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_r[0] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->betaparticleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->betaparticleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->betaparticleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->betaparticleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->betaparticleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[5] += energyKaon / GeV;
        CreateTree::Instance()->betaparticleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->betaparticleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_absorb_r_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->betaparticleID[7] += energyIonProton / GeV;
      }
      if (thePrePVName.contains("ecalCrystalP_r_fiber_scinti"))
      {
        CreateTree::Instance()->depositedEnergyECAL_r[1] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_r[1] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_r[1] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[5] += energyKaon / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_scinti_r_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_scinti_r_particleID[7] += energyIonProton / GeV;
      }

      if (thePrePVName.contains("ecalCrystalP_r_fiber_cheren"))
      {
        CreateTree::Instance()->depositedEnergyECAL_r[2] += energy / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_r[2] += energyIon / GeV;
        CreateTree::Instance()->depositedElecEnergyECAL_r[2] += energyElec / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[0] += energyPion_n / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[0] += energyIonPion_n / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[1] += energyPositron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[1] += energyIonPositron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[2] += energyElectron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[2] += energyIonElectron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[3] += energyPhoton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[3] += energyIonPhoton / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[4] += energyPion_p / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[4] += energyIonPion_p / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[5] += energyKaon / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[5] += energyIonKion / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[6] += energyNeutron / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[6] += energyIonNeutron / GeV;

        CreateTree::Instance()->depositedEnergyECAL_cheren_r_particleID[7] += energyProton / GeV;
        CreateTree::Instance()->depositedIonEnergyECAL_cheren_r_particleID[7] += energyIonProton / GeV;
      }
    }
    //    if (thePrePVName.contains("ecalCrystalL_f_absorb") || thePrePVName.contains("ecalCrystalL_r_absorb"))
    if (thePrePVName.contains("ecalCrystalP_r_fiber_scinti") || thePrePVName.contains("ecalCrystalP_f_fiber_scinti"))
    {

      int iiZ, iiT;
      if (global_z < 0)
        iiZ = 0;
      if (global_z > 5000)
        iiZ = 2499;
      if (global_z >= 0 && global_z < 5000)
        iiZ = (global_z / 2);

      if (((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3) < 0)
        iiT = 0;
      if (((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3) > 5)
        iiT = 2499;
      if (((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3) >= 0 && ((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3) < 5)
        iiT = (((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3) * 500);

      //cout<<"iiZ = "<<iiZ<<" iiT = "<<iiT<<endl;

      CreateTree::Instance()->E_Zdep_0to5000mm_total[iiZ] += energyIon / GeV;
      CreateTree::Instance()->E_Zdep_0to5000mm_Pion_n[iiZ] += energyIonPion_n / GeV;
      CreateTree::Instance()->E_Zdep_0to5000mm_Positron[iiZ] += energyIonPositron / GeV;
      CreateTree::Instance()->E_Zdep_0to5000mm_Electron[iiZ] += energyIonElectron / GeV;
      CreateTree::Instance()->E_Zdep_0to5000mm_Photon[iiZ] += energyIonPhoton / GeV;
      CreateTree::Instance()->E_Zdep_0to5000mm_Pion_p[iiZ] += energyIonPion_p / GeV;
      CreateTree::Instance()->E_Zdep_0to5000mm_Neutron[iiZ] += energyIonNeutron / GeV;
      CreateTree::Instance()->E_Zdep_0to5000mm_Kion[iiZ] += energyIonKion / GeV;
      CreateTree::Instance()->E_Zdep_0to5000mm_Proton[iiZ] += energyIonProton / GeV;

      CreateTree::Instance()->E_Tdep_0to5ns_total[iiT] += energyIon / GeV;
      CreateTree::Instance()->E_Tdep_0to5ns_Pion_n[iiT] += energyIonPion_n / GeV;
      CreateTree::Instance()->E_Tdep_0to5ns_Positron[iiT] += energyIonPositron / GeV;
      CreateTree::Instance()->E_Tdep_0to5ns_Electron[iiT] += energyIonElectron / GeV;
      CreateTree::Instance()->E_Tdep_0to5ns_Photon[iiT] += energyIonPhoton / GeV;
      CreateTree::Instance()->E_Tdep_0to5ns_Pion_p[iiT] += energyIonPion_p / GeV;
      CreateTree::Instance()->E_Tdep_0to5ns_Kion[iiT] += energyIonKion / GeV;
      CreateTree::Instance()->E_Tdep_0to5ns_Neutron[iiT] += energyIonNeutron / GeV;
      CreateTree::Instance()->E_Tdep_0to5ns_Proton[iiT] += energyIonProton / GeV;

    }
*/
    //hcal
    if (thePrePVName.contains("hcalTile_layer"))
    {
      CreateTree::Instance()->depositedEnergyHCALAct += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyHCALAct += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyHCALAct += energyElec / GeV;
    }
    if (thePrePVName.contains("hcalAbs"))
    {
      CreateTree::Instance()->depositedEnergyHCALPas += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyHCALPas += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyHCALPas += energyElec / GeV;
    }

    if (thePrePVName.contains("world"))
    {
      CreateTree::Instance()->depositedEnergyWorld += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyWorld += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyWorld += energyElec / GeV;
    }

    if (thePrePVName.contains("services"))
    {
      CreateTree::Instance()->depositedEnergyServices += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyServices += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyServices += energyElec / GeV;
    }

    if (thePrePVName.contains("TimingGap"))
    {
      CreateTree::Instance()->depositedEnergyTimingGap += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTimingGap += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTimingGap += energyElec / GeV;
    }

    if (thePrePVName.contains("ecalGap"))
    {
      CreateTree::Instance()->depositedEnergyEcalGap += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyEcalGap += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyEcalGap += energyElec / GeV;
    }

    if (thePrePVName.contains("ecalDet"))
    {
      CreateTree::Instance()->depositedEnergyEcalDet += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyEcalDet += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyEcalDet += energyElec / GeV;
    }

    if (thePrePVName.contains("solenoid"))
    {
      CreateTree::Instance()->depositedEnergySolenoid += energy / GeV;
      CreateTree::Instance()->depositedIonEnergySolenoid += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergySolenoid += energyElec / GeV;
    }

    //G4cout << ">>> end non optical photon" << G4endl;
  } // non optical photon

  return;
}
