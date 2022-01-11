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
  G4double energyOther = 0.;
  //ion energy by particle types
  G4double energyIonPion_n = 0.;
  G4double energyIonPositron = 0.;
  G4double energyIonElectron = 0.;
  G4double energyIonPhoton = 0.;
  G4double energyIonPion_p = 0.;
  G4double energyIonKaon = 0.;
  G4double energyIonNeutron = 0.;
  G4double energyIonProton = 0.;
  G4double energyIonOther =0.;

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
    energyIonKaon = energyIon;
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
  else {
    energyOther = energy;
    energyIonOther =energyIon;
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
  float aabc2 = thePostPoint->GetGlobalTime() / ns - (global_z+10)/300;
  float aabc3=aabc2;
  if(aabc2>4.999) aabc2=4.999;
  if(aabc3>499) aabc3=499;
  CreateTree::Instance()->pdg_time->Fill(TrPDGid,aabc2,energyIon / GeV);
  CreateTree::Instance()->pdg_time2->Fill(TrPDGid,aabc3,energyIon / GeV);

 
  int itime=aabc2/0.25;
  if(itime<0) itime=10;
  if(itime>79) itime=10;

  if(TrPDGid==-211) {
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[0][itime] += energyIon / GeV;
  } 
  else if(TrPDGid==-11) {
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[1][itime] += energyIon / GeV;
  }
  else if(TrPDGid==11) {
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[2][itime] += energyIon / GeV;
  }
  else if(TrPDGid==22) {
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[3][itime] += energyIon / GeV;
  }
  else if(TrPDGid==211) {
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[4][itime] += energyIon / GeV;
  }
  else if(TrPDGid==321) {
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[5][itime] += energyIon / GeV;
  }
  else if(TrPDGid==2112) {
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[6][itime] += energyIon / GeV;
  }
  else if(TrPDGid==2212) {
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[7][itime] += energyIon / GeV;
  }
  else {
    if(fabs(TrPDGid)<1000000000) std::cout<<" other pid is "<<TrPDGid<<std::endl;
    CreateTree::Instance()->depositedIonEnergyECAL_pidtime[8][itime] += energyIon / GeV;
  }



  if(nStep==1) {
    CreateTree::Instance()->nNeutrons++;
    CreateTree::Instance()->h_keneutrons->Fill(theStep->GetPreStepPoint()->GetKineticEnergy()/GeV);
  }

//if(TrPDGid>=111) cout<<energyIon / GeV<<endl;

  if(thePostPoint->GetProcessDefinedStep()->GetProcessName().contains("Inelast")){
    CreateTree::Instance()->Ninelastic++;
    CreateTree::Instance()->inelasticEK->Fill((theStep->GetPreStepPoint()->GetKineticEnergy() - theStep->GetPostStepPoint()->GetKineticEnergy())/GeV);

  }


  if(abs(TrPDGid)==22){
    CreateTree::Instance()->ion_z->Fill(global_z, energyIon);
  }


  float aabc = thePostPoint->GetGlobalTime() / ns - (global_z+10)/300;
  if(abs(TrPDGid)==22 || abs(TrPDGid)==11){
    CreateTree::Instance()->h_time_z_egamma->Fill(global_z, aabc,energy / GeV);

    
  } else if(abs(TrPDGid)==2212) {
    CreateTree::Instance()->h_time_z_proton->Fill(global_z,aabc ,energy / GeV);
  } else if(abs(TrPDGid)==2112) {
    CreateTree::Instance()->h_time_z_neutron->Fill(global_z,aabc ,energy / GeV);
  } else{
    CreateTree::Instance()->h_time_z_othernotpn->Fill(global_z,aabc ,energy / GeV);
    //    if(aabc>50) std::cout<<"late time "<<aabc<<" for pid "<<TrPDGid<<" with energy deposit "<<energy<<std::endl;
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

  bool outworld = ((theStep->GetPostStepPoint())->GetStepStatus()) == fWorldBoundary;
  if (outworld)
    {
      CreateTree::Instance()->kineticEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetKineticEnergy()/GeV;
      CreateTree::Instance()->TotalEnergyEscapeWorld+=(theStep->GetPostStepPoint())->GetTotalEnergy()/GeV;
      CreateTree::Instance()->pdgid_escape->push_back( TrPDGid);
      CreateTree::Instance()->KineticEnergy_escape->push_back((theStep->GetPostStepPoint())->GetKineticEnergy() / GeV);
      CreateTree::Instance()->positionx_escape->push_back(thePostPosition.x() / mm);
      CreateTree::Instance()->positiony_escape->push_back(thePostPosition.y() / mm);
      CreateTree::Instance()->positionz_escape->push_back(thePostPosition.z() / mm);
    }


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

        if (thePrePVName.contains("world"))
        {
          //theTrack->SetTrackStatus(fKillTrackAndSecondaries);
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
  else  // not optical photons
  {
    if(thePrePVName!="ecalCrystalP_f_0"&&thePrePVName!="worldPV"&&
       thePrePVName!="ecalWrapperP_f_0")
    std::cout<<"prepvname is "<<thePrePVName<<std::endl;
    //count tracks before SCEPCAL at the tracker layers

    if ((thePrePVName.contains("world") || thePrePVName.contains("ecalGapP_f") || thePrePVName.contains("ecalDetP_f")) && thePostPVName.contains("ecalCrystalP_f") // interface between world and E1
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

      CreateTree::Instance()->depositedEnergyECAL_f += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyECAL_f += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyECAL_f += energyElec / GeV;

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
      CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[5] += energyIonKaon / GeV;

      CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[6] += energyNeutron / GeV;
      CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[6] += energyIonNeutron / GeV;

      CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[7] += energyProton / GeV;
      CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[7] += energyIonProton / GeV;


      CreateTree::Instance()->depositedEnergyECAL_absorb_f_particleID[8] += energyOther / GeV;
      CreateTree::Instance()->depositedIonEnergyECAL_absorb_f_particleID[8] += energyIonOther / GeV;
    }
    


    if (thePrePVName.contains("world"))
    {
      CreateTree::Instance()->depositedEnergyWorld += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyWorld += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyWorld += energyElec / GeV;
    }



    if (thePrePVName.contains("Wrap"))
    {
      CreateTree::Instance()->depositedEnergyWrap += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyWrap += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyWrap += energyElec / GeV;
    }




    //G4cout << ">>> end non optical photon" << G4endl;
  } // non optical photon

  return;
}
