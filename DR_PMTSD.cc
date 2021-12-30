
#include "DR_PMTHit.hh"
#include "DR_PMTSD.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "CreateTree.hh"

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"


DR_PMTSD::DR_PMTSD(G4String name)
    : G4VSensitiveDetector(name),
      fPMTHitsCollection(nullptr), HCID(-1)
{
  collectionName.insert("PMTHitsCollection");
}

DR_PMTSD::~DR_PMTSD() {}

void DR_PMTSD::Initialize(G4HCofThisEvent *HCE)
{
  fPMTHitsCollection = new DR_PMTHitsCollection(GetName(), collectionName[0]);
  if (HCID < 0)
  {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(fPMTHitsCollection);
  }
  HCE->AddHitsCollection(HCID, fPMTHitsCollection);

}

G4bool
DR_PMTSD::ProcessHits( G4Step*       aStep,
                                   G4TouchableHistory*  )
{
    return false;
}
G4bool
DR_PMTSD::ProcessHits_constStep( const G4Step* aStep, G4TouchableHistory* )
{

  G4Track *theTrack = aStep->GetTrack();
  float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy() / CLHEP::eV);

  if (theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
  {
    return false;
  }
  G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
  if (true /*processName == "Cerenkov"*/)
  {
    //std::cout<<"SD detected"<<std::endl;
    G4StepPoint *thePrePointsd = aStep->GetPreStepPoint();
    G4VPhysicalVolume *thePrePVsd = thePrePointsd->GetPhysicalVolume();
    G4String thePrePVNamesd = thePrePVsd->GetName();

    G4StepPoint *thePostPointsd = aStep->GetPostStepPoint();
    G4VPhysicalVolume *thePostPVsd = thePostPointsd->GetPhysicalVolume();
    G4String thePostPVNamesd = thePostPVsd->GetName();

    if (thePostPVNamesd.contains("ecalDetP_"))
    {
      DR_PMTHit* hit = new DR_PMTHit(0);
      hit->AddEnergy(theTrack->GetTotalEnergy());
      hit->IncPhotonCount(); // increment hit for the selected pmt
      hit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
      hit->SetLength(theTrack->GetTrackLength());
      hit->SetBounceCount(theTrack->GetCurrentStepNumber());
      fPMTHitsCollection->insert(hit);
   }
    else if (thePostPVNamesd.contains("ecalDetP_fr"))
    {
      DR_PMTHit* hit = new DR_PMTHit(1);
      hit->AddEnergy(theTrack->GetTotalEnergy());
      hit->IncPhotonCount(); // increment hit for the selected pmt
      hit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
      hit->SetLength(theTrack->GetTrackLength());
      hit->SetBounceCount(theTrack->GetCurrentStepNumber());
      fPMTHitsCollection->insert(hit);
    }
    else if (thePostPVNamesd.contains("ecalDetP_rf"))
    {
      DR_PMTHit* hit = new DR_PMTHit(2);
      hit->AddEnergy(theTrack->GetTotalEnergy());
      hit->IncPhotonCount(); // increment hit for the selected pmt
      hit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
      hit->SetLength(theTrack->GetTrackLength());
      hit->SetBounceCount(theTrack->GetCurrentStepNumber());
      fPMTHitsCollection->insert(hit);
    }
    else if (thePostPVNamesd.contains("ecalDetP_rr"))
    {
      DR_PMTHit* hit = new DR_PMTHit(3);
      hit->AddEnergy(theTrack->GetTotalEnergy());
      hit->IncPhotonCount(); // increment hit for the selected pmt
      hit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
      hit->SetLength(theTrack->GetTrackLength());
      hit->SetBounceCount(theTrack->GetCurrentStepNumber());
      fPMTHitsCollection->insert(hit);
     }else{
      const G4ThreeVector &thePrePositionsd = thePrePointsd->GetPosition();
      const G4ThreeVector &thePostPositionsd = thePostPointsd->GetPosition();
      G4Navigator* navigator=G4TransportationManager::GetTransportationManager()->GetNavigator( "worldPV" );
      G4VPhysicalVolume* volumepre = navigator->LocateGlobalPointAndSetup( thePrePositionsd);
      G4VPhysicalVolume* volumepost = navigator->LocateGlobalPointAndSetup( thePostPositionsd );
      std::cout<<"wierd SD place, Pre: "<<volumepre->GetName()<<" ("<<thePrePositionsd.x()<<", "<<thePrePositionsd.y()<<", "<<thePrePositionsd.z()<<")"<<" post: "<<volumepost->GetName()<<" ("<<thePostPositionsd.x()<<", "<<thePostPositionsd.y()<<", "<<thePostPositionsd.z()<<")"<<std::endl;
      return false;
    }
  return true;
  }else{
    return false;
  }
}



void DR_PMTSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int nHits = fPMTHitsCollection->entries();
  if (verboseLevel >= 1)
  {
    G4cout << "[DR_] PMT collection: " << nHits << " hits" << G4endl;
    if (verboseLevel >= 2)
    {
      fPMTHitsCollection->PrintAllHits();
    }
  }
  fPMTHitsCollection->PrintAllHits();
}

void DR_PMTSD::clear() {}

void DR_PMTSD::DrawAll() {}

void DR_PMTSD::PrintAll() {}
