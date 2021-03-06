// ====================================================================
//   CH1SD.cc
//
// ====================================================================
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "CH1SD.hh"
#include "CH1Hit.hh"

////////////////////////////////////////////////
  CH1SD::CH1SD( G4String name)
: G4VSensitiveDetector(name)
  ////////////////////////////////////////////////
{
  collectionName.insert("Ch1Collection");
}

/////////////////////////////////
CH1SD::~CH1SD()
  /////////////////////////////////
{
}


////////////////////////////////////////////////
void CH1SD::Initialize(G4HCofThisEvent* HCTE)
  ////////////////////////////////////////////////
{

  static int HCID = -1;
  Ch1Collection = new CH1HitsCollection( SensitiveDetectorName,
                                         collectionName[0] );
  if( HCID<0 )
    {
      HCID = GetCollectionID(0);
    }
  HCTE->AddHitsCollection( HCID, Ch1Collection );

}


///////////////////////////////////////////////////////////
G4bool CH1SD::ProcessHits(G4Step* aStep,
    G4TouchableHistory* ROhist)
  ///////////////////////////////////////////////////////////
{

  // get step information from "PreStepPoint"
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  G4String particleName;

  //if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
  //  return false;

  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();
  /*
  // e+/e- rejection
  if( particleName == "e-")
  return false;
  if( particleName == "e+")
  return false;
  */
  const G4Track* aTrack = aStep->GetTrack();
  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

  //if(particleType == "lepton")
  //  return false;
  //
  //  if( (particleName != "kaon+")){}
  //    return false;
  //  {}

  //  if( (particleName != "pi-") && (particleName != "pi+")){}
  //    return false;
  //  {}
  //  if((particleName != "pi+")&& (particleName != "pi-")
  //     && (particleName != "proton")){}
  //    return false;
  //  {}
  if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();


  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack -> GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m

  G4ThreeVector pos= preStepPoint-> GetPosition();
  G4ThreeVector mom= preStepPoint-> GetMomentum();
  G4double tof= preStepPoint-> GetGlobalTime();
  G4int tid =  aStep-> GetTrack()-> GetTrackID();
  G4int pid =  aStep-> GetTrack()-> GetDefinition() -> GetPDGEncoding();
  G4double mass =  aStep-> GetTrack()-> GetDynamicParticle() -> GetMass();
  G4int qq =  aStep-> GetTrack()-> GetDynamicParticle() -> GetCharge();
  G4double tlength = aStep->GetTrack()-> GetTrackLength();
  //  G4double slength = aStep->GetTrack()-> GetStepLength();
  //  G4cout<<mass<<":"<<tlength<<":"<<slength<<G4endl;
  G4int parentID =  aStep-> GetTrack()-> GetParentID();
  //    G4cout <<"test: "<<parentID << G4endl;
  // Get Pad number
  G4int iDet;
  G4String name = physVol->GetName();
  //  G4cout << name << G4endl;

  //  sscanf(name,"CH1PV%d",&iDet);
  G4int copyNo = preStepPoint -> GetPhysicalVolume()->GetCopyNo();
  iDet=copyNo;

  // create a new hit and push them to "Hit Coleltion"
  CH1Hit* ahit= new CH1Hit(pos, mom, tof, tid, pid, iDet,mass,qq,parentID, VertexPosition, VertexMomentum, VertexEnergy,tlength);
  Ch1Collection-> insert(ahit);

  return true;
}

////////////////////////////////////////////////////
void CH1SD::EndOfEvent(G4HCofThisEvent* HCTE)
  ////////////////////////////////////////////////////
{
}

////////////////////////////
void CH1SD::DrawAll()
  ////////////////////////////
{
}

/////////////////////////////
void CH1SD::PrintAll()
  /////////////////////////////
{
  Ch1Collection-> PrintAllHits();
}
