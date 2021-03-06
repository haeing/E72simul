/*
  TargetHit.cc from SKS
  2017/10  Yang
*/

#include "TargetHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>

TargetHit::TargetHit()
  : layerID_(0), segID_(0), time_(-100.0), edep_(0.0),
    trackNo_(0), fSignal_(false), hitdiff_(0.)
{}

G4Allocator<TargetHit> TargetHitAllocator;

void TargetHit::Draw()
{
  G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();
  if( pVisManager )
    {
      G4Circle circle( pos_ );
      circle.SetScreenSize( 3.0 );              // in pixels
      circle.SetFillStyle( G4Circle::filled );
      G4Colour colour( 0., 1., 0. );            // green
      G4VisAttributes attribs( colour );
      circle.SetVisAttributes( attribs );
      pVisManager->Draw( circle );
    }
}

void TargetHit::Print()
{
  G4cout << "TargetHit " << std::setw(2) << layerID_ << "-" 
	 << std::setw(2) << segID_ << " > Time: " 
	 << time_/ns << " [ns] Edep.: " << edep_/MeV
	 << " [MeV]" << G4endl;
}

