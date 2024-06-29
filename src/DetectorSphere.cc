#include "DetectorSphere.hh"
#include <G4Step.hh>
#include <fstream>
#include <iostream>
#include <G4EventManager.hh>
#include <G4SDManager.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4SystemOfUnits.hh>
#include <G4RunManager.hh>
#include <iomanip>
#include <sstream>

using namespace std;

DetectorSphere::DetectorSphere(const G4String& name) : G4VSensitiveDetector(name)
{

}

DetectorSphere::~DetectorSphere()
{

}

G4bool DetectorSphere::ProcessHits(G4Step* step, G4TouchableHistory*) 
{
 
    // G4Track* track = step->GetTrack();
    // const G4ParticleDefinition* particle = track->GetDefinition();
    // G4String particleName =  step->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
    // G4ThreeVector momentumDirection = track->GetMomentumDirection();
    // if (particleName == "gamma" && track->GetKineticEnergy() >425* keV)
    // {
    //     // if (momentumDirection == G4ThreeVector(0, -1, 0) && momentumDirection == G4ThreeVector(0, 1, 0))
    //     // {

    //     // }
        
    //     track->SetTrackStatus(fStopAndKill);
    // }
    return true;
}



