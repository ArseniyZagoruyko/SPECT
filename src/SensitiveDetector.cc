#include "SensitiveDetector.hh"
#include <G4Step.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

SensitiveDetector::SensitiveDetector(G4String name, G4int id) : G4VSensitiveDetector(name), detectorID(id) 
{

}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) 
{


    return true;
}

G4int SensitiveDetector::GetDetectorID() const
{
    return detectorID;
}

SensitiveDetector::~SensitiveDetector() 
{

}