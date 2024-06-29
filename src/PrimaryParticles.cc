#include "PrimaryParticles.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "time.h"
#include "G4Geantino.hh"
#include "G4IonTable.hh"
#include "G4ChargedGeantino.hh"
#include "G4NuclideTable.hh"
#include "G4VIsotopeTable.hh"
#include "G4Ions.hh"
#include <cmath> 

using namespace std;

PrimaryParticles::PrimaryParticles() 
{
    gun = new G4ParticleGun(1);
    gun->SetParticleEnergy(140 * keV);

    auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    gun->SetParticleDefinition(particleDefinition);
}

PrimaryParticles::~PrimaryParticles() 
{
  delete gun;
}

void PrimaryParticles::GeneratePrimaries(G4Event* event) 
{
  G4int seed = clock();


  G4double DirPhi = 2.0 * M_PI * G4UniformRand();
  G4double DirTheta = acos(1.0 - 2.0 * G4UniformRand());


  G4double Dirx = sin(DirTheta) * cos(DirPhi);
  G4double Diry = sin(DirTheta) * sin(DirPhi);
  G4double Dirz = cos(DirTheta);

  G4double R = 0;

  G4double PosPhi = 2.0 * M_PI * G4UniformRand();
  G4double PosTheta = acos(1.0 - 2.0 * G4UniformRand());  

  G4double Posx = R*sin(PosTheta)*cos(PosPhi);
  G4double Posy = R*sin(PosTheta)*sin(PosPhi);
  G4double Posz = R*cos(PosTheta);

  // gun->SetParticlePosition(G4ThreeVector(Posx, Posy, Posz));
  // gun->SetParticleMomentumDirection(G4ThreeVector(Dirx,Diry,Dirz));
  gun->SetParticlePosition(G4ThreeVector(0, 0, 0));
  gun->SetParticleMomentumDirection(G4ThreeVector(0,1.0,0.0));

  gun->GeneratePrimaryVertex(event);

}