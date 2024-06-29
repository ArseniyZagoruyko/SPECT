#include "MyParallelWorld.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include <G4SystemOfUnits.hh>
#include<G4SDManager.hh>
#include <G4NistManager.hh>
#include <G4PVReplica.hh>
#include "G4VSolid.hh"
#include "G4Torus.hh"
#include "G4SubtractionSolid.hh"

using namespace std;


MyParallelWorld::MyParallelWorld(G4String worldName):G4VUserParallelWorld(worldName)
{


}
MyParallelWorld::~MyParallelWorld()
{;}


void MyParallelWorld::Construct()
{

    G4VPhysicalVolume* ghostWorld = GetWorld();
    G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

}


void MyParallelWorld::ConstructSD()
{





}
