#include "globals.hh"
#include "G4VUserParallelWorld.hh"
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4NistManager.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>

class MyParallelWorld : public G4VUserParallelWorld
{

public:

   G4NistManager* nist1;

   MyParallelWorld(G4String worldName);
   ~MyParallelWorld();

   virtual void Construct();
   virtual void ConstructSD();
};
