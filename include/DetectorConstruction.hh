#include <G4VUserDetectorConstruction.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4NistManager.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>
#include <G4AssemblyVolume.hh>
#include <G4SubtractionSolid.hh>
#include <G4SDManager.hh>
#include "SensitiveDetector.hh"

class SteppingAction;


class DetectorConstruction : public G4VUserDetectorConstruction
{
private:
    G4Box* world;
    G4LogicalVolume* log_world;
    G4VPhysicalVolume* phys_world;

    G4NistManager* nist;
    G4Material* world_mat;
    G4Material* water_mat;

    G4double distance1;
    G4double distance2;

    G4double radiatorThin1;
    G4double radiatorThin2;
    G4double radiatorLength;
    G4double radiatorHeight;

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    void ConstructWaterPhantom();

    void ConstructSphere();
    

    SensitiveDetector* sensitiveDetector1;
    SensitiveDetector* sensitiveDetector2;
    SensitiveDetector* sensitiveDetector3;

    G4LogicalVolume* thinVolumeLogical1;
    G4LogicalVolume* thinVolumeLogical2;
    G4LogicalVolume* radiatorLogical1;
    G4LogicalVolume* radiatorLogical2;
    G4LogicalVolume* log_Sphere;



    // G4SDManager для чувствительных детекторов 
    G4SDManager* SDman;

public:
    void ConstructScintillator();

    DetectorConstruction();
    ~DetectorConstruction();
};