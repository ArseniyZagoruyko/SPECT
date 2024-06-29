#include "G4UserSteppingAction.hh"
#include <G4SystemOfUnits.hh>
#include <stdio.h>
#include <string>
#include <G4NistManager.hh>
#include <G4Material.hh>
#include <SensitiveDetector.hh>


using namespace std;


class DetectorConstruction;



class SteppingAction : public G4UserSteppingAction
{
public:
	SteppingAction(DetectorConstruction* detectorConstruction);
	virtual ~SteppingAction();

	virtual void UserSteppingAction(const G4Step* step);

private:
	
	ofstream outputFile;
	DetectorConstruction* fDetConstruction;

};

