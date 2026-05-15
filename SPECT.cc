#include "DetectorConstruction.hh"
#include <QGSP_BIC.hh>
#include "LBE.hh"
#include <G4UImanager.hh>
#include <G4UIExecutive.hh>
#include <G4EmStandardPhysics.hh>
#include <G4EmLowEPPhysics.hh>
#include "Action.hh"
#include "G4UIcommand.hh"
#include "Randomize.hh"

#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4VisExecutive.hh"
#include "SensitiveDetector.hh"

#include <FTFP_BERT.hh>
#include "G4OpticalPhysics.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "MyParallelWorld.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include <unistd.h>
#include <ctime>

int main(int argc, char** argv)
{

    G4long seed = static_cast<G4long>(std::time(nullptr));
    CLHEP::RanecuEngine* engine = new CLHEP::RanecuEngine;
    engine->setSeed(seed);
    G4Random::setTheEngine(engine);

    G4RunManager* run = new G4RunManager;

    DetectorConstruction* geom = new DetectorConstruction();
    G4String paraWorldName = "ParallelWorld";
    geom->RegisterParallelWorld(new MyParallelWorld(paraWorldName));
    run->SetUserInitialization(geom);

    QGSP_BIC* pl = new QGSP_BIC;
    G4OpticalPhysics* OpticalPhysics = new G4OpticalPhysics;
    // pl->RegisterPhysics(OpticalPhysics);
    pl->RegisterPhysics(new G4ParallelWorldPhysics(paraWorldName));
    // pl->RegisterPhysics(new G4EmLowEPPhysics());
    run->SetUserInitialization(pl);

    Action* act = new Action;
    run->SetUserInitialization(act);

    G4VisExecutive* vis = new G4VisExecutive;
    vis->Initialize();
    run->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/run/setCut 0.01 mm");
    UImanager->ApplyCommand("/run/setCutForRegion RadiatorRegion 0.001 mm");
    run->Initialize();

    G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    UImanager->ApplyCommand("/control/verbose 0");
    UImanager->ApplyCommand("/run/verbose 0");

    bool lowelectro = false;

    if (lowelectro) {
        UImanager->ApplyCommand("/process/em/fluo true");
        UImanager->ApplyCommand("/process/em/auger true");
        UImanager->ApplyCommand("/process/em/augerCascade true");
        UImanager->ApplyCommand("/process/em/pixe true");
    }

    UImanager->ApplyCommand("/run/initialize");

    bool visual = true;

    if (visual) {
        UImanager->ApplyCommand("/vis/verbose 2");
        UImanager->ApplyCommand("/tracking/verbose 0");
        UImanager->ApplyCommand("/vis/scene/create");
        UImanager->ApplyCommand("/vis/open VRML2FILE");
        UImanager->ApplyCommand("/vis/drawVolume");
        UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 90. 90");
        UImanager->ApplyCommand("/vis/scene/add/trajectories");
        UImanager->ApplyCommand("/vis/scene/add/hits");
        UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByCharge");
        UImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2");
        UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
    }


    int number = 1000 + 0*static_cast<int>(G4UniformRand() * 1000);
    UImanager->ApplyCommand("/run/beamOn " + std::to_string(number));

    delete ui;
    delete vis;
    delete run;
}
