#include "SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "DetectorConstruction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include <fstream>
#include <iostream>
#include <G4SystemOfUnits.hh>
#include <string>
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include <vector>
#include <cmath>

using namespace std;

vector<float> R1(100000000, 0.0);
vector<float> R2(100000000, 0.0);
// vector<float> T1(100000000, 0.0);
// vector<float> T2(100000000, 0.0);

vector<pair<int, int>> photonCoords1;
vector<pair<int, int>> photonCoords2;



// vector<vector<string>> interactionHistory(100000000, vector<string>());

// счетчик для отслеживания количества событий в текущем файле
int eventsInCurrentFile = 0;
int fileIndex = 0; // индекс текущего файла

SteppingAction::SteppingAction(DetectorConstruction* detectorConstruction)
    : G4UserSteppingAction(),
      fDetConstruction(detectorConstruction)
{
    outputFile1.open("deposits_FirstDet.txt", ios::app);
    if (!outputFile1.is_open()) {
        G4cerr << "Error opening deposits.txt!" << G4endl;
    }

    outputFile2.open("GammaCoordinates_FirstDet.txt", ios::app);
    if (!outputFile2.is_open()) {
        G4cerr << "Error opening GammaCoordinates.txt!" << G4endl;
    }

    // outputFile3.open("deposits_SecondDet.txt", ios::app);
    // if (!outputFile3.is_open()) {
    //     G4cerr << "Error opening deposits.txt!" << G4endl;
    // }

    // outputFile4.open("GammaCoordinates_SecondDet.txt", ios::app);
    // if (!outputFile4.is_open()) {
    //     G4cerr << "Error opening GammaCoordinates.txt!" << G4endl;
    // }

    // backscatterFile.open("backscattering_events.txt", ios::app);
    // if (!backscatterFile.is_open()) {
    //     G4cerr << "Error opening backscattering_events.txt!" << G4endl;
    // }

    // spectreFile.open("spectreFile.txt", ios::app);
    // if (!spectreFile.is_open()){
    //     G4cerr<<"Error open spectre"<<G4endl;
    // }


    // OpenNewBinaryFile();
}

// void SteppingAction::OpenNewBinaryFile() {

//     if (binaryOutputFile.is_open()) {
//         binaryOutputFile.close();
//     }


//     string fileName = "photon_coordinates_" + to_string(fileIndex) + ".bin";
//     binaryOutputFile.open(fileName, ios::binary | ios::app);
//     if (!binaryOutputFile.is_open()) {
//         G4cerr << "Error opening " << fileName << "!" << G4endl;
//     }

//     // Сбрасываем счетчик событий
//     eventsInCurrentFile = 0;
// }

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();

    int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    G4LogicalVolume* currentVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4StepPoint* prePoint = step->GetPreStepPoint();
    G4TouchableHandle touch1 = prePoint->GetTouchableHandle();
    G4VPhysicalVolume* volume = touch1->GetVolume();
    G4String vname = volume->GetName();
    G4double energyDeposit = step->GetTotalEnergyDeposit() / keV;
    G4ThreeVector position = step->GetPostStepPoint()->GetPosition();

    G4StepPoint* postPoint = step->GetPostStepPoint();
    G4VPhysicalVolume* preVolume = prePoint->GetTouchableHandle()->GetVolume();
    G4VPhysicalVolume* postVolume = postPoint->GetTouchableHandle()->GetVolume();

    G4String preVolumeName = preVolume ? preVolume->GetName() : "";
    G4String postVolumeName = postVolume ? postVolume->GetName() : "";

    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();



    // ////Backscattering
    
    // if (particleName == "gamma") {
    //     if (vname.contains("radiator1")) {
    //         interactionHistory[eventID].push_back("radiator1");
    //     } else if (vname.contains("radiator2")) {
    //         interactionHistory[eventID].push_back("radiator2");
    //     }
    // }

    // if (interactionHistory[eventID].size() >= 3) {
    //     // последовательность: radiator1 -> radiator2 -> radiator1
    //     if (interactionHistory[eventID][interactionHistory[eventID].size() - 3] == "radiator1" &&
    //         interactionHistory[eventID][interactionHistory[eventID].size() - 2] == "radiator2" &&
    //         interactionHistory[eventID][interactionHistory[eventID].size() - 1] == "radiator1") {
    //         backscatterFile << "Event ID: " << eventID << "\tBackscattering" << G4endl;
    //     }
    // }

    // ////





    if (vname.contains("radiator1") && track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
        R1[eventID] += energyDeposit; 
    }

    if (vname.contains("radiator2") && track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
        R2[eventID] += energyDeposit;
    }

    // if (vname.contains("radiator3") && track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
    //     T1[eventID] += energyDeposit;
    // }

    // if (vname.contains("radiator4") && track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
    //     T2[eventID] += energyDeposit;
    // }



    if (vname.contains("radiator1") && particleName == "gamma") {
        Rx1[eventID] = position.x() / mm;
        Ry1[eventID] = position.y() / mm;
        Rz1[eventID] = position.z() / mm;
    }

    if (vname.contains("radiator2") && particleName == "gamma") {
        Rx2[eventID] = position.x() / mm;
        Ry2[eventID] = position.y() / mm;
        Rz2[eventID] = position.z() / mm;
    }

    if (vname.contains("radiator3") && particleName == "gamma") {
        Tx1[eventID] = position.x() / mm;
        Ty1[eventID] = position.y() / mm;
        Tz1[eventID] = position.z() / mm;
    }

    if (vname.contains("radiator4") && particleName == "gamma") {
        Tx2[eventID] = position.x() / mm;
        Ty2[eventID] = position.y() / mm;
        Tz2[eventID] = position.z() / mm;
    }

    // потом доделать для опт фотонов на 2 паре детекторов 

    // if (vname.contains("radiator1") && position.y() == 104 && particleName == "opticalphoton") {
    //     double Px1 = position.x() / mm;
    //     double Pz1 = position.z() / mm;

    //     double cellSize = 5.0; 
    //     Px1 = round(Px1 / cellSize) * cellSize;
    //     Pz1 = round(Pz1 / cellSize) * cellSize;

    //     photonCoords1.push_back(make_pair(Px1, Pz1));
    // }

    // if (vname.contains("radiator2") && position.y() == 138 && particleName == "opticalphoton") {
    //     double Px2 = position.x() / mm;
    //     double Pz2 = position.z() / mm;

    //     double cellSize = 5.0; 
    //     Px2 = round(Px2 / cellSize) * cellSize;
    //     Pz2 = round(Pz2 / cellSize) * cellSize;

    //     photonCoords2.push_back(make_pair(Px2, Pz2));
    // }


    // if (vname.contains("phys") && particleName == "opticalphoton")
    // {
    //     track->SetTrackStatus(fStopAndKill);
    // }
    

    if (eventID != tempId) {



        // /// SPECTRE

        // spectreFile << eventID << "\t" << R1[eventID - 1] << "\t" << R2[eventID - 1] << endl;

        // /////




        //первая пара детекторов
        if ((R1[eventID - 1] * R2[eventID - 1]) != 0) {

            outputFile2 << eventID << "\t" << Rx1[eventID - 1] << "\t" << Ry1[eventID - 1] << "\t" << Rz1[eventID - 1]
                        << "\t" << Rx2[eventID - 1] << "\t" << Ry2[eventID - 1] << "\t" << Rz2[eventID - 1] << endl;
            outputFile1 << eventID << "\t" << R1[eventID - 1] << "\t" << R2[eventID - 1] << endl;

            //Запись оптических фотонов

            // for (const auto& coord : photonCoords1) {
            //     int radiatorID = 1;
            //     binaryOutputFile.write(reinterpret_cast<const char*>(&eventID), sizeof(int));       // Event ID
            //     binaryOutputFile.write(reinterpret_cast<const char*>(&radiatorID), sizeof(int));    // Radiator ID
            //     binaryOutputFile.write(reinterpret_cast<const char*>(&coord.first), sizeof(double)); // X coordinate
            //     binaryOutputFile.write(reinterpret_cast<const char*>(&coord.second), sizeof(double)); // Z coordinate
            // }

            // for (const auto& coord : photonCoords2) {
            //     int radiatorID = 2;
            //     binaryOutputFile.write(reinterpret_cast<const char*>(&eventID), sizeof(int));       // Event ID
            //     binaryOutputFile.write(reinterpret_cast<const char*>(&radiatorID), sizeof(int));    // Radiator ID
            //     binaryOutputFile.write(reinterpret_cast<const char*>(&coord.first), sizeof(double)); // X coordinate
            //     binaryOutputFile.write(reinterpret_cast<const char*>(&coord.second), sizeof(double)); // Z coordinate
            // }

            // // Увеличиваем счетчик событий
            // eventsInCurrentFile++;

            // if (eventsInCurrentFile >= 100) {
            //     fileIndex++;
            //     OpenNewBinaryFile();
            // }
        }

        //вторая пара детекторов

        // if ((T1[eventID - 1] * T2[eventID - 1]) != 0) {


        //     outputFile4 << eventID << "\t" << Tx1[eventID - 1] << "\t" << Ty1[eventID - 1] << "\t" << Tz1[eventID - 1]
        //     << "\t" << Tx2[eventID - 1] << "\t" << Ty2[eventID - 1] << "\t" << Tz2[eventID - 1] << endl;
        //     outputFile3 << eventID << "\t" << T1[eventID - 1] << "\t" << T2[eventID - 1] << endl;



        // }

        photonCoords1.clear();
        photonCoords2.clear();
    }



    tempId = eventID;
}

SteppingAction::~SteppingAction()
{   
    if (outputFile1.is_open()) {
        outputFile1.close();
    }

    if (outputFile2.is_open()) {
        outputFile2.close();
    }

    // if (outputFile3.is_open()) {
    //     outputFile3.close();
    // }

    // if (outputFile4.is_open()) {
    //     outputFile4.close();
    // }

    // if (binaryOutputFile.is_open()) {
    //     binaryOutputFile.close();
    // }

    

    // if (backscatterFile.is_open()) {
    //     backscatterFile.close();
    // }

    // if (spectreFile.is_open())
    // {
    //     spectreFile.close();
    // }
    
}