#ifndef STEPPING_ACTION_HH
#define STEPPING_ACTION_HH

#include "G4UserSteppingAction.hh"
#include <G4SystemOfUnits.hh>
#include <fstream>
#include <string>
#include <G4NistManager.hh>
#include <G4Material.hh>
#include "G4UserEventAction.hh"
#include <vector>

using namespace std;

class DetectorConstruction;

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(DetectorConstruction* detectorConstruction);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step* step);

  
  void OpenNewBinaryFile();

  double x1;
  double y1;       
  double z1;
  double x2 = 0;
  double y2 = 0;       
  double z2 = 0;

  double Xp1 = 0;
  double Zp1 = 0;
  double Xp2 = 0;
  double Zp2 = 0;
   
  G4double angle;

  int tempId = 0;

  std::vector<int> eventIDs;


  float Rx1[100000000] = {0};
  float Ry1[100000000] = {0};
  float Rz1[100000000] = {0};
  float Rx2[100000000] = {0};
  float Ry2[100000000] = {0};
  float Rz2[100000000] = {0};

  float Tx1[100000000] = {0};
  float Ty1[100000000] = {0};
  float Tz1[100000000] = {0};
  float Tx2[100000000] = {0};
  float Ty2[100000000] = {0};
  float Tz2[100000000] = {0};

private:
  ofstream outputFile1;
  ofstream outputFile2;
  ofstream outputFile3;
  ofstream outputFile4;
  ofstream binaryOutputFile; 
  ofstream backscatterFile;
  ofstream spectreFile;


  //  для управления разделением файлов
  int eventsInCurrentFile = 0; // Счетчик событий в текущем файле
  int fileIndex = 0;           // Индекс текущего файла

  DetectorConstruction* fDetConstruction;
};

#endif