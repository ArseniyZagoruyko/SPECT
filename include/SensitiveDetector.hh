#ifndef SENSITIVEDETECTOR_HH
#define SENSITIVEDETECTOR_HH

#include <G4VSensitiveDetector.hh>
#include <fstream>

class G4Step;

class SensitiveDetector : public G4VSensitiveDetector 
{

public:

    SensitiveDetector(G4String name, G4int id);
    virtual ~SensitiveDetector();

    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* );

// геттер
    G4int GetDetectorID() const;

private:

// индификатор чувств детектора 
    G4int detectorID;
    
    std::ofstream hitsFile; 
};


#endif