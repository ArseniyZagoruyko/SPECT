#ifndef DETECTORSPHERE_HH
#define DETECTORSPHERE_HH

#include <G4VSensitiveDetector.hh>


class G4Step;

class DetectorSphere : public G4VSensitiveDetector 
{
public:

    DetectorSphere (const G4String& name);
    virtual ~DetectorSphere();

    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* );

    

private:

};


#endif