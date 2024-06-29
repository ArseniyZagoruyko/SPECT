class G4Run;

#include "G4Types.hh"

class G4UserRunAction
{
  public:
    G4UserRunAction();
    virtual ~G4UserRunAction();

  public:

    virtual void BeginOfRunAction(const G4Run* aRun);


};




