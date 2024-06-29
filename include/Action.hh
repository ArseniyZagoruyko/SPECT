#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class Action : public G4VUserActionInitialization{

    public:
        Action();
        ~Action();

        virtual void Build() const;

    private:
        DetectorConstruction* fDetConstruction;

};