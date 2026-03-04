#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include <fstream>
#include <vector>

class PrimaryParticles : public G4VUserPrimaryGeneratorAction{
public:

    PrimaryParticles();
    ~PrimaryParticles();
    int  num = 0;

    virtual void GeneratePrimaries(G4Event* event);

private:

    std::vector<double> energySpectrum; // Спектр энергий
    double CalculateEnergy(int channel); // Функция для расчета энергии по каналу
    G4ParticleGun* gun;
};