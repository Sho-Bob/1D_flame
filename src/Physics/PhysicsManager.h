#ifndef PHYSICSMANAGER_H
#define PHYSICSMANAGER_H

#include "pengrobinson.h"

#include "IO/input.h"

class Physics {
  public:
    Physics();
    ~Physics();
    void initialize();

    void set_input_file(Input* input_file) { this->input = input_file; }

    // \brief: set mixture state using T, P, Y; used in initialization
    void SetMixture_TPY(const double T, const double p, const double* Y) {
      this->eos->SetMixture_TPY(T, p, Y);
    }

    void SetMixture_PRY(const double p, const double rho, const double* Y, const double T_guess = 300.0, const bool flash = true) {
      this->eos->SetMixture_PRY(p, rho, Y, T_guess, flash);
    }

    void SetMixture_ERY(const double e, const double rho, const double* Y, const double T_guess = 300.0, const double P_guess = 1e5, const bool flash = true) {
      this->eos->SetMixture_ERY(e, rho, Y, T_guess, P_guess, flash);
    }

    double GetRho() { return *rho; }
    double GetT() { return *T; }
    double GetP() { return *p; }
    double GetE() { return *e; } // specific internal energy
    double GetSos() { return *sos; } // speed of sound

  private:
    PengRobinson* eos; // for now
    std::vector<std::string> species_names;
    Input* input;

    // thermodynamic state
    double* T;
    double* p;
    std::vector<double>* Y; // mass fractions
    double* X; // mole fractions
    double* rho;
    double* e; // specific internal energy
    double* sos; // speed of sound


};

#endif // PHYSICSMANAGER_H
