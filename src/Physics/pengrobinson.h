#ifndef PENGROBINSON_H
#define PENGROBINSON_H

#include <vector>
#include "Common/common.h"

class PengRobinson {
 public:
  PengRobinson();
  ~PengRobinson();

  virtual void init(const std::vector<std::string>& speciesNames, const double Pamb, std::string mixing_rule= "MR3");

  // Friend insertion to ostream
  friend std::ostream& operator<<(std::ostream& os, const PengRobinson& thermal);

  // These functions will NOT result in a consistent class state for all member variables
  //! \brief set mixture based on T, P, and rhoYF; used in initialization
  void SetMixture_TPRYF(const double T, const double P, const double rhoYF, const double YFGuess= 0.5, const bool flash=false);
  virtual double GetRho_SetMixture_TPY(const double T, const double P, const double* Y);
  virtual double GetRho_SetMixture_TPX(const double T, const double P, const double* X);
  double GetRho_SetMixture_TPRY(const double T, const double P, double RYF, double rhoGuess);
  double helperFun(double T, double R, double RYF);

  void GetTdew_PXvapor(const double P, const double* Xvapor, const double TGuess, const double* XliqGuess, double& Tdew, double* Xliq);
  void GetPsat_TXvapor(const double T, const double* Xvapor, const double PGuess, const double* XliqGuess, double& Psat, double* Xliq);
  virtual void TPFlash(const double T, const double P, const double* XvGuess, const double* XlGuess, double* Xv, double* Xl, double& Tout);
  // \brief: Flash calculation for two-phase equilibrium, force to get the two-phase and return the lowest value of X[0] that gives a two-phase
  virtual double TPFlashForce(const double T, const double P, double* Xv, double* Xl, double& Tout);

  // These functions will result in a consistent class state for all member variables
  virtual void SetMixture_TRY(const double T, const double R, const double* Y);
  virtual void SetMixture_TRX(const double T, const double R, const double* X);

  // @brief: Temperature, pressure and mass fraction; typically used in initial conditions
  void SetMixture_TPY(const double T, const double P, const double* Y);
  // @brief: (mass-based) energy, density, and species; used after conservative variables are set
  void SetMixture_ERY(const double e, const double rho, const double* Y, const double T_guess = 300.0, const double P_guess = 1e5, const bool flash = false);
  // @brief: pressure, density, and species; used after using quasi-conservative methods such as double-flux
  void SetMixture_PRY(const double P, const double rho, const double* Y, const double T_guess = 300.0, const bool flash = true);
  void SetMixture_EKorRY(const double E, const double R, const double *Y, const double *gradRhogradRho, const double T_guess);

  void SetMixture_HPY(const double H, const double P, const double* Y, const double T_guess);

  inline void SetMolarFraction(const double *scalars_in) {
    // set molar fraction from input
    double M = 0.0;
    LOOP_l_N(this->num_species) M += scalars_in[l] * this->MW[l];
    LOOP_l_N(this->num_species) {
      this->X[l] = scalars_in[l];
      this->Y[l] = scalars_in[l] * this->MW[l] / M;
    }
    this->MW_M = M;
  }

  inline void SetMolarFraction() {
    LOOP_l_N(this->num_species) this->X[l] =
                                    this->Y[l] * this->MW_M / this->MW[l];
  };
  inline void SetMassFraction(const double *scalars_in) {
    //LOOP_l_N(this->num_species) this->Y[l] = scalars_in[l];
    this->Y[this->num_species - 1] = 1.0;
    LOOP_l_N(this->num_species - 1) {
      //this->Y[l] = max(min(scalars_in[l] * rho_in, rho_in), 0.0) / rho_in;
      //this->Y[l] = max(min(scalars_in[l], 1.0), 0.0); // we clip here
      this->Y[l] = scalars_in[l];
      this->Y[this->num_species - 1] -= this->Y[l];
    }
  };
  // Set Pamb
  inline void SetPamb(const double Pamb) { this->Pamb = Pamb; }

  // Get methods
  inline double GetUniversalGasCosntant() const { return this->gas_constant_universal; }

  inline int GetNSpecies() const { return this->n_species; }
  inline const std::vector<double>& GetMW() const { return this->MW; }
  inline const std::vector<double>& GetMolarDensities() const { return this->molar_density; }

  inline double GetT() const { return this->T; }
  inline double* GetT_ptr() { return &T; }
  inline double GetP() const { return this->P; }
  inline double* GetP_ptr() { return &P; }
  inline double GetPGD() const { return this->P_GD; }
  inline double GetRho() const { return this->rho; }
  inline double* GetRho_ptr() { return &rho; }
  inline double* GetE_ptr() { return &E_mass; }
  inline double GetZ() const { return this->Z; }
  inline double GetBv() const { return this->expansivity; }
  inline double GetCp(const bool molar=false) const { 
    if (molar) {
      return this->cp; 
    } else {
      return this->cp / this->MW_M; // J/kg/K
    }
  }
  inline const std::vector<double> GetY() const { return this->Y; }
  inline std::vector<double>* GetY_ptr() { return &Y; }
  inline const std::vector<double> GetX() const { return this->X; }
  inline std::vector<double> GetCSpec() { return this->c_spec; };
  inline const double GetX_TP(const int ids) const { return this->x_TP[ids]; }
  inline const double GetY_TP(const int ids) const { return this->y_TP[ids]; }

  inline double GetMWm() const { return this->MW_M; }
  inline double GetInternalEnergyPhase(const int phase = 0) const {
    if (phase == this->iVAPPH_TP) {
      return this->E_V; // J/mol
    } else if (phase == this->iLIQPH_TP) {
      return this->E_L; // J/mol
    } else {
      return this->E; // J/mol
    }
  }

  inline double GetDepartureInternalEnergy(const int phase = 0) const {
    if (phase == this->iVAPPH_TP) {
      return this->depE_V; // J/mol
    } else if (phase == this->iLIQPH_TP) {
      return this->depE_L; // J/mol
    } else {
      return 0.0; // No departure energy for mixture
    }
  }

  inline double GetDepartureInternalEnergydT(const int phase = 0) const {
    if (phase == this->iVAPPH_TP) {
      return this->dDepEdT_V;
    } else if (phase == this->iLIQPH_TP) {
      return this->dDepEdT_L;
    } else {
      return 0.0; // No departure energy for mixture
    }
  }



  inline double GetInternalEnergy(const bool molar = false) const { 
    if (molar) {
      return this->E; // J/kmol
    } else {
      return this->E_mass; // J/kg
    }
  }
  inline double GetdEdT(const int phase = 0) const { 
    if (phase == this->iVAPPH_TP) {
      return this->dEdT_V; // J/mol/K
    } else if (phase == this->iLIQPH_TP) {
      return this->dEdT_L; // J/mol/K
    } else {
      return this->dEdT; 
    }
  }

  inline double GetdEdp(const int phase = 0) const { 
    if (phase == this->iVAPPH_TP) {
      return this->dEdp_V; // J/mol/K
    } else if (phase == this->iLIQPH_TP) {
      return this->dEdp_L; // J/mol/K
    } else {
      return this->dEdp; 
    }
  }
  inline void GetdEdX(double* dEdX, const int phase = 0) { 
    if (phase == this->iVAPPH_TP) {
      LOOP_k_N(this->num_species) dEdX[k] = this->dEdchi_V[k];
    } else if (phase == this->iLIQPH_TP) {
      LOOP_k_N(this->num_species) dEdX[k] = this->dEdchi_L[k];
    } else {
      LOOP_k_N(this->num_species) dEdX[k] = 0.0; // No derivative for mixture
    }
  }

  inline double GetBeta() const { return this->betaV; }
  inline double GetInternalEnergyKor(const double *gradRhogradRho) {
    double EKorOnly = 0.;
    LOOP_k_N(this->n_species) {
      LOOP_l_N(this->n_species) {
        int apos = k * this->n_species + l;
        EKorOnly += 0.5 * this->kappa[apos] * gradRhogradRho[apos] / this->MW[k] / this->MW[l];
      }
    }
    EKorOnly /= this->rho;
    return this->E + EKorOnly;
  }
  double GetIdealGasCp(double T_in, const bool molar = false);
  double GetIdealGasEnthalpy(double T_in, const bool molar = false);
  double GetEnergyFromTemperature(const double T_in, const bool molar = false);
  inline const std::vector<double> GetPartialEnthalpies() const { return this->hspecies; }
  inline const std::vector<double> GetPartialChemEnergies() const { return this->muspecies; }
  inline const std::vector<std::vector<double>> GetPartialChemEnergiesGradient() const { return this->dmudn_species; }
  inline const std::vector<double> GetPartialEntropies() const { return this->sspecies; }
  inline const std::vector<double> GetMassDiffusionFactors() const { return this->alphaD_IJ; }
  inline const std::vector<double> GetThermalDiffusionFactors() const { return this->alphaT; }
  inline const std::vector<double> GetCapillaryCoefficient() const { return this->kappa; }
  inline const std::vector<double> GetBinDiffusivity() const { return this->Dm; }
  inline double GetViscosity() const { return this->mu; }
  inline double GetThermalConducvitity() const { return this->lambda; }
  inline const std::vector<double> GetLij() const { return this->Lij; }
  inline const std::vector<double> GetLqi() const { return this->Lqi; }
  inline double getSos() const { return this->sos; }
  inline double* getSos_ptr() { return &sos; }
  inline double getSoS2() const { return this->sos2; }
  // Estimate interface thickness and time scale
  inline void GetInterfaceThicknessEstimate(double* estimates) const {
    LOOP_k_N(this->n_species) {
      LOOP_l_N(this->n_species) {
        int apos = k * this->n_species + l;
        estimates[apos] = std::sqrt(this->kappa[apos] / this->Am);
      }
    }
  }

  inline void GetTimeScaleEstimate(double* estimates) const {
    LOOP_k_N(this->n_species) {
      LOOP_l_N(this->n_species) {
        int apos = k * this->n_species + l;
        estimates[apos] = std::sqrt(this->kappa[apos] / this->Am)
            / std::sqrt(this->Am / this->Bm / this->MW_M);
      }
    }
  }

  double getdmuF_muOdrhoYF();
  void getdmuidrhoj(double ** dmuidrhoj);
  void getdmuidnj(double ** dmuidnj);
  virtual void GetRhoLRhoV(const double T, const double P, double & , double &);
  inline void get_rhoYFL_rhoYFV(double& rho_yfl, double& rho_yfv) {
    rho_yfl = this->rhoYF_L;
    rho_yfv = this->rhoYF_V;
  }
  inline void get_rhoL_rhoV(double& rho_l, double& rho_v) {
    rho_l = this->rho_L;
    rho_v = this->rho_V;
  }
  double GetRhoFromN(const double *);
  double GetHelmholtzEnergy() {return this->helmholtz_energy;};
  inline double GetPhi() { return this->Phi; };
  inline double GetZPhi() const { return this->Z_phi; }
  inline double GetPhiRho() const { return this->phi_rho; }

  inline void GetdChidT(const int phase, double* dchidT) {
    if (phase == this->iLIQPH_TP) {
      LOOP_k_N(this->num_species) dchidT[k] = this->dxdT[k];
    } else if (phase == this->iVAPPH_TP) {
      LOOP_k_N(this->num_species) dchidT[k] = this->dydT[k];
    } else {
      throw std::runtime_error("Invalid phase for dChidT");
    }
  }

  inline void GetdChidp(const int phase, double* dchidp) {
    if (phase == this->iLIQPH_TP) {
      LOOP_k_N(this->num_species) dchidp[k] = this->dxdp[k];
    } else if (phase == this->iVAPPH_TP) {
      LOOP_k_N(this->num_species) dchidp[k] = this->dydp[k];
    } else {
      throw std::runtime_error("Invalid phase for dChidp");
    }
  }

 protected:
  void ReadNasaPolynomials();
  void ReadCriticalProperties();
  void SetRealFluidConstants();

  inline void SetMassFractionFromY(const double* Y) {
    this->Y[this->n_species - 1] = 1.;
    LOOP_l_N(this->n_species - 1) {
      this->Y[l] = Y[l];
      this->Y[this->n_species - 1] -= this->Y[l];
    }
  }

  inline void SetMolecularWeightMixture() { this->SetMolecularWeightMixtureFromY(); }
  inline void SetMolecularWeightMixtureFromY() {
    this->MW_M = 0.;
    LOOP_l_N(this->n_species)
        this->MW_M += this->Y[l] / MW[l];
    this->MW_M = 1. / this->MW_M;

    this->gas_constant = this->gas_constant_universal / this->MW_M;
  }

  inline void SetMolarFractionFromY() {
    LOOP_l_N(this->n_species)
        this->X[l] = this->Y[l] * this->MW_M / this->MW[l];
  }

  inline void SetMolarFractionFromX(const double* X) {
    this->X[this->n_species - 1] = 1.;
    LOOP_l_N(this->n_species - 1) {
      this->X[l] = X[l];
      this->X[this->n_species - 1] -= this->X[l];
    }
  }

  inline void SetMolecularWeightMixtureFromX() {
    this->MW_M = 0.;
    LOOP_l_N(this->n_species)
        this->MW_M += this->X[l] * this->MW[l];

    this->gas_constant = this->gas_constant_universal / this->MW_M;
  }

  inline void SetMassFractionFromX() {
    LOOP_l_N(this->n_species)
        this->Y[l] = this->X[l] * this->MW[l] / this->MW_M;
  }

  void SyncRealFluidThermodynamicsFromTemperatureDensity();
  void SyncIdealFluidThermodynamicsFromTemperature();

  void SyncRZFromPressureTemperature();
  void SyncPZFromTemperatureDensity();
  void SyncEFromTemperatureDensity();
  void SyncPartialPropertiesFromTemperatureDensity();
  void SyncTransportFromTemperatureDensity();
  void SyncMolarDensities();
  virtual void SyncHelmholtzEnergy();

  // Transport models
  void SetTransportConstants();

  // For iterating to find T from E
  void PrepareRealFluidThermodynamicsForTFromE();
  double GetTemperatureFromEnergyKor(const double E_in, const double *gradRhogradRho, double T_guess);
  double EFromTFunc(double T, const double *gradRhogradRho);

  // Utilities
  // Test if current T and P can give more than one solution
  bool CheckThreeRoots();
  void ComputeMixtureCrit(double& Tcmix, double& Pcmix, double& wcmix, double& Zcmix);

  // For dew point calulcation only
  void getPhifromTandPandX(const double T_in, const double P_in, const double* X_in, double* phiOut);
  virtual void getPhifromTandP_SS(const double T_in, const double P_in, double* phiOut) {};

 protected:
  double gas_constant_universal;
  double N_Av;
  double Boltzman;

  double delta1, delta2, delta1_p_delta2, delta1_m_delta2, delta1delta2;

  int n_species;
  int num_species;
  std::vector<std::string> species;
  double Pamb;

  double T;
  double P;
  double P_GD;
  double PForDiff;
  double MW_M;
  double gas_constant;

  double rho;
  double Z;
  double ZForDiff;

  double E; // J/mol, per unit mole, to be consistent with the literature
  double E_mass; // J/kg

  double cp, cv, gamma;
  double sos;
  double sos2;

  std::vector<double> A_IJ;
  double Am;
  double Bm;
  double dAmdT;
  double d2AmdT2;
  std::vector<double> dAmdN;
  std::vector<double> d2AmdTdN;
  double K1;
  std::vector<double> dK1dN;

  double BT;
  double betaT;
  double expansivity;
  double dPdT;
  double dPdV;
  std::vector<double> dPdN;
  std::vector<double> dVdN;

  std::vector<double> hspecies;
  std::vector<double> sspecies;
  std::vector<double> muspecies;
  std::vector<double> dmudn;
  std::vector<std::vector<double>> dmudn_species;
  std::vector<double> muspecies_molar;
  std::vector<double> edspecies;

  std::vector<double> alphaD_IJ;

  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> MW;
  std::vector<double> molar_density;
  double molar_density_total;

  std::vector<double> Tcrit;
  std::vector<double> Pcrit;
  std::vector<double> rhocrit;
  std::vector<double> Vcrit;
  std::vector<double> Zcrit;
  std::vector<double> omega;
  std::vector<double> sigma;
  std::vector<double> dipole;
  std::vector<double> epsOverk;

  std::string mixing_rule;

  int NasaCoef;
  std::vector<double> nasa_poly_coeff;
  std::vector<std::pair<double, double>> nasa_poly_bounds;
  std::vector<double> cpspecies_ig;
  std::vector<double> hspecies_ig;
  std::vector<double> sspecies_ig;
  double h_ig;
  double cp_ig;
  double cv_ig;

  std::vector<double> Tcrit_IJ;
  std::vector<double> Vcrit_IJ;
  std::vector<double> Zcrit_IJ;
  std::vector<double> Pcrit_IJ;
  std::vector<double> omega_IJ;

  std::vector<double> cst_a;
  std::vector<double> cst_b;
  std::vector<double> cst_c; // k

  // Transport models
  std::vector<double> sigma_IJ;
  std::vector<double> epsOverk_IJ;
  std::vector<double> omega_visc_IJ;
  std::vector<double> MW_IJ;

  std::vector<double> Dm;
  double mu;
  double lambda;
  std::vector<double> alphaT;

  std::vector<double> m_IJ;
  std::vector<double> sigma_diff_IJ;

  // Onsager
  std::vector<double> Lij;
  std::vector<double> Lqi;
  std::vector<double> kT;

  // Capillary coefficient
  std::vector<double> k0;
  std::vector<double> k1;
  std::vector<double> k2;
  std::vector<double> Akappa;
  std::vector<double> Bkappa;
  std::vector<double> kappa; // J m^5 kmol^-2
                             //
  // Helmholtz
  double helmholtz_energy;
  double Phi;
  std::vector<double> kappa_spec;
  std::vector<double> c_spec;
 public:
  // additional real fluid methods
  double GetAm(const double T_in, const bool compute_gradients = false);
  double GetBm();
  // \brief: compute dAmdT, d2AmdT2, and other quantities that are needed for them
  void computeBasicDerivatives(const int phase = 0);
  // \brief: compute fugacity coefficient and fugacities
  void computeFugacityCoefficient(double *, double *, int phase = 0);
  // TWO-PHASE METHODS
  // \brief: compute densities, energies, and fugacities of a phase: Algorithm 1
  void generateVLstructure(const double T, const double P, const double* chi, const int phase);
  // \brief: VLE calculation
  int VLE(const double T_in, const double P_in, const double* X_in,
      double& beta_out, double* x_out, double* y_out);
  // \brief: SSI to get the phase equilibrium
  void SSI(const double T_in, const double P_in, const double* X_in,
      const double* K_in, double& beta_out, double* x_out, double* y_out, int& flag,
      const double eps = 1e-6, const int max_iter = 1000);
  // \brief: StabilityAnalysis
  void StabilityAnalysis(const double T_in, const double P_in, const double* X_in,
      double* K, const double* d_in, int& nphases, const double eps = 1e-6);
  // \brief: UpdateTrialComposition
  int UpdateTrialComposition(const double T_in, const double P_in, const double* X_in,
      const double* d, int phase, double* chi, 
      double& tpd_phase,
      const double eps_tpd=1e-6);
  // \brief: Get two-phase properties based on e and rho
  void solveERHO(const double e_in, const double rho_in, const double* X_in,
      double& T_out, double& P_out, double& beta_out, double* x_out, double* y_out,
      const double eps_ERHO = 1e-8, const double deltaT_max = 1.0, const double deltaP_max = 1.0e3,
      const double Nmax = 1000,
      const double T_0 = 300.0, const double P_0 = 1.0e5, const double lambda_0 = 1.0, const bool accelerate = true);
  // \brief: Get two-phase properties based on e and rho, but with mass fractions
  void solveERHOY(const double e_in, const double rho_in, const double* Y_in,
      double& T_out, double& P_out, double& beta_out, double* x_out, double* y_out,
      const double eps_ERHO = 1e-8, const double deltaT_max = 1.0, const double deltaP_max = 1.0e3,
      const double Nmax = 1000,
      const double T_0 = 300.0, const double P_0 = 1.0e5, const double lambda_0 = 1.0, const bool accelerate = true);
  // \brief: Get two-phase properties based on P and rho
  void solvePRHO(const double P_in, const double rho_in, const double* X_in,
      double& T_out, double& beta_out, double* x_out, double* y_out,
      const double eps_PRHO = 1e-8, const double deltaT_max = 1.0, const double N_max = 1000,
      const double T_0 = 300.0, const double lambda_0 = 1.0);
  // \brief: Get two-phase properties based on P and rho, but with mass fractions
  void solvePRHOY(const double P_in, const double rho_in, const double* Y_in,
      double& T_out, double& beta_out, double* x_out, double* y_out,
      const double eps_PRHO = 1e-8, const double deltaT_max = 1.0, const double N_max = 1000,
      const double T_0 = 300.0, const double lambda_0 = 1.0);
  // \brief: Solve VLE problem and set energies
  void solveTP(const double T_in, const double P_in, const double* X_in, const bool calc_derivatives = false);
  // \brief:basically like above but with mass fractions
  void solveTPY(const double T_in, const double P_in, const double* Y_in, const bool calc_derivatives = false);
  // \brief: compute additional derivatives for PRHO and ERHO problems
  void computeAdditionalDerivatives(const double T_in, const double P_in, const double* X_in,
      const double beta_in, const double* x_in, const double* y_in);
  // \brief: variations in T, P, x, y
  void computeTPvariations(const double T_in, const double P_in, const double* X_in,
      const double beta_in, const double* x_in, const double* y_in);
  // \brief: compute d_logPhi*
  void computeFugacityDerivatives(const double* chi_in, const int phase);
  // \brief: compute phase derivatives of phase
  void computeAdditionalPhaseDerivatives(const double T_in, const double P_in, const double *X_in, const int phase, const double * chi_in);
  double GetDensityFromPressureTemperature_TP(const double p_in,
                                                   const double T_in, const int phase = 0);
  inline double* GetFugacityCoefficient(int phase = 0) {
    if (phase == this->iVAPPH_TP) {
      return this->phi_V;
    } else if (phase == this->iLIQPH_TP) {
      return this->phi_L;
    } else {
      return this->phi_fugacity;
    }
  }

  inline void SyncZPhi() { 
    if (this->betaV <= 0.0) {
      this->Z_phi = 0.0;
      this->phi_rho = 0.0;
    } else if (this->betaV >= 1.0) {
      this->Z_phi = 1.0;
      this->phi_rho= 1.0;
    } else {
      this->Z_phi = (this->rho * Y[0] - this->rhoYF_V) / (this->rhoYF_L - this->rhoYF_V); 
      this->phi_rho = (this->rho - this->rho_V) / (this->rho_L - this->rho_V);
    }
  }

  double solveRachfordRice(const double* X_in, const double* x_in, const double* y_in);
  inline void GetFugacityDerivatives(const int phase, double* d_logPhi_dX, double* d_logPhi_dT, double* d_logPhi_dp) {
    if (phase == this->iVAPPH_TP) {
      LOOP_k_N(this->num_species) {
        LOOP_l_N(this->num_species) {
          int apos = k * this->num_species + l;
          d_logPhi_dX[apos] = this->d_logPhiV_dy[apos];
        }
        d_logPhi_dT[k] = this->d_logPhiV_dT[k];
        d_logPhi_dp[k] = this->d_logPhiV_dp[k];
      }
    } else if (phase == this->iLIQPH_TP) {
      LOOP_k_N(this->num_species) {
        LOOP_l_N(this->num_species) {
          int apos = k * this->num_species + l;
          d_logPhi_dX[apos] = this->d_logPhiL_dx[apos];
        }
        d_logPhi_dT[k] = this->d_logPhiL_dT[k];
        d_logPhi_dp[k] = this->d_logPhiL_dp[k];
      }
    } else {
      throw std::runtime_error("Invalid phase for fugacity derivatives");
    }
  }

  inline int GetLiquidPhaseFlag() const { return this->iLIQPH_TP; }
  inline int GetVaporPhaseFlag() const { return this->iVAPPH_TP; }
  inline void GetAmStuff(double& Am, double& dAmdT, double& d2AmdT2) {
    Am = this->Am;
    dAmdT = this->dAmdT;
    d2AmdT2 = this->d2AmdT2;
  }
  inline double GetCompressibilityFactor(const int phase = 0) { 
    if (phase == this->iVAPPH_TP) {
      return this->Z_V;
    } else if (phase == this->iLIQPH_TP) {
      return this->Z_L;
    } else {
      return this->P * this->MW_M / this->rho / this->gas_constant_universal / this->T; 
    }
  }

  inline double GetdZdT(const int phase = 0) {
    if (phase == this->iVAPPH_TP) {
      return this->dZdT_V;
    } else if (phase == this->iLIQPH_TP) {
      return this->dZdT_L;
    } else {
      return this->dZdT;
    }
  }

  inline void GetdZdChi(double* dZdX, const int phase = 0) {
    if (phase == this->iVAPPH_TP) {
      LOOP_k_N(this->num_species) dZdX[k] = this->dZVdy[k];
    } else if (phase == this->iLIQPH_TP) {
      LOOP_k_N(this->num_species) dZdX[k] = this->dZLdx[k];
    } else {
      LOOP_k_N(this->num_species) dZdX[k] = this->dZdXk[k];
    }
  }

  inline void GetdAmdX(double* dAmdX, double* d2AmdXT = nullptr) {
    LOOP_k_N(this->num_species) dAmdX[k] = this->dAmdXk[k];
    if (d2AmdXT) {
      LOOP_k_N(this->num_species) d2AmdXT[k] = this->d2AmdTdXk[k];
    }
  }

  inline double GetDBetaDT() { return this->dBetadT; }
  inline double GetDRhoDT(const int phase = 0) {
    if (phase == this->iVAPPH_TP) {
      return this->drhodT_V;
    } else if (phase == this->iLIQPH_TP) {
      return this->drhodT_L;
    } else {
      return this->drhodT;
    }
  }
  inline double GetDRhoDp(const int phase = 0) {
    if (phase == this->iVAPPH_TP) {
      return this->drhodp_V;
    } else if (phase == this->iLIQPH_TP) {
      return this->drhodp_L;
    } else {
      return this->drhodp;
    }
  }
  inline double GetMWPhase(const int phase = 0) {
    if (phase == this->iVAPPH_TP) {
      return this->MW_V;
    } else if (phase == this->iLIQPH_TP) {
      return this->MW_L;
    } else {
      return this->MW_M;
    }
  }

  inline double GetdMWdT(const int phase = 0) {
    if (phase == this->iVAPPH_TP) {
      return this->dMWVdT;
    } else if (phase == this->iLIQPH_TP) {
      return this->dMWLdT;
    } else {
      return 0.0; // No derivative for mixture
    }
  }

  inline void GetC2C3(double& c2_out, double& c3_out, double* dc2dXk_out = nullptr, double* dc3dXk_out = nullptr) {
    c2_out = this->c2;
    c3_out = this->c3;
    if (dc2dXk_out) {
      LOOP_k_N(this->num_species) dc2dXk_out[k] = this->dc2dXk[k];
    }
    if (dc3dXk_out) {
      LOOP_k_N(this->num_species) dc3dXk_out[k] = this->dc3dXk[k];
    }
  }

 private:
  // flags, they are constant once initialized
  int iTWOPH_TP;
  int iLIQPH_TP;
  int iVAPPH_TP;
  int iMINGIBBSPH_TP;
  int iSINGLEPH_TP;
  int iSOLIDPH_TP;
  int iFAKEPH_TP;

  int phase; // phase of cell based on thermodynamic conditions.

  // Two-phase properties
  double* x_TP; // liquid mole fractions
  double* y_TP; // vapor mole fractions
  double betaV; // vapor quality

  // I'm so desperate...
  double c2;
  double c3;
  double* dc2dXk;
  double* dc3dXk;

  // molecular weights
  double MW_L; // liquid molecular weight
  double MW_V; // vapor molecular weight
  double dMWLdT; // derivative of liquid molecular weight with respect to T
  double dMWLdp; // derivative of liquid molecular weight with respect to p
  double dMWVdT; // derivative of vapor molecular weight with respect to T
  double dMWVdp; // derivative of vapor molecular weight with respect to p
  // compressibility factors
  double Z_L; // liquid compressibility factor
  double Z_V; // vapor compressibility factor
  // densities
  double rho_L; // liquid density
  double rho_V; // vapor density
  double rhoYF_L; // liquid density
  double rhoYF_V; // vapor density
  double Z_phi;
  double phi_rho; // for rho reg
  // energies
  double E_L;
  double E_V;
  double depE_L, depE_V;
  // fugacity coefficients
  double* phi_L; // liquid fugacity coefficient
  double* phi_V; // vapor fugacity coefficient
  // fugacities
  double* f_L; // liquid fugacity
  double* f_V; // vapor fugacity

  // fugacities
  double* f_fugacity;
  double* phi_fugacity;

  // derivatives
  double* d_logPhiV_dy; // num_species * num_species;
  double* d_logPhiL_dx; // num_species * num_species;
  double* d_logPhiV_dT; // num_species
  double* d_logPhiL_dT; // num_species
  double* d_logPhiV_dp; // num_species
  double* d_logPhiL_dp; // num_species

  double *alpha; // derivative of alpha with respect to T
  double *dalphadT; // derivative of alpha with respect to T
  double *d2alphadT2; // derivative of alpha with respect to T

  double *dSqAlphaAlphadT_IJ;
  double *d2SqAlphaAlphadT2_IJ;
  double* dBmdXk;
  double* dAmdXk;
  double* d2AmdTdXk;
  double* dZdXk;

  double dZdT, dZdp;
  double dBmdT, dBmdp;
  double dAmdp;

  double dZdT_L;
  double dZdT_V;
  double dZdp_L;
  double dZdp_V;
  double* dZLdx;
  double* dZVdy;

  double dBetadT;
  double dBetadp;

  double dpdT;
  double dpdT_L;
  double dpdT_V;

  double dpdrho;
  double dpdrho_L;
  double dpdrho_V;

  double* drhodchi;
  double* drhoVdy;
  double* drhoLdx;

  double dDepEdT;
  double dDepEdp;
  double dDepEdrho;
  double* dDepEdXk;

  double dDepEdT_L;
  double dDepEdp_L;
  double dDepEdrho_L;
  double* dDepEdXk_L;

  double dDepEdT_V;
  double dDepEdp_V;
  double dDepEdrho_V;
  double* dDepEdXk_V;

  double dEdT;
  double dEdT_V;
  double dEdT_L;

  double dEdp;
  double dEdp_V;
  double dEdp_L;
  double* dEdchi_V;
  double* dEdchi_L;

  double drhodT;
  double drhodT_V;
  double drhodT_L;
  double drhodp;
  double drhodp_V;
  double drhodp_L;

  double *dydT;
  double *dxdT;
  double *dydp;
  double *dxdp;

  double *dKdp;
  double *dKdT;

  double v;
  double vg;
  double vl;


  double alpha_p;
  double kappa_T;
  double kappa_S;

  double e;

  double Y_clip_low;
  double Y_clip_high;

  bool trivial_vle;
  bool m_includeChem = false;
  int niter_ssi;
  int flag_vle;
 public:
  inline int GetFlagVLE() const { return flag_vle; }
  inline int GetNiterSSI() const { return niter_ssi; }
};



#endif // PENGROBINSON_H
