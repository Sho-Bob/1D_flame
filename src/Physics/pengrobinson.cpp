#include "pengrobinson.h"

#include "Math/poly34.h"

// #include "commonTestVariables.h"
#include "pugixml.hpp"
#include <iostream>

#include "boost/math/tools/roots.hpp"
#include <Eigen/Dense>

namespace Physics {

#define SAFEGUARD(val, name) if (std::isnan(val) || std::isinf(val)) { std::cerr << "Critical error: " << name << " is not valid" << std::endl; throw(-1); }
#define SAFEGUARDINFO(val, name) if (std::isnan(val) || std::isinf(val)) { printf("T, P, rho, Y[0] = %g, %g, %g, %g\n", this->T, this->P, this->rho, this->Y[0]); std::cerr << "Critical error: " << name << " is not valid" << std::endl; throw(-1); }
#define REPLACEIFWEIRD(val) if (std::isnan(val) || std::isinf(val)) { val = 0.0; }

PengRobinson::PengRobinson() = default;

void PengRobinson::init(const std::vector<std::string>& speciesNames, const double Pamb, std::string mixing_rule) {
  this->gas_constant_universal = 8314.4621; // [J kmol^-1 K^-1]
  this->delta1 = 1. + sqrt(2.0); this->delta2 = 1. - sqrt(2.0); // PR constants
  this->delta1_p_delta2 = 2.; this->delta1delta2 = -1.0;
  this->delta1_m_delta2 = 2. * sqrt(2.0); // PR constants
  this->N_Av = 6.022e26; // kmol^-1
  this->Boltzman = 1.380649e-23;  // J/K
  this->mixing_rule = mixing_rule;

  this->n_species = speciesNames.size();
  this->num_species = this->n_species; // just to make things easier
  printf("Initializing PengRobinson with %d species\n", this->n_species);
  this->Pamb = Pamb;

  // phase flags
  this->iTWOPH_TP = 0;
  this->iLIQPH_TP = 1;
  this->iVAPPH_TP = 2;

  // array initialization
  this->species.resize(this->n_species);

  this->X.resize(this->n_species);
  this->Y.resize(this->n_species);
  this->MW.resize(this->n_species);
  this->molar_density.resize(this->n_species);

  this->Tcrit.resize(this->n_species);
  this->Pcrit.resize(this->n_species);
  this->rhocrit.resize(this->n_species);
  this->Vcrit.resize(this->n_species);
  this->Zcrit.resize(this->n_species);
  this->omega.resize(this->n_species);
  this->sigma.resize(this->n_species);
  this->dipole.resize(this->n_species);
  this->epsOverk.resize(this->n_species);

  this->NasaCoef = 7;
  this->nasa_poly_coeff.resize(this->NasaCoef * this->n_species * 2);
  this->nasa_poly_bounds.resize(this->n_species * 2);
  this->hspecies_ig.resize(this->n_species);
  this->cpspecies_ig.resize(this->n_species);
  this->sspecies_ig.resize(this->n_species);

  this->Tcrit_IJ.resize(this->n_species * this->n_species);
  this->Vcrit_IJ.resize(this->n_species * this->n_species);
  this->Zcrit_IJ.resize(this->n_species * this->n_species);
  this->Pcrit_IJ.resize(this->n_species * this->n_species);
  this->omega_IJ.resize(this->n_species * this->n_species);

  this->cst_b.resize(this->n_species);
  this->cst_a.resize(this->n_species * this->n_species);
  this->cst_c.resize(this->n_species * this->n_species);

  this->A_IJ.resize(this->n_species * this->n_species);
  this->dAmdN.resize(this->n_species);
  this->d2AmdTdN.resize(this->n_species);
  this->dK1dN.resize(this->n_species);
  this->dPdN.resize(this->n_species);
  this->dVdN.resize(this->n_species);

  this->hspecies.resize(this->n_species);
  this->sspecies.resize(this->n_species);
  this->muspecies.resize(this->n_species);
  this->muspecies_molar.resize(this->n_species);
  this->edspecies.resize(this->n_species);
  this->alphaD_IJ.resize(this->n_species * this->n_species);

  // Transport models
  this->sigma_IJ.resize(this->n_species * this->n_species);
  this->epsOverk_IJ.resize(this->n_species * this->n_species);
  this->omega_visc_IJ.resize(this->n_species * this->n_species);
  this->MW_IJ.resize(this->n_species * this->n_species);

  this->Dm.resize(this->n_species * this->n_species);
  this->alphaT.resize(this->n_species * this->n_species);

  this->m_IJ.resize(this->n_species * this->n_species);
  this->sigma_diff_IJ.resize(this->n_species * this->n_species);

  this->k0.resize(this->n_species);
  this->k1.resize(this->n_species);
  this->k2.resize(this->n_species);
  this->Akappa.resize(this->n_species);
  this->Bkappa.resize(this->n_species);
  this->kappa.resize(this->n_species * this->n_species);
  this->kappa_spec.resize(this->n_species);
  this->c_spec.resize(this->n_species);

  // Onsager
  this->Lij.resize((this->n_species-1)*(this->n_species-1));
  this->Lqi.resize(this->n_species-1);
  this->kT.resize(this->n_species-1);

  // derivatives and TP
  this->x_TP = new double [this->num_species];
  this->y_TP = new double [this->num_species];

  this->dc2dXk = new double[this->num_species];
  this->dc3dXk = new double[this->num_species];

  this->f_L = new double[num_species];
  this->f_V = new double[num_species];
  this->phi_L = new double[num_species];
  this->phi_V = new double[num_species];

  this->alpha = new double[num_species];
  this->dalphadT = new double[num_species];
  this->d2alphadT2 = new double[num_species];

  this->dSqAlphaAlphadT_IJ = new double[num_species * num_species];
  this->d2SqAlphaAlphadT2_IJ = new double[num_species * num_species];
  this->phi_fugacity = new double[num_species];
  this->f_fugacity = new double[num_species];
  this->d_logPhiV_dy = new double[num_species * num_species];
  this->d_logPhiL_dx = new double[num_species * num_species];
  this->d_logPhiV_dT = new double[num_species];
  this->d_logPhiL_dT = new double[num_species];
  this->d_logPhiV_dp = new double[num_species];
  this->d_logPhiL_dp = new double[num_species];
  this->dAmdXk = new double[num_species];
  this->d2AmdTdXk = new double[num_species];
  this->dBmdXk = new double[num_species];
  this->dZdXk = new double[num_species];
  this->dZVdy = new double[num_species];
  this->dZLdx = new double[num_species];

  this->dydT = new double[num_species];
  this->dxdT = new double[num_species];
  this->dydp = new double[num_species];
  this->dxdp = new double[num_species];

  this->dKdp = new double[num_species];
  this->dKdT = new double[num_species];

  this->drhodchi = new double[num_species];
  this->drhoVdy = new double[num_species];
  this->drhoLdx = new double[num_species];

  this->dDepEdXk   = new double[num_species];
  this->dDepEdXk_L = new double[num_species];
  this->dDepEdXk_V = new double[num_species];
  this->dEdchi_L = new double[num_species];
  this->dEdchi_V = new double[num_species];

  // init constants
  LOOP_l_N(this->n_species)
    this->species[l] = speciesNames[l];

  this->ReadNasaPolynomials();
  this->ReadCriticalProperties();
  this->SetRealFluidConstants();
  this->SetTransportConstants();
}

//----------------------------------------------------------------------------

PengRobinson::~PengRobinson() = default;

//----------------------------------------------------------------------------

void PengRobinson::ReadNasaPolynomials() {
  // std::string nasa_poly_dbpath = common::testVariables::CharlesX_Data_Root
  //     + "/Data/Physics/RealFluid/nasa_poly.xml";
  std::string nasa_poly_dbpath = "/Users/benkris/Documents/research/multiphase/nanoscale/src/NASAPOLY/nasa_poly.xml";

  std::vector<double> readFromXMLFile_low;
  std::vector<double> readFromXMLFile_high;

  pugi::xml_document xmlDoc;
  pugi::xml_parse_result result = xmlDoc.load_file(nasa_poly_dbpath.c_str());

  LOOP_l_N(this->n_species) {
    // search for species
    pugi::xml_node
      species = xmlDoc.child("ctml").child("speciesData").first_child();
    while (species.attribute("name").value() != this->species[l]) {
      species = species.next_sibling();
    }

    pugi::xml_node lowTemp = species.child("thermo").first_child();
    pugi::xml_node highTemp = lowTemp.next_sibling();

    // read temperature ranges
    this->nasa_poly_bounds[0 + 2 * l].first =
      lowTemp.attribute("Tmin").as_double();
    this->nasa_poly_bounds[0 + 2 * l].second =
      highTemp.attribute("Tmin").as_double();
    this->nasa_poly_bounds[1 + 2 * l].first =
      lowTemp.attribute("Tmax").as_double();
    this->nasa_poly_bounds[1 + 2 * l].second =
      highTemp.attribute("Tmax").as_double();

    readFromXMLFile_low =
      common::Tokenize(lowTemp.child("floatArray").child_value(), " ,\n");
    readFromXMLFile_high =
      common::Tokenize(highTemp.child("floatArray").child_value(), " ,\n");

    LOOP_k_N(this->NasaCoef) {
      this->nasa_poly_coeff[k + this->NasaCoef * l] = readFromXMLFile_low[k];
      this->nasa_poly_coeff[k + this->NasaCoef * l
        + this->NasaCoef * this->n_species] = readFromXMLFile_high[k];
    }
  }
}

//----------------------------------------------------------------------------

void PengRobinson::ReadCriticalProperties() {
  //! \brief Reads the XML data-base and sets the coefficients
  LOOP_l_N(this->n_species) {
    if (this->species[l] == "O2") {
      this->MW[l] = 31.9988; //O2
      this->Tcrit[l] = 154.5800;
      this->Pcrit[l] = 5.0430e+6;
      this->rhocrit[l] = 436.140;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.0222;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "H2") {
      this->MW[l] = 2.01588; //H2
      this->Tcrit[l] = 33.1450;
      this->Pcrit[l] = 1.2964e+6;
      this->rhocrit[l] = 31.262;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = -0.219;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "H2O") {
      this->MW[l] = 18.01528;
      this->Tcrit[l] = 647.096;
      this->Pcrit[l] = 22.0640e+6;
      this->rhocrit[l] = 322.0;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.3443;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "CO2") {
      this->MW[l] = 44.01;
      this->Tcrit[l] = 304.2;
      this->Pcrit[l] = 7.3765e+6;
      this->Zcrit[l] = 0.274;
      this->Vcrit[l] = this->Zcrit[l] * this->gas_constant_universal * this->Tcrit[l] / this->Pcrit[l];
      this->rhocrit[l] = this->MW[l] / this->Vcrit[l];
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.225;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "OH") {
      this->MW[l] = 17.0073; //OH kg/kmol
      this->Tcrit[l] = 0.0;
      this->Pcrit[l] = 0.0;
      this->rhocrit[l] = 0.0;
      this->Vcrit[l] = 0.0;
      this->Zcrit[l] = 0.0;
      this->omega[l] = 0.0;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "N2") {
      this->MW[l] = 28.0134; //N2 kg/kmol
      this->Tcrit[l] = 126.1900;
      this->Pcrit[l] = 3.3958e+6;
      this->rhocrit[l] = 313.3;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.03720;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "NH3") {
      this->MW[l] = 17.031; //NH3 kg/kmol
      this->Tcrit[l] = 405.40;
      this->Pcrit[l] = 11.3330e+6;
      this->rhocrit[l] = 225.000;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.25601;
      this->dipole[l] = 1.47;
    } else if (this->species[l] == "C8H18,isooctane") {
      this->MW[l] = 114.23;
      this->Tcrit[l] = 543.9;
      this->Pcrit[l] = 25.7e5;
      this->rhocrit[l] = 244.4522;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.394;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "NC12H26") {
      this->MW[l] = 170.33484;
      this->Tcrit[l] = 658.10;
      this->Pcrit[l] = 1.817e+6;
      this->rhocrit[l] = 226.5453372;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.574;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "CH4") {
      this->MW[l] = 16.04;
      this->Tcrit[l] = 190.6;
      this->Pcrit[l] = 46.1e5;
      this->rhocrit[l] = 162.0;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.011;
      this->dipole[l] = 0.0;
    } else if (this->species[l] == "C7H16,n-heptane") {
      this->MW[l] = 100.2019;
      this->Tcrit[l] = 540;
      this->Pcrit[l] = 27.4e5;
      this->rhocrit[l] = 2.35*100.2019;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
        (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.349;
      this->dipole[l] = 0.0;
    } else {
      std::cout << " WARNING -> Unknown species :[" << this->species[l]
        << "]. No critical properties found." << std::endl;
    }
  }
  // Corellation for sigma and epsOverk
  LOOP_l_N(this->n_species) {
    this->sigma[l] = 0.809 * std::pow(1000.0 * this->Vcrit[l], 1.0 / 3.0);
    this->epsOverk[l] = (this->Tcrit[l] / 1.2593);
  }

  return;
}

//----------------------------------------------------------------------------

// void PengRobinson::SetRealFluidConstants() {
//   LOOP_k_N(this->n_species) {
//     LOOP_l_N(this->n_species) {
//       double tmp_k = 0; // zero binary-interaction parameter
//                         // if (k != l)
//                         //   tmp_k = 0.156;
//
//       int apos = k * this->n_species + l;
//
//       this->Tcrit_IJ[apos] = std::sqrt(this->Tcrit[l] * this->Tcrit[k]) * (1.0 - tmp_k);
//       this->Vcrit_IJ[apos] = std::pow(std::pow(this->Vcrit[l], 1.0 / 3.0) + std::pow(this->Vcrit[k], 1.0 / 3.0), 3.0) / 8.0;
//       this->Zcrit_IJ[apos] = 0.5 * (this->Zcrit[l] + this->Zcrit[k]);
//       this->Pcrit_IJ[apos] = this->Zcrit_IJ[apos] * this->gas_constant_universal * this->Tcrit_IJ[apos] / this->Vcrit_IJ[apos];
//       this->omega_IJ[apos] = 0.5 * (this->omega[l] + this->omega[k]);
//     }
//   }
//
//   LOOP_k_N(this->n_species) {
//     this->cst_b[k] = 0.077796 * this->gas_constant_universal * this->Tcrit[k] / this->Pcrit[k];
//
//     LOOP_l_N(this->n_species) {
//       int apos = k * n_species + l;
//       this->cst_a[apos] = 0.457236 * std::pow(this->gas_constant_universal * this->Tcrit_IJ[apos], 2.0)  / this->Pcrit_IJ[apos];
//       if (this->omega_IJ[apos] <=0.49) // this is m_i
//         this->cst_c[apos] = 0.37464 + 1.54226 * this->omega_IJ[apos] - 0.26992 * std::pow(this->omega_IJ[apos], 2);
//       else
//         this->cst_c[apos] = 0.379642 + 1.485030*this->omega_IJ[apos] - 0.164423*std::pow(this->omega_IJ[apos], 2) + 0.016666*std::pow(this->omega_IJ[apos], 3);
//     }
//   }
// }

//----------------------------------------------------------------------------

void PengRobinson::SetTransportConstants() {
  // Set constants for Chung viscosity/conductivity model and Chapman-Enskog model
  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * n_species + l;
      this->sigma_IJ[apos] = std::sqrt(this->sigma[l] * this->sigma[k]);
      this->epsOverk_IJ[apos] = std::sqrt(this->epsOverk[l] * this->epsOverk[k]);
      this->omega_visc_IJ[apos] =
        0.5 * (std::max(this->omega[l], 0.0) + std::max(this->omega[k], 0.0));
      this->MW_IJ[apos] =
        2.0 * this->MW[l] * this->MW[k] / (this->MW[l] + this->MW[k]);

      double mi = this->MW[k] / this->N_Av;
      double mj = this->MW[l] / this->N_Av;
      this->m_IJ[apos] = mi * mj / (mi + mj);
      this->sigma_diff_IJ[apos] = 0.5 * (this->sigma[l] + this->sigma[k])* 1e-10;

      // thermal diffusion can be set here as well
      this->alphaT[apos] = 0.0567 * (mi - mj) / (mi + mj);
    }
    // For Capillary coefficient
    this->k0[k] = -2.985 + 4.332 * this->Zcrit[k] + 10.859 * this->Zcrit[k] * this->Zcrit[k] - 1.990 * this->omega[k] + 1.798 * this->omega[k] * this->omega[k];
    this->k1[k] = -0.965 + 1.405 * this->Zcrit[k] + 2.764 * this->Zcrit[k] * this->Zcrit[k] - 0.963 * this->omega[k] + 1.346 * this->omega[k] * this->omega[k];
    this->k2[k] = -0.121 + 0.156 * this->Zcrit[k] + 0.540 * this->Zcrit[k] * this->Zcrit[k] - 0.123 * this->omega[k] + 0.090 * this->omega[k] * this->omega[k];
    LOOP_k_N(this->n_species) {
      if (this->species[k] == "N2") {
        this->Akappa[k] = -0.8488e-16;// * std::pow(this->N_Av, 2./3.);
        this->Bkappa[k] = 1.2363e-16;// * std::pow(this->N_Av, 2./3.);
      } else if (this->species[k] == "NC12H26") {
        this->Akappa[k] = -0.4698e-16;// * std::pow(this->N_Av, 2./3.);
        this->Bkappa[k] = 0.5531e-16;// * std::pow(this->N_Av, 2./3.);
      } else if (this->species[k] == "C7H16,n-heptane") {
        this->Akappa[k] = -0.5874e-16;// * std::pow(this->N_Av, 2./3.);
        this->Bkappa[k] = 0.7019e-16;// * std::pow(this->N_Av, 2./3.);
      } else if (this->species[k] == "H2O") {
        this->Akappa[k] = 1e-16; // not gonna be used anyways...
        this->Bkappa[k] = 1e-16; // not gonna be used anyways...
      } else if (this->species[k] == "CH4") {
        this->Akappa[k] = 1e-16; // not gonna be used anyways...
        this->Bkappa[k] = 1e-16; // not gonna be used anyways...
      } else if (this->species[k] == "CO2") {
        this->Akappa[k] = 1e-16; // not gonna be used anyways...
        this->Bkappa[k] = 1e-16; // not gonna be used anyways...
      } else {
        std::cout << "Species " << this->species[k] << " not supported.";
        throw -1;
      }
    }
  }
}

//----------------------------------------------------------------------------

double PengRobinson::GetRho_SetMixture_TPY(const double T, const double P, const double *Y) {
  this->T = T;
  this->P = P;
  this->SetMassFractionFromY(Y);

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncRZFromPressureTemperature();

  // molar density
  this->SyncMolarDensities();

  return this->rho;
}

double PengRobinson::GetRho_SetMixture_TPX(const double T, const double P, const double *X) {
  this->T = T;
  this->P = P;
  this->SetMolarFractionFromX(X);

  this->SetMolecularWeightMixtureFromX();
  this->SetMassFractionFromX();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncRZFromPressureTemperature();
  this->SyncTransportFromTemperatureDensity();
  // molar density
  this->SyncMolarDensities();

  this->rho = this->GetDensityFromPressureTemperature_TP(P, T, 0);

  return this->rho;
}

double PengRobinson::GetRho_SetMixture_TPRY(const double T, const double P, const double RYF, double rhoGuess) {
  assert(this->n_species==2);
  double R, nextR, error, nextError;
  R = rhoGuess;
  nextR = rhoGuess*1.00001;
  do {
    error = (P - this->helperFun(T, R, RYF)) / P;
    nextError = (P - this->helperFun(T, nextR, RYF)) / P;
    double temp = R;
    R = nextR;
    nextR = temp + (nextR - temp) * (0. - error) / (nextError - error);
  } while (std::fabs(error) > 1e-12);
  // while (std::fabs((R-nextR)/R) > 1e-10);
  return R;
}

double PengRobinson::helperFun(double T, double R, double RYF) {
  assert(this->n_species==2);
  double YY[2] = {RYF/R, 1.-RYF/R};
  this->T = T;
  this->rho = R;
  this->SetMassFractionFromY(YY);

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncPZFromTemperatureDensity();
  return this->GetP();
}

//----------------------------------------------------------------------------

/*
 * @brief: Given T and P, compute the equilibrium XL and XV using a flash solver using a VLE model
 */
void PengRobinson::TPFlash(const double T, const double P, const double *XvGuess, const double *XlGuess, double *Xv, double *Xl, double& Tout) {

  // instead, use VLE
  double beta_out;
  this->VLE(T, P, XvGuess, beta_out, Xv, Xl);
  printf("VLE done: beta_out = %.4f\n", beta_out);
  LOOP_k_N (this->num_species) {
    printf("xv[%d] = %.4f, ", k, Xv[k]);
    printf("xl[%d] = %.4f\n", k, Xl[k]);
  }
  Tout = T;
  return; 
  // printf("T = %.4f K, P = %.4f MPa\n", T, P/1e6);
  double Tguess = T;
  if (this->n_species == 2) { // only if 2-species...
    double XvFGuess = XvGuess[0];
    printf("XvGuess = %.4f\n", XvFGuess);
    // Lambda of the objective fun
    auto objFun = [&](double x) {
      double Xvin[2] = {x, 1.-x};
      double Xlout[2];
      double ret;
      this->GetTdew_PXvapor(P, Xvin, Tguess, XlGuess, ret, Xlout);
      // printf("ret = %.4f K\n", ret);
      return ret - T;
    };

    double step = 1e-3;
    bool isRising = (objFun(XvFGuess*(1+step))-objFun(XvFGuess)) > 0;

    // Max iterations
    const boost::uintmax_t maxit = 500;
    boost::uintmax_t it = maxit;

    double factor = 1.1;
    boost::math::tools::eps_tolerance<double> tol(36);

    std::pair<double, double> sol = boost::math::tools::bracket_and_solve_root(objFun, XvFGuess, factor, isRising, tol, it);

    if (it >= maxit)
      throw std::runtime_error("PengRobinson::TPFlash(): max iteration achieved.");

    double XvF = 0.5 * (sol.first + sol.second);
    Xv[0] = XvF;
    Xv[1] = 1.-XvF;

    this->GetTdew_PXvapor(P, Xv, Tguess, XlGuess, Tout, Xl);
  } else if (this->n_species == 1) { // for single-species, the result if obviously unity molar fractions
    Xv[0] = 1.0;
    Xl[0]= 1.0;
    Tout = Tguess;
  } else {
    std::cout << this->n_species << " species simulations for TPFlash are not implemented yet. Exiting!" << std::endl;
    throw(-1);
  }
}

double PengRobinson::TPFlashForce(const double T, const double P, double *Xv, double *Xl, double& Tout) {
  // Basically TPFlash, but force to get the vapor-liquid equilibrium
  double beta_out; double XvGuess[this->num_species];
  int nphase;
  double XF = 0.05;

  do {
    XvGuess[0] = XF; XvGuess[1] = 1.0 - XF;
    nphase = this->VLE(T, P, XvGuess, beta_out, Xl, Xv);
    XF += 0.05; // increase the vapor fraction guess
    if (XF > 0.95) {
      std::cout << "TPFlashForce: Could not find a VLE with 2 phases." << std::endl;
      throw -1;
    }
  } while (nphase != 2); // loop until we have a VLE with 2 phases
  printf("VLE done: beta_out = %.4f\n", beta_out);
  LOOP_k_N (this->num_species) {
    printf("xv[%d] = %.4f, ", k, Xv[k]);
    printf("xl[%d] = %.4f\n", k, Xl[k]);
  }
  Tout = T;
  return XF;
}

/*
 * brief: Calculate Dew Temperature given pressure P and vapor fraction Xvapor
 */
void PengRobinson::GetTdew_PXvapor(const double P, const double *Xvapor, const double TGuess, const double *XliqGuess, double &Tdew_out, double *Xliq_out) {
  assert(this->n_species == 2);
  // std::cout << "Starting GetTdew_PXvapor..." << std::endl;
  // First iter
  double T = TGuess;
  double Xliq [2] = {XliqGuess[0], XliqGuess[1]};

  double phil[2], phiv[2];
  this->getPhifromTandPandX(T, P, Xliq, phil);
  this->getPhifromTandPandX(T, P, Xvapor, phiv);
  // printf("phi_liquid = %.4f, %.4f\n", phil[0], phil[1]);
  // printf("phi_vapor = %.4f, %.4f\n", phiv[0], phiv[1]);
  double K[2] = {phil[0] / phiv[0], phil[1] / phiv[1]};
  double xT = Xvapor[0] / K[0] + (1.-Xvapor[0]) / K[1];
  // printf("xT = %.4f\n", xT);
  // printf("T = %.4f\n", T);

  double xT_oldold = xT;
  double T_oldold = T;

  // Second iter
  T = T * xT;
  Xliq[0] = Xvapor[0] / K[0] / xT; Xliq[1] = 1. - Xliq[0];

  this->getPhifromTandPandX(T, P, Xliq, phil);
  this->getPhifromTandPandX(T, P, Xvapor, phiv);
  K[0] = phil[0] / phiv[0]; K[1] = phil[1] / phiv[1];
  xT = Xvapor[0] / K[0] + (1.-Xvapor[0]) / K[1];

  double xT_old = xT;
  double T_old = T;

  // If we need more iter
  while (std::fabs(xT-1) > 1e-12) {
    T = T_oldold + (1-xT_oldold)/(xT_old-xT_oldold)*(T_old-T_oldold);
    Xliq[0] = Xvapor[0] / K[0] / xT; Xliq[1] = 1. - Xliq[0];

    this->getPhifromTandPandX(T, P, Xliq, phil);
    this->getPhifromTandPandX(T, P, Xvapor, phiv);
    K[0] = phil[0] / phiv[0]; K[1] = phil[1] / phiv[1];
    xT = Xvapor[0] / K[0] + (1.-Xvapor[0]) / K[1];

    xT_oldold = xT_old;
    T_oldold = T_old;
    xT_old = xT;
    T_old = T;
  }
  Tdew_out = T;
  Xliq_out[0] = Xliq[0]; Xliq_out[1] = Xliq[1];
}


void PengRobinson::getPhifromTandPandX(const double T_in, const double P_in, const double *X_in, double *phiOut) {
  double B = 0.;
  LOOP_k_N(this->n_species)
    B += this->cst_b[k] * X_in[k];

  double A = 0.;
  double A_IJ [this->n_species * this->n_species];
  double dAdN [this->n_species];
  LOOP_k_N(this->n_species) {
    dAdN[k] = 0.;
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      double X_X = X_in[l] * X_in[k];
      A_IJ[apos] = this->cst_a[apos] * std::pow(1.0 + this->cst_c[apos] * (1.0 - std::sqrt(T_in / this->Tcrit_IJ[apos])), 2);
      A += X_X * A_IJ[apos];
      dAdN[k] += X_in[l] * A_IJ[apos];
    }
    dAdN[k] *= 2.;
  }

  double Amix = A * P_in / std::pow(this->gas_constant_universal, 2) / std::pow(T_in, 2);
  double Bmix = B * P_in / this->gas_constant_universal / T_in;

  double a0 = -(Amix * Bmix - Bmix * Bmix - Bmix * Bmix * Bmix);
  double a1 = Amix - 3 * Bmix * Bmix - 2 * Bmix;
  double a2 = -(1. - Bmix);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; i++)
    if (xZ[i] > Bmix)
      Z.push_back(xZ[i]);
  double Zout;
  if (Z.size() == 1) {
    Zout = Z[0];
  } else {
    std::vector<double> lnPhi(Z.size());

    int iMin = 0;
    for (int i = 0; i < Z.size(); i++) {
      lnPhi[i] = -std::log(Z[i] - Bmix)
        - Amix / Bmix / std::sqrt(8) * std::log((Z[i] + (1 + std::sqrt(2)) * Bmix) / (Z[i] + (1 - std::sqrt(2)) * Bmix))
        + Z[i] - 1;
      if (lnPhi[i] < lnPhi[iMin]) iMin = i;
    }
    Zout = Z[iMin];
  }

  LOOP_k_N(this->n_species) {
    double Bi = this->cst_b[k] * P_in / this->gas_constant_universal / T_in;
    phiOut[k] = std::exp(Bi / Bmix * (Zout - 1.) - std::log(Zout - Bmix)
        - Amix / Bmix / 2.8284 * std::log((Zout+2.4142*Bmix)/(Zout-.4142*Bmix)) * (dAdN[k] / A - Bi/Bmix));
  }

}

//----------------------------------------------------------------------------

void PengRobinson::SetMixture_TRY(const double T, const double R, const double *Y) {
  this->T = T;
  this->rho = R;
  this->SetMassFractionFromY(Y);

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncPZFromTemperatureDensity();

  this->SyncIdealFluidThermodynamicsFromTemperature();
  this->SyncEFromTemperatureDensity();

  return; 
  this->SyncPartialPropertiesFromTemperatureDensity();
  this->SyncTransportFromTemperatureDensity();
  // molar density
  this->SyncMolarDensities();

}

void PengRobinson::SetMixture_TRX(const double T, const double R, const double *X) {
  this->T = T;
  this->rho = R;
  this->SetMolarFractionFromX(X);

  this->SetMolecularWeightMixtureFromX();
  this->SetMassFractionFromX();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncPZFromTemperatureDensity();

  this->SyncIdealFluidThermodynamicsFromTemperature();
  this->SyncEFromTemperatureDensity();

  this->SyncPartialPropertiesFromTemperatureDensity();
  this->SyncTransportFromTemperatureDensity();
  // molar density
  this->SyncMolarDensities();

}

void PengRobinson::SetMixture_EKorRY(const double E, const double R, const double *Y, const double *gradRhogradRho, const double T_guess) {
  this->rho = R;
  this->SetMassFractionFromY(Y);

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->PrepareRealFluidThermodynamicsForTFromE();
  this->T = GetTemperatureFromEnergyKor(E, gradRhogradRho, T_guess);

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncPZFromTemperatureDensity();

  this->SyncIdealFluidThermodynamicsFromTemperature();
  this->SyncEFromTemperatureDensity();

  this->SyncPartialPropertiesFromTemperatureDensity();
  this->SyncTransportFromTemperatureDensity();
}

void PengRobinson::SetMixture_HPY(const double H, const double P_in, const double *Y, const double T_guess) {
  double convCrit = 1e-8;

  this->SetMassFractionFromY(Y);

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  auto objectiveFun = [this, H, P_in, Y](double T_in) {
    double r = this->GetRho_SetMixture_TPY(T_in, P_in, Y);
    this->SetMixture_TRY(T_in, r, Y);

    double Hout = this->GetInternalEnergy() + this->GetP()/this->GetRho();

    return (Hout - H) / std::fabs(H);
  };

  double step = 1e-3;
  bool isRising = true;

  // Max iterations
  const boost::uintmax_t maxit = 100;
  boost::uintmax_t it = maxit;

  double factor = 1.01;
  boost::math::tools::eps_tolerance<double> tol(24);

  std::pair<double, double> sol = boost::math::tools::bracket_and_solve_root(objectiveFun, T_guess, factor, isRising, tol, it);

  if (it >= maxit) {
    std::cout << "---------------------------------" << '\n'
      << "Over 100 iters: " << H << " " << T_guess << std::endl;
  }

  double T_1 = (sol.first + sol.second)/2.;

  double r = this->GetRho_SetMixture_TPY(T_1, P_in, Y);
  this->SetMixture_TRY(T_1,r,Y);
}

//----------------------------------------------------------------------------

void PengRobinson::PrepareRealFluidThermodynamicsForTFromE() {
  double v = this->MW_M / this->rho;

  this->Bm = 0.;

  LOOP_k_N(this->n_species) {
    this->Bm += this->X[k] * this->cst_b[k];
  }
  this->K1 = 1.0 / (std::sqrt(8.0) * this->Bm)
    * std::log((v + (1 - std::sqrt(2.0)) * this->Bm) / (v + (1 + std::sqrt(2.0)) * this->Bm));
}

double PengRobinson::GetTemperatureFromEnergyKor(const double E_in, const double *gradRhogradRho, double T_guess) {
  double T_0 = 300.;
  if (T_guess > 50. && T_guess < 4000.)
    T_0 = T_guess;
  double E_0 = this->EFromTFunc(T_0, gradRhogradRho);

  double T_1 = (E_0 > E_in) ? 0.999 * T_0 : 1.001 * T_0;
  double E_1 = this->EFromTFunc(T_1, gradRhogradRho);

  int counter = 0;
  //  while (std::fabs(E_1 - E_in) > 1e-8) {
  while (std::fabs(T_1 - T_0) > 1e-10) {
    double T_tmp = T_1;
    double E_tmp = E_1;

    T_1 -= (E_1 - E_in) * (T_1 - T_0) / (E_1 - E_0);
    T_1 = std::max(T_1, 20.0);
    E_1 = this->EFromTFunc(T_1, gradRhogradRho);

    T_0 = T_tmp;
    E_0 = E_tmp;
    counter++;
    if (counter > 1000) {
      std::cout << "Over 1000 iter" << std::endl;
      return T_guess;
    }
  }

  return T_1;
}

/*
 * @brief: Given T and \nabla \rho Y_i \nabla \rho Y_j, compute the energy E_{GD} (including gradient terms!)
 */
double PengRobinson::EFromTFunc(double T_in, const double *gradRhogradRho) {
  // Update Pengrobinson Factors
  this->Am = 0.;
  this->dAmdT = 0.;
  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      double X_X = this->X[l] * this->X[k];
      this->A_IJ[apos] = this->cst_a[apos] * std::pow(1.0 + this->cst_c[apos] * (1.0 - std::sqrt(T_in / this->Tcrit_IJ[apos])), 2);
      double G = this->cst_c[apos] * std::sqrt(T_in / this->Tcrit_IJ[apos])
        / (1.0 + cst_c[apos] * (1.0 - std::sqrt(T_in / this->Tcrit_IJ[apos])));
      this->Am += X_X * this->A_IJ[apos];
      this->dAmdT -= X_X * this->A_IJ[apos] * G;
    }
  }
  this->dAmdT /= T_in;

  // Find Korteweg tensor and Korteweg energy contribution
  double kappa_species[this->n_species];
  LOOP_k_N(this->n_species) {
    double t = 1. - T_in / this->Tcrit[k];
    //    if (t<0.01)
    //      t=0.01;
    if (t <= 0.) {
      kappa_species[k] = 0.;
    } else {
      double temp = std::log(t);
      double nonDimKappa = std::exp(this->k0[k] + this->k1[k] * temp + this->k2[k] * temp * temp);
      //      if (nonDimKappa > 1.) {
      //        std::cout << "Error in nonDimKappa: " << nonDimKappa << std::endl;
      //        throw -1;
      //      }
      int apos = k * this->n_species + k;
      kappa_species[k] = nonDimKappa * this->A_IJ[apos] * std::pow(this->cst_b[k] / this->N_Av, 2./3.);
    }
    double nonDimKappa = (this->Akappa[k] * t + this->Bkappa[k])/100.*0.5*(1.+std::tanh((t+1.)/0.01));
    int apos = k * this->n_species + k;
    kappa_species[k] = nonDimKappa * this->A_IJ[apos] * std::pow(this->cst_b[k], 2./3.);
  }
  kappa_species[1] = 0.;
  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      this->kappa[apos] = std::sqrt(kappa_species[k] * kappa_species[l]);
    }
  }

  double EKorOnly = 0.;
  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      EKorOnly += 0.5 * this->kappa[apos] * gradRhogradRho[apos] / this->MW[k] / this->MW[l];
    }
  }
  EKorOnly /= this->rho;

  // Find Ideal and real energy contribution
  double T1 = T_in;
  double T2 = T_in * T_in;
  double T3 = T2 * T_in;
  double T4 = T3 * T_in;

  this->h_ig = 0.;
  LOOP_l_N(this->n_species) {
    if (T_in < 1000.0) {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1 / 2.0;
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 3.0;
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 4.0;
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 5.0;
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 5] / T1;
    } else {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
        + this->NasaCoef * this->n_species];
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
        + this->NasaCoef * this->n_species] * T1 / 2.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
        + this->NasaCoef * this->n_species] * T2 / 3.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
        + this->NasaCoef * this->n_species] * T3 / 4.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
        + this->NasaCoef * this->n_species] * T4 / 5.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 5
        + this->NasaCoef * this->n_species] / T1;
    }
    this->h_ig += this->X[l] * (hspecies_ig[l]);
  }
  this->h_ig *= (T_in * this->gas_constant_universal) / this->MW_M;

  double dep = (this->Am - T_in * this->dAmdT) * this->K1 / this->MW_M;
  this->E = this->h_ig - this->gas_constant * T_in + dep;

  return this->E + EKorOnly;
}

//----------------------------------------------------------------------------

void PengRobinson::SyncRealFluidThermodynamicsFromTemperatureDensity() {
  double v = this->MW_M / this->rho;

  this->Am = 0.;
  this->Bm = 0.;
  this->dAmdT = 0.;
  this->d2AmdT2 = 0.;

  LOOP_k_N(this->n_species) {
    this->dAmdN[k] = 0.;
    this->d2AmdTdN[k] = 0.;

    this->Bm += this->X[k] * this->cst_b[k];

    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      double X_X = this->X[l] * this->X[k];
      this->A_IJ[apos] = this->cst_a[apos] * std::pow(1.0 + this->cst_c[apos] * (1.0 - std::sqrt(this->T / this->Tcrit_IJ[apos])), 2);
      double G = this->cst_c[apos] * std::sqrt(this->T / this->Tcrit_IJ[apos])
        / (1.0 + cst_c[apos] * (1.0 - std::sqrt(this->T / this->Tcrit_IJ[apos])));
      double D = this->cst_c[apos] * (1.0 + cst_c[apos]) * this->Tcrit_IJ[apos]
        / this->Pcrit_IJ[apos] * std::sqrt(this->Tcrit_IJ[apos] / this->T);
      this->Am += X_X * this->A_IJ[apos];
      this->dAmdT -= X_X * this->A_IJ[apos] * G;
      this->d2AmdT2 += X_X * D;

      this->dAmdN[k] += this->X[l] * this->A_IJ[apos];
      this->d2AmdTdN[k] += this->X[l] * this->A_IJ[apos] * G;
    }
    this->dAmdN[k] *= 2.0;
    this->d2AmdTdN[k] *= -2.0 / this->T;
  }
  this->dAmdT /= this->T;
  this->d2AmdT2 *= 0.457236 * std::pow(this->gas_constant_universal, 2) / (2.0 * this->T);
  this->dPdT = this->gas_constant_universal / (v - this->Bm)
    - this->dAmdT / (std::pow(v, 2) + 2.0 * v * this->Bm - std::pow(this->Bm, 2));
  double arg = this->gas_constant_universal * this->T * (v + this->Bm)
    * std::pow((v / (v - this->Bm) + this->Bm / (v + this->Bm)), 2);
  this->dPdV = -this->gas_constant_universal * this->T / std::pow((v - this->Bm), 2)
    * (1.0 - 2.0 * this->Am / arg);
  // Regularization of expansivity and partial volume (Jofre & Urzay 2020)
  this->BT = 1. - 2. * this->Am
    / ((v+this->Bm) * this->gas_constant_universal * this->T * std::pow(v/(v-this->Bm)+this->Bm/(v+this->Bm), 2.));
  this->BT = std::max(this->BT, 0.5);
  this->betaT = 1. / v * std::pow(v-this->Bm, 2.) / this->gas_constant_universal / this->T / this->BT;
  this->expansivity = this->betaT * this->dPdT;
  //  this->expansivity = -this->dPdT / (v * this->dPdV); //ideal gas: equal to 3.34E-3 (1/K)
  this->K1 = 1.0 / (std::sqrt(8.0) * this->Bm)
    * std::log((v + (1 - std::sqrt(2.0)) * this->Bm) / (v + (1 + std::sqrt(2.0)) * this->Bm));

  double temp = v * v + 2.0 * this->Bm * v - this->Bm * this->Bm;
  LOOP_k_N(this->n_species) {
    this->dPdN[k] = this->gas_constant_universal * this->T / (v - this->Bm) +
      this->gas_constant_universal * this->T * this->cst_b[k]
      / std::pow((v - this->Bm), 2) - this->dAmdN[k] / temp
      + 2.0 * this->Am * this->cst_b[k] * (v - this->Bm) / std::pow(temp, 2);
    // Define dVdN through betaT for regularization
    this->dVdN[k] = v * this->betaT * this->dPdN[k];
    // this->dVdN[k] = -this->dPdN[k] / this->dPdV;
    this->dK1dN[k] = 1.0 / temp * this->dVdN[k] - this->cst_b[k] / this->Bm *
      (this->K1 + v / temp);
  }
}

//----------------------------------------------------------------------------
void PengRobinson::SyncIdealFluidThermodynamicsFromTemperature() {
  double T1 = this->T;
  double T2 = this->T * this->T;
  double T3 = T2 * this->T;
  double T4 = T3 * this->T;

  this->h_ig = 0.;
  LOOP_l_N(this->n_species) {
    if (this->T < 1000.0) {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1 / 2.0;
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 3.0;
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 4.0;
      this->hspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 5.0;
      //      this->hspecies_ig[l] +=
      //          this->nasa_poly_coeff[l * this->NasaCoef + 5] / T1;
    } else {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
        + this->NasaCoef * this->n_species];
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
        + this->NasaCoef * this->n_species] * T1 / 2.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
        + this->NasaCoef * this->n_species] * T2 / 3.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
        + this->NasaCoef * this->n_species] * T3 / 4.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
        + this->NasaCoef * this->n_species] * T4 / 5.0;
      //      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 5
      //          + this->NasaCoef * this->n_species] / T1;
    }
    this->h_ig += this->X[l] * (hspecies_ig[l]);
  }
  this->h_ig *= (this->T * this->gas_constant_universal) / this->MW_M;

  this->cp_ig = 0.;
  this->cv_ig = 0.;
  LOOP_l_N(this->n_species) {
    if (this->T <= 1000.0) {
      this->cpspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->cpspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1;
      this->cpspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2;
      this->cpspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3;
      this->cpspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4;
    } else {
      this->cpspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
        + this->NasaCoef * this->n_species];
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
        + this->NasaCoef * this->n_species] * T1;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
        + this->NasaCoef * this->n_species] * T2;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
        + this->NasaCoef * this->n_species] * T3;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
        + this->NasaCoef * this->n_species] * T4;
    }
    double tmp = cpspecies_ig[l] * this->gas_constant_universal;
    this->cp_ig += this->X[l] * tmp; //cp_mol
    this->cv_ig += this->X[l] * (tmp - this->gas_constant_universal);
  }
  this->cp_ig /= this->MW_M;
  this->cv_ig /= this->MW_M;

  LOOP_l_N(this->n_species) {
    if (this->T < 1000.0) {
      this->sspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0] * std::log(T1);
      this->sspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1;
      this->sspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 2.0;
      this->sspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 3.0;
      this->sspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 4.0;
      this->sspecies_ig[l] +=
        this->nasa_poly_coeff[l * this->NasaCoef + 6];
    } else {
      this->sspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
        + this->NasaCoef * this->n_species] * std::log(T);
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
        + this->NasaCoef * this->n_species] * T1;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
        + this->NasaCoef * this->n_species] * T2 / 2.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
        + this->NasaCoef * this->n_species] * T3 / 3.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
        + this->NasaCoef * this->n_species] * T4 / 4.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 6
        + this->NasaCoef * this->n_species];
    }
    this->sspecies_ig[l] -= std::log(this->X[l]*this->PForDiff/1e5);
    if (this->X[l] < 1e-10)
      sspecies_ig[l] = 0.;
    if (std::isnan(sspecies_ig[l])) {
      std::cout << "Something wrong Ideal. T = " << this->T << ", rho = " << this->rho
        << ", YF = " << this->Y[0] << " , P = " << this->P << std::endl;
      throw -1;
    }
  }
}

//----------------------------------------------------------------------------

void PengRobinson::SyncRZFromPressureTemperature() {
  // instead, use GetDensityFromPressureTemperature_TP to set rho and Z
  // this->rho = this->GetDensityFromPressureTemperature_TP(this->P, this->T, 0);
  // this->Z = this->P / this->gas_constant_universal / this->T / this->rho;
  // if (this->rho < 0 || std::isnan(this->rho) || std::isinf(this->rho)) {
  //   std::cout << "Impossible density: " << this->rho << ", P = " << this->P
  //     << ", T = " << this->T << std::endl;
  //   printf("Y[0] = %f", this->Y[0]);
  //   std::cout << std::endl;
  //   throw -1;
  // }
  // return;
  // this->Am = this->GetAm(this->T); this->Bm = this->GetBm();
  double a = this->Am;
  double b = this->Bm;

  double A = a * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
  double B = b * this->P / this->gas_constant_universal / this->T;

  double a0 = -(A * B - B * B - B * B * B);
  double a1 = A - 3 * B * B - 2 * B;
  double a2 = -(1 - B);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; i++)
    if (xZ[i] > B)
      Z.push_back(xZ[i]);

  /*
     std::cout << "xZ = ";
     for(int idz = 0; idz < xZ.size(); idz++) {
     std::cout << std::scientific << xZ[idz];
     if (idz != xZ.size()-1)
     std::cout << ", ";
     }
     std::cout << std::endl;
     */

  if (Z.size() == 1) {
    this->Z = Z[0];
    this->rho = this->P / this->gas_constant / this->T / this->Z;
    return;
  } else {
    std::vector<double> lnPhi(Z.size());

    int iMin = 0;
    for (int i = 0; i < Z.size(); i++) {
      lnPhi[i] = -std::log(Z[i] - B)
        - A / B / std::sqrt(8) * std::log((Z[i] + (1 + std::sqrt(2)) * B) / (Z[i] + (1 - std::sqrt(2)) * B))
        + Z[i] - 1;
      if (lnPhi[i] < lnPhi[iMin]) iMin = i;
    }
    this->Z =  Z[iMin];
    this->rho = this->P / this->gas_constant / this->T / this->Z;
    return;
  }

}

//----------------------------------------------------------------------------

void PengRobinson::SyncPZFromTemperatureDensity() {
  double v = this->MW_M / this->rho; // m^3 / kmolpengro

  this->P = (this->gas_constant_universal * this->T) / (v - this->Bm)
    - this->Am / (std::pow(v, 2) + 2.0 * v * this->Bm - pow(this->Bm, 2));
  // SAFEGUARD(this->P, "P");

  this->Z = this->P * v / this->gas_constant_universal / this->T;

  double Tcmix, Pcmix, wcmix, Zcmix;
  this->ComputeMixtureCrit(Tcmix, Pcmix, wcmix, Zcmix);
  // Test if we might be in vapor dome
  this->PForDiff = this->Pamb;
  //  if (PForDiff < 0) {
  //    // Psat estimate
  //    double Psat = Pcmix * std::pow(10.0, (7. / 3. * (1. + wcmix)) * (1. - Tcmix / this->T));

  //    double a = this->Am;
  //    double b = this->Bm;

  //    double A = a * Psat / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
  //    double B = b * Psat / this->gas_constant_universal / this->T;

  //    double a0 = -(A * B - B * B - B * B * B);
  //    double a1 = A - 3 * B * B - 2 * B;
  //    double a2 = -(1 - B);

  //    std::vector<double> Z(3);
  //    int n = SolveP3(&Z[0], a2, a1, a0);
  //    if (n == 1)
  //      PForDiff = Psat; // short-cut estimate wrong with PR
  //    else {
  //      std::sort(Z.begin(), Z.end());

  //      double rhoL = Psat / this->gas_constant / this->T / Z[0];
  //      double rhoV = Psat / this->gas_constant / this->T / Z[2];

  //      if (PForDiff < 0 || (this->rho < rhoL && this->rho > rhoV))
  //        PForDiff = Psat;
  //    }
  //  }
  if (PForDiff < 0 || std::isnan(PForDiff)) {
    std::cout << "Impossible pressure for diffusion: " << PForDiff << ", " << this->P << std::endl;
    throw -1;
  }
  this->ZForDiff = this->PForDiff * v / this->gas_constant_universal / this->T;
}

//----------------------------------------------------------------------------

void PengRobinson::SyncEFromTemperatureDensity() {
  double dep = (this->Am - this->T * this->dAmdT) * this->K1 / this->MW_M;
  this->E = this->h_ig - this->gas_constant * this->T + dep;

  double departureCv = (-this->T * this->d2AmdT2 * this->K1) / this->MW_M;
  double departureCp = departureCv
    + (-this->T * std::pow(this->dPdT, 2) / this->dPdV - this->gas_constant_universal)
    / this->MW_M;
  this->cv = this->cv_ig + departureCv;
  this->cp = this->cp_ig + departureCp;
  this->gamma = this->cp / this->cv;
  double v = this->MW_M / this->rho;
  double isocompressibility = -1. / (v * this->dPdV);
  this->sos = std::sqrt(this->gamma / (this->rho * std::fabs(isocompressibility)));
  this->sos2 = this->gamma / (this->rho * (isocompressibility));
}

//----------------------------------------------------------------------------
double PengRobinson::getdmuF_muOdrhoYF() {
  assert(this->n_species == 2);
  double XF = this->X[0];
  double XO = this->X[1];
  double MF = this->MW[0];
  double MO = this->MW[1];

  double YF = this->Y[0];

  double dYFdXF = MF/(XF*MF+XO*MO) - XF*MF/(XF*MF+XO*MO)/(XF*MF+XO*MO)*(MF-MO);

  double v = this->MW_M / this->rho;
  double drhodNF = 1./v * MF - this->MW_M/v/v*this->dVdN[0];
  double drhodNO = 1./v * MO - this->MW_M/v/v*this->dVdN[1];

  double drhodXF = drhodNF - drhodNO;

  double drhodYF = drhodXF / dYFdXF;

  // Got here


  double denom = this->rho + YF*drhodYF;

  double rt = this->gas_constant_universal * this->T;
  //TODO: test these
  double dmuKdNL [2][2];
  for (int k = 0; k < 2; ++k) {
    for (int l = 0; l < 2; ++l) {
      int apos = k*2 + l;
      dmuKdNL[k][l] = 1./this->MW[k]*rt * (-this->cst_b[k] * this->cst_b[l] / this->Bm / this->Bm * (this->Z - 1.)
          + this->cst_b[k] / this->Bm * this->P / rt * this->dVdN[l]
          - 1./(this->ZForDiff - this->Bm*this->PForDiff/rt)*(this->PForDiff / rt * this->dVdN[l] - this->cst_b[l]*this->PForDiff/rt)
          - (this->dAmdN[l]/this->Bm - this->Am*this->cst_b[l]/this->Bm/this->Bm)/rt/2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm)
          - this->Am / this->Bm / rt / 2.8284 * ((this->PForDiff/rt*this->dVdN[l] + 2.4142*this->cst_b[l]*this->PForDiff/rt)/(this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt) - (this->PForDiff/rt*this->dVdN[l] - 0.4142*this->cst_b[l]*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm)
          - this->Am / this->Bm / rt / 2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (2.*this->A_IJ[apos]/this->Am - this->dAmdN[k]/this->Am/this->Am*this->dAmdN[l] + this->cst_b[k]/this->Bm/this->Bm*this->cst_b[l]));
    }
    dmuKdNL[k][k] += rt  / this->X[k] / this->MW[k];

  }

  double dmuFdXF = dmuKdNL[0][0] - dmuKdNL[0][1];
  double dmuOdXF = dmuKdNL[1][0] - dmuKdNL[1][1];

  double dmuFdYF = dmuFdXF / dYFdXF;
  double dmuOdYF = dmuOdXF / dYFdXF;

  return dmuFdYF - dmuOdYF;

  double dmuFdrhoYF = dmuFdYF / denom;
  double dmuOdrhoYF = dmuOdYF / denom;

  return dmuFdrhoYF - dmuOdrhoYF;

}

void PengRobinson::getdmuidrhoj(double** dmuidrhoj) {
  assert(this->n_species == 2);
  double XF = this->X[0];
  double XO = this->X[1];
  double MF = this->MW[0];
  double MO = this->MW[1];

  double YF = this->Y[0];
  double YO = this->Y[1];

  double dYFdXF = MF/(XF*MF+XO*MO) - XF*MF/(XF*MF+XO*MO)/(XF*MF+XO*MO)*(MF-MO);
  double dYOdXO = MO/(XF*MF+XO*MO) - XO*MO/(XF*MF+XO*MO)/(XF*MF+XO*MO)*(MO-MF);

  double v = this->MW_M / this->rho;
  double drhodNF = 1./v * MF - this->MW_M/v/v*this->dVdN[0];
  double drhodNO = 1./v * MO - this->MW_M/v/v*this->dVdN[1];

  double drhodXF = drhodNF - drhodNO;
  double drhodXO = drhodNO - drhodNF;

  double drhodYF = drhodXF / dYFdXF;
  double drhodYO = - drhodYF;
  // Got here


  // denom is drhodYFdYF
  double denom = this->rho + YF*drhodYF;
  double drhoYOdYO = this->rho + YO * drhodYO;


  double rt = this->gas_constant_universal * this->T;
  double dmuKdNL [2][2];
  for (int k = 0; k < 2; ++k) {
    for (int l = 0; l < 2; ++l) {
      int apos = k*2 + l;
      // note that I deleted the 1/MW[k] term, so this would be in J/mol
      /*
         dmuKdNL[k][l] = 1./this->MW[k]*rt * (-this->cst_b[k] * this->cst_b[l] / this->Bm / this->Bm * (this->Z - 1.)
         + this->cst_b[k] / this->Bm * this->P / rt * this->dVdN[l]
         - 1./(this->ZForDiff - this->Bm*this->PForDiff/rt)*(this->PForDiff / rt * this->dVdN[l] - this->cst_b[l]*this->PForDiff/rt)
         - (this->dAmdN[l]/this->Bm - this->Am*this->cst_b[l]/this->Bm/this->Bm)/rt/2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm)
         - this->Am / this->Bm / rt / 2.8284 * ((this->PForDiff/rt*this->dVdN[l] + 2.4142*this->cst_b[l]*this->PForDiff/rt)/(this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt) - (this->PForDiff/rt*this->dVdN[l] - 0.4142*this->cst_b[l]*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm)
         - this->Am / this->Bm / rt / 2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (2.*this->A_IJ[apos]/this->Am - this->dAmdN[k]/this->Am/this->Am*this->dAmdN[l] + this->cst_b[k]/this->Bm/this->Bm*this->cst_b[l]));
         */
      dmuKdNL[k][l] = rt * (-this->cst_b[k] * this->cst_b[l] / this->Bm / this->Bm * (this->Z - 1.)
          + this->cst_b[k] / this->Bm * this->P / rt * this->dVdN[l]
          - 1./(this->Z - this->Bm*this->P/rt)*(this->P / rt * this->dVdN[l] - this->cst_b[l]*this->P/rt)
          - (this->dAmdN[l]/this->Bm - this->Am*this->cst_b[l]/this->Bm/this->Bm)/rt/2.8284 * std::log((this->Z+2.4142*this->Bm*this->P/rt)/(this->Z-.4142*this->Bm*this->P/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm)
          - this->Am / this->Bm / rt / 2.8284 * ((this->P/rt*this->dVdN[l] + 2.4142*this->cst_b[l]*this->P/rt)/(this->Z+2.4142*this->Bm*this->P/rt) - (this->P/rt*this->dVdN[l] - 0.4142*this->cst_b[l]*this->P/rt)/(this->Z-.4142*this->Bm*this->P/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm)
          - this->Am / this->Bm / rt / 2.8284 * std::log((this->Z+2.4142*this->Bm*this->P/rt)/(this->Z-.4142*this->Bm*this->P/rt)) * (2.*this->A_IJ[apos]/this->Am - this->dAmdN[k]/this->Am/this->Am*this->dAmdN[l] + this->cst_b[k]/this->Bm/this->Bm*this->cst_b[l]));
    }
    dmuKdNL[k][k] += rt  / this->X[k];
  }
  // comment: what dimension is dmuKdNL?
  // answer: J/kg
  // can I change it to J/mol, which is the true definition of the chemical potential?
  // DEBUG
  {
    LOOP_k_N(this->n_species) {
      LOOP_l_N(this->n_species) {
        printf("dmu%ddN%d = %.4e\n", k, l, dmuKdNL[k][l]);
      }
    }
    throw(-1);
  }
  /*
     (-2.060658076895155*this->Am*this->nTotal*(-0.48528137423698625))/ ((-2.4142135622162813 + this->Bm* n0 + this->Bm *n1)*(0.4142135623684789 + this->Bm*n0 + this->Bm*n1)) - 
     (this->Bm*(n0 + n1) * rt)/(-1 + this->Bm*(n0 + n1)) + rt *std::log(n0) - rt *Log(1 - this->Bm*(n0 + n1)) + (0.3535530311869129*this->Am*std::log((1 - 0.4142135624*this->Bm*(n0 + n1))/(1 + 2.4142135624*this->Bm*(n0 + n1))))/this->Bm
     */


  double dmuFdXF = dmuKdNL[0][0] - dmuKdNL[0][1];
  double dmuOdXF = dmuKdNL[1][0] - dmuKdNL[1][1];
  double dmuFdXO = - dmuFdXF;
  double dmuOdXO = - dmuOdXF;

  double dmuFdYF = dmuFdXF / dYFdXF;
  double dmuOdYF = dmuOdXF / dYFdXF;
  double dmuFdYO = - dmuFdYF;
  double dmuOdYO = - dmuOdYF;
  double dmuFdYO_test = dmuFdXO / dYOdXO;
  double dmuOdYO_test = dmuOdXO / dYOdXO;

  if (std::fabs((dmuFdYO - dmuFdYO_test) / dmuFdYO) > 1e-6) {
    std::cout << "Error 21: dmufdyo is weird" << std::endl;
    throw(-1);
  }
  if (std::fabs((dmuOdYO - dmuOdYO_test) / dmuOdYO) > 1e-6) {
    std::cout << "Error 22: dmuodyo is weird" << std::endl;
    throw(-1);
  }

  double dmuFdrhoYF = dmuFdYF / denom;
  double dmuOdrhoYF = dmuOdYF / denom;
  double dmuFdrhoYO = dmuFdYO / drhoYOdYO;
  double dmuOdrhoYO = dmuOdYO / drhoYOdYO;

  // Let's do finite difference here.
  double curr_rho = this->rho;
  double curr_rhoYF = this->rho * this->Y[0];
  double delta_rhoYF = 0.01 * curr_rhoYF;
  double curr_rhoYO = this->rho * this->Y[1];
  double delta_rhoYO = 0.001 * curr_rhoYO;
  double rhoYF_plus = curr_rhoYF + delta_rhoYF;
  double rhoYF_minus = curr_rhoYF - delta_rhoYF;
  double rhoYO_plus = curr_rhoYO + delta_rhoYO;
  double rhoYO_minus = curr_rhoYO - delta_rhoYO;

  // (1) d(rhoYF)

  this->rho = this->GetRho_SetMixture_TPRY(this->T, this->Pamb, rhoYF_plus, this->rho);
  this->Y[0] = rhoYF_plus / this->rho;
  this->Y[1] = 1.0 - this->Y[0];
  double Y_plus [2] = {this->Y[0], this->Y[1]};
  this->SetMixture_TRY(this->T, this->rho, Y_plus);
  double muF_plus = this->muspecies[0];
  double muO_plus = this->muspecies[1];

  this->rho = this->GetRho_SetMixture_TPRY(this->T, this->Pamb, rhoYF_minus, this->rho);
  this->Y[0] = rhoYF_minus / this->rho;
  this->Y[1] = 1.0 - this->Y[0];
  double Y_minus [2] = {this->Y[0], this->Y[1]};
  this->SetMixture_TRY(this->T, this->rho, Y_minus);
  double muF_minus = this->muspecies[0];
  double muO_minus = this->muspecies[1];

  // fd 
  double fd_dmuFdrhoYF = (muF_plus - muF_minus) / (rhoYF_plus - rhoYF_minus);
  double fd_dmuOdrhoYF = (muO_plus - muO_minus) / (rhoYF_plus - rhoYF_minus);
  /*
     printf("true dmuFdrhoYF = %.4e\n", dmuFdrhoYF);
     printf("fd dmuFdrhoYF = %.4e\n", fd_dmuFdrhoYF);
     printf("true dmuOdrhoYF = %.4e\n", dmuOdrhoYF);
     printf("fd dmuOdrhoYF = %.4e\n", fd_dmuOdrhoYF);
     */

  // (2) d(rhoYO)

  // plus
  this->rho = curr_rhoYF + rhoYO_plus;
  this->Y[1] = rhoYO_plus / this->rho;
  this->Y[0] = 1.0 - this->Y[1];
  Y_plus[0] = this->Y[0];
  Y_plus[1] = this->Y[1];
  this->SetMixture_TRY(this->T, this->rho, Y_plus);
  muF_plus = this->muspecies[0];
  muO_plus = this->muspecies[1];

  // minus
  this->rho = curr_rhoYF + rhoYO_minus;
  this->Y[1] = rhoYO_minus / this->rho;
  this->Y[0] = 1.0 - this->Y[1];
  Y_minus[0] = this->Y[0];
  Y_minus[1] = this->Y[1];
  this->SetMixture_TRY(this->T, this->rho, Y_minus);
  muF_minus = this->muspecies[0];
  muO_minus = this->muspecies[1];

  // fd 
  double fd_dmuFdrhoYO = (muF_plus - muF_minus) / (rhoYO_plus - rhoYO_minus);
  double fd_dmuOdrhoYO = (muO_plus - muO_minus) / (rhoYO_plus - rhoYO_minus);

  /*
     printf("true dmuFdrhoYO = %.4e\n", dmuFdrhoYO);
     printf("fd dmuFdrhoYO = %.4e\n", fd_dmuFdrhoYO);
     printf("true dmuOdrhoYO = %.4e\n", dmuOdrhoYO);
     printf("fd dmuOdrhoYO = %.4e\n", fd_dmuOdrhoYO);
  // TODO: FD and calculation of dmuFdrhoYO is not the same...
  throw(-1);
  */

  dmuidrhoj[0][0] = dmuFdrhoYF;
  dmuidrhoj[0][1] = dmuFdrhoYO;
  dmuidrhoj[1][0] = dmuOdrhoYF;
  dmuidrhoj[1][1] = dmuOdrhoYO;
}

void PengRobinson::getdmuidnj(double** dmuidnj) {
  double** dmuidrhoj = new double* [2];
  dmuidrhoj[0] = new double[2];
  dmuidrhoj[1] = new double[2];
  this->getdmuidrhoj(dmuidrhoj);
  // dimension is in J/mol
  dmuidnj[0][0] = dmuidrhoj[0][0] * (this->MW[0] * 1e-3);
  dmuidnj[0][1] = dmuidrhoj[0][1] * (this->MW[1] * 1e-3);
  dmuidnj[1][0] = dmuidrhoj[1][0] * (this->MW[0] * 1e-3);
  dmuidnj[1][1] = dmuidrhoj[1][1] * (this->MW[1] * 1e-3);
  delete[] dmuidrhoj;
}

//----------------------------------------------------------------------------

void PengRobinson::SyncPartialPropertiesFromTemperatureDensity() {
  // hspecies_ig is molar, without rt
  double temp = this->Am - this->T * this->dAmdT;
  double R_u = this->gas_constant_universal;
  double m_P = (this->P < 0) ? this->Pamb : this->P;
  double Bm_tilde = this->Bm * m_P / R_u / this->T;
  double Z = m_P / (this->rho * this->gas_constant * this->T);
  double rt = this->gas_constant_universal * this->T;

  LOOP_k_N(this->n_species) {
    this->edspecies[k] = this->dK1dN[k] * temp +
      this->K1 * (this->dAmdN[k] - this->T * this->d2AmdTdN[k]);
    this->hspecies[k] = this->hspecies_ig[k] * rt - rt + this->edspecies[k] +
      m_P * dVdN[k];
    this->hspecies[k] /= this->MW[k];
    this->edspecies[k] /= this->MW[k];
  }

  // muspecies, sspecies
  double temp1;
  LOOP_k_N(this->n_species) {
    int i = k;
    double sum = 0;
    LOOP_l_N (this->num_species) {
      int j = l;
      int apos = j * this->num_species + i;
      double Aji = this->cst_a[apos] * std::pow(1.0 + this->cst_c[apos] * (1.0 - std::sqrt(this->T / this->Tcrit_IJ[apos])), 2);
      sum += this->X[j] * Aji;
    }
    // (M.3.18)
    // double logphik = - log(Z - Bm_tilde)
    //   + this->cst_b[k] / this->Bm * (Z-1.0)
    //   - log((this->delta1 * Bm_tilde +  Z)/(this->delta2 * Bm_tilde + Z))
    //   * this->Am / this->Bm / R_u / this->T / (this->delta1 - this->delta2)
    //   * (2. * sum / this->Am - this->cst_b[k] / this->Bm);
    // if (std::isnan(logphik) || std::isinf(logphik)) {
    //   std::cerr << "Error: PengRobinson::SyncPartialPropertiesFromTemperatureDensity: "
    //             << "lnPhi[" << k << "] is NaN or Inf." << std::endl;
    //   printf("k = %d, rho = %.4f, Z = %.4f, Bm_tilde = %.4f, MW = %.4f, T = %.4f, P = %.4f, m_P = %.4f\n",
    //          k, this->rho, Z, Bm_tilde, this->MW_M, this->T, this->P, m_P);
    //   exit(-1);
    // }
    double logphik = this->cst_b[k] / this->Bm * (this->Z - 1.) - std::log(this->ZForDiff - this->Bm*this->PForDiff/rt)
      - this->Am / this->Bm / rt / 2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm);
    /*
       double logphik = this->cst_b[k] / this->Bm * (this->Z - 1.) - std::log(this->Z - this->Bm*this->P/rt)
       - this->Am / this->Bm / rt / 2.8284 * std::log((this->Z+2.4142*this->Bm*this->P/rt)/(this->Z-.4142*this->Bm*this->P/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm);
       */

    this->muspecies[k] = this->hspecies_ig[k]*rt - this->sspecies_ig[k]*rt + rt * logphik;
    this->muspecies[k] /= this->MW[k];
    this->sspecies[k] = (this->hspecies[k] - this->muspecies[k]) / this->T;
    this->muspecies_molar[k] = 0.0;
  }
  if (std::isnan(muspecies[0])) {
    std::cout << "Something wrong Partial. T = " << this->T << ", rho = " << this->rho
      << ", YF = " << this->Y[0] << " , P = " << this->P << ", " << this->ZForDiff - this->Bm*this->PForDiff/rt << ", " << (this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt) << ", " << temp1
      << ", " << this->cst_b[0] / this->Bm * (this->Z - 1.) - std::log(this->ZForDiff - this->Bm*this->PForDiff/rt)
      << ", " << this->Am / this->Bm / rt / 2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[0]/this->Am - this->cst_b[0]/this->Bm)
      << ", " << this->cst_b[0] / this->Bm * (this->Z - 1.) - std::log(this->ZForDiff - this->Bm*this->PForDiff/rt)
      - this->Am / this->Bm / rt / 2.8284 * std::log((this->ZForDiff+2.4142*this->Bm*this->PForDiff/rt)/(this->ZForDiff-.4142*this->Bm*this->PForDiff/rt)) * (this->dAmdN[0]/this->Am - this->cst_b[0]/this->Bm) << std::endl;
    throw -1;
  }
  return;


  // alpha_D
  if (this->n_species == 2) {
    int i = 0, j = 0; // We can use a relation to speed up, so only need to calculate alpha_D[0,0]
    int apos = i * this->n_species + j;

    double A = this->Am * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
    double dAdXi = this->dAmdN[i] * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
    double dAdXj = this->dAmdN[j] * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
    double B = this->Bm * this->P / this->gas_constant_universal / this->T;
    double dBdXj = this->cst_b[j] * this->P / this->gas_constant_universal / this->T;
    double Z = this->Z;
    double dZdXj = this->dVdN[j] * this->P / this->gas_constant_universal / this->T;

    double A_IJ[this->n_species * this->n_species];
    LOOP_k_N(this->n_species  * this->n_species)
      A_IJ[k] = this->A_IJ[k] * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
    double B_I[this->n_species];
    LOOP_k_N(this->n_species)
      B_I[k] = this->cst_b[k] * this->P / this->gas_constant_universal / this->T;

    double dlogPhiidXj = B_I[i] / B * dZdXj - B_I[i] * Z / B / B * dBdXj
      - 1. / (Z - B) * (dZdXj - dBdXj)
      - (dAdXj / B / std::sqrt(8.) * std::log((Z + (1.+std::sqrt(2.)) * B) / (Z + (1.-std::sqrt(2.)) * B)) * (dAdXi / A - B_I[i] / B)
          - A / B / B * dBdXj / std::sqrt(8.) * std::log((Z + (1.+std::sqrt(2.)) * B) / (Z + (1.-std::sqrt(2.)) * B)) * (dAdXi / A - B_I[i] / B)
          + A / B / std::sqrt(8.) / ((Z + (1.+std::sqrt(2.)) * B) / (Z + (1.-std::sqrt(2.)) * B)) * ((dZdXj + (1.+std::sqrt(2.)) * dBdXj) / (Z + (1.-std::sqrt(2.)) * B)
            - (Z + (1.+std::sqrt(2.)) * B)
            / (Z + (1.-std::sqrt(2.)) * B) / (Z + (1.-std::sqrt(2.)) * B)
            * (dZdXj + (1.-std::sqrt(2.)) * dBdXj))
          + A / B / std::sqrt(8.) * std::log((Z + (1.+std::sqrt(2.)) * B) / (Z + (1.-std::sqrt(2.)) * B)) * (2. * A_IJ[apos] / A
            - dAdXi / A / A * dAdXj
            + B_I[i] / B / B * dBdXj));

    this->alphaD_IJ[0] = dlogPhiidXj * this->X[i] + 1;
    // The Gibbs-Duhem relation for binary mixture
    this->alphaD_IJ[3] = this->alphaD_IJ[0];
    this->alphaD_IJ[1] = -this->alphaD_IJ[0];
    this->alphaD_IJ[2] = -this->alphaD_IJ[0];
  } else if (this->n_species == 1) {

  } else {
    std::cout << "Not implemented: alpha_D for more than two species" << std::endl;
    throw -1;
  }
}

//----------------------------------------------------------------------------

void PengRobinson::SyncTransportFromTemperatureDensity() {
  // ChungHP for viscosity and conductivity
  double sigma_M = 0.0;
  double epsOverk_M = 0.0;
  double omega_M = 0.0;
  double omega_M_visc = 0.0;
  double mu_M = 0.0;
  double kappa_M = 0.0;
  double MW_M_visc = 0.0;

  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * n_species + l;
      double X_X = this->X[l] * this->X[k];
      sigma_M += X_X * std::pow(this->sigma_IJ[apos], 3.0);
      MW_M_visc += X_X * epsOverk_IJ[apos] * std::pow(this->sigma_IJ[apos], 2)
        * sqrt(this->MW_IJ[apos]);
      epsOverk_M += X_X * epsOverk_IJ[apos] * std::pow(this->sigma_IJ[apos], 3);
      omega_M += X_X * this->omega_IJ[apos] * std::pow(this->sigma_IJ[apos], 3.0);
      omega_M_visc +=
        X_X * this->omega_visc_IJ[apos] * std::pow(this->sigma_IJ[apos], 3.0);
      if (fabs(this->sigma_IJ[apos]) > 1.0e-20) {
        mu_M += X_X * std::pow(this->dipole[k] * this->dipole[l], 2.0)
          / std::pow(this->sigma_IJ[apos], 3.0);
      }
    }
  }
  sigma_M = std::pow(sigma_M, 1.0 / 3.0);
  epsOverk_M = epsOverk_M / std::pow(sigma_M, 3);
  omega_M = omega_M / std::pow(sigma_M, 3.0);
  omega_M_visc = omega_M_visc / std::pow(sigma_M, 3.0);
  mu_M = std::pow(mu_M * std::pow(sigma_M, 3.0), 0.25);
  MW_M_visc = std::pow(MW_M_visc / (epsOverk_M * std::pow(sigma_M, 2)), 2);

  double VCrit_M = std::pow((sigma_M / 0.809), 3.0);
  double TCrit_M = 1.2593 * epsOverk_M;
  double TStar_M = this->T / epsOverk_M;
  double mu_RM = 131.3 * mu_M / std::sqrt(VCrit_M * TCrit_M);
  double Neufeld =
    1.16145 * std::pow(TStar_M, -0.14874) + 0.52487 * std::exp(-0.7732 * TStar_M)
    + 2.16178 * std::exp(-2.43787 * TStar_M); // Eq (9-4-3) in Poling
  double
    Fc = 1.0 - 0.2756 * omega_M_visc + 0.059035 * std::pow(mu_RM, 4) + kappa_M;

  double mu = 40.785 * Fc * std::sqrt(MW_M_visc * this->T)
    / (Neufeld * std::pow(VCrit_M, 2.0 / 3.0))
    / 10.0e6; // used for Chung_HP conductivity

  const double a1[10] =
  {6.324, 1.210e-3, 5.283, 6.623, 19.745, -1.900, 24.275, 0.7972, -0.2382,
    0.06863};
  const double b1[10] =
  {50.412, -1.154e-3, 254.209, 38.096, 7.630, -12.537, 3.450, 1.1170,
    0.0677, 0.3479};
  const double c1[10] =
  {-51.680, -6.257e-3, -168.480, -8.464, -14.354, 4.985, -11.291, 0.01235,
    -0.8163, 0.5926};
  const double d1[10] =
  {1189.000, 0.03728, 3898.0, 31.42, 31.53, -18.15, 69.35, -4.117, 4.025,
    -0.727};
  double E[10];
  LOOP_k_N(10)
    E[k] = a1[k] + b1[k] * omega_M_visc + c1[k] * std::pow(mu_RM, 4.0)
    + d1[k] * kappa_M;

  double y = this->rho / (1000 * MW_M_visc) * VCrit_M / 6.0;

  double G1 = (1.0 - 0.5 * y) / std::pow((1.0 - y), 3);
  double G2 = (E[0] / y * (1.0 - std::exp(-E[3] * y)) + E[1] * G1 * std::exp(E[4] * y)
      + E[2] * G1) / (E[0] * E[3] + E[1] + E[2]);

  double mustarstar = E[6] * std::pow(y, 2) * G2
    * std::exp(E[7] + E[8] / TStar_M + E[9] / std::pow(TStar_M, 2));
  double mustar =
    std::pow(TStar_M, 0.5) / Neufeld * (Fc * (1.0 / G2 + E[5] * y)) + mustarstar;

  double mixViscosity =
    (36.644 * mustar * std::sqrt(MW_M_visc * TCrit_M) / std::pow(VCrit_M, 2.0 / 3.0))
    / 10.0e6;
  // conductivity
  const double
    a2[7] = {2.4166, -0.50924, 6.6107, 14.543, 0.79274, -5.8634, 91.089};
  const double
    b2[7] = {0.74824, -1.5094, 5.6207, -8.9139, 0.82019, 12.801, 128.11};
  const double
    c2[7] = {-0.91858, -49.991, 64.760, -5.6379, -0.69369, 9.5893, -54.217};
  const double
    d2[7] = {121.72, 69.983, 27.039, 74.344, 6.3173, 65.529, 523.81};
  // compute B_i coefficients
  double B_coef[7];
  LOOP_k_N(7)
    B_coef[k] = a2[k] + b2[k] * omega_M + c2[k] * std::pow(mu_RM, 4.0)
    + d2[k] * kappa_M;

  G2 = ((B_coef[0] / y) * (1.0 - exp(-B_coef[3] * y))
      + B_coef[1] * G1 * exp(B_coef[4] * y) + B_coef[2] * G1)
    / (B_coef[0] * B_coef[3] + B_coef[1] + B_coef[2]);

  double cv_ig = this->cp_ig * this->MW_M
    - this->gas_constant_universal;

  double alpha = cv_ig / this->gas_constant_universal - 3.0 / 2.0;
  double beta = 0.7862 - 0.7109 * omega_M + 1.3168 * std::pow(omega_M, 2);
  double Zr = 2.0 + 10.5 * std::pow(this->T / TCrit_M, 2);
  double q = 0.003586 * std::sqrt(TCrit_M / (MW_M_visc / 1000.0))
    / std::pow(VCrit_M, 2.0 / 3.0);
  double psi = 1.0 + alpha
    * ((0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * Zr)
        / (0.6366 + beta * Zr + 1.061 * alpha * beta));

  double mixConductivity =
    31.2 * mu * psi / (MW_M_visc / 1000) * (1.0 / G2 + B_coef[5] * y)
    + q * B_coef[6] * std::pow(y, 2) * sqrt(this->T / TCrit_M) * G2;

  this->mu = mixViscosity;
  this->lambda = mixConductivity;


  // diffusion

  double a[4] = {1.0548, 0.1550, 0.55909, 2.1705};

  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * n_species + l;
      if (k == l)
        this->Dm[apos] = 0.;
      else {
        double T_star = this->T / this->epsOverk_IJ[apos];
        double Omega = a[0] * std::pow(T_star, -a[1]) + std::pow(T_star + a[2], -a[3]);

        this->Dm[apos] = 3. / 16. * std::sqrt(2. * 3.1415 * std::pow(this->Boltzman * this->T, 3.) / this->m_IJ[apos]) / (
            this->PForDiff * 3.1415 * std::pow(this->sigma_diff_IJ[apos], 2) * Omega);
      }
    }
  }

  // Onsager
  if (this->n_species > 2) {
    std::cout << "Onsager relations not implemented" << std::endl;
    throw -1;
  } else if (this->n_species == 2) {
    this->Lij[0] = this->rho * this->Dm[1] * this->MW[0] * this->MW[1] * this->X[0] * (1.-this->X[0])
      / this->gas_constant_universal / this->MW_M;
    double lv = 4.5;
    double Qstar[this->n_species];

    double sum1 = 0.;
    double totalV = 0.;
    LOOP_l_N(this->n_species) {
      double dVdN_t = this->dVdN[l];
      sum1 += this->X[l] * this->edspecies[l]*this->MW[l] / lv;
      totalV += this->X[l] * dVdN_t;
    }

    LOOP_k_N(this->n_species) {
      double dVdN_t = this->dVdN[k];
      Qstar[k] = -this->edspecies[k]*this->MW[k] / lv + sum1 * dVdN_t / totalV;
    }
    this->Lqi[0] = this->Lij[0] * (Qstar[0]/this->MW[0] - Qstar[1]/this->MW[1]);
    this->kT[0] = this->Lqi[0] / this->rho / this->T / this->Dm[1];
    // From here assume 2 species
    // Capillary coefficient
    double kappa_species[this->n_species];

    LOOP_k_N(this->n_species) {
      double tr = this->T / this->Tcrit[k];
      double t = 1. - tr;
      //    if (t<0.01)
      //      t=0.01;
      if (t <= 0.) {
        kappa_species[k] = 0.;
      } else {
        double temp = std::log(t);
        double nonDimKappa = std::exp(this->k0[k] + this->k1[k] * temp + this->k2[k] * temp * temp);
        //      if (nonDimKappa > 1.) {
        //        std::cout << "Error in nonDimKappa: " << nonDimKappa << std::endl;
        //        throw -1;
        //      }
        int apos = k * this->n_species + k;
        kappa_species[k] = nonDimKappa * this->A_IJ[apos] * std::pow(this->cst_b[k] / this->N_Av, 2./3.);
      }
      double nonDimKappa = (this->Akappa[k] * t + this->Bkappa[k])/100.*0.5*(1.+std::tanh((t+1.)/0.01));
      int apos = k * this->n_species + k;
      kappa_species[k] = nonDimKappa * this->A_IJ[apos] * std::pow(this->cst_b[k], 2./3.);

      /*
       * [2.7959074445220687e-20, 4.8782874932788375e-20]
       * [1.634272050813004e-20, 3.1826209821406494e-21]
       */
      // Special cases
      if (species[k] == "H2O") {
        // kappa_species[k] = 1.3989e-14; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
        // kappa_species[k] = 1.3989e-20; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
        // kappa_species[k] = 1.36e-14;
        // kappa_species[k] = 1.634272e-14;
        double alpha = - 1e-16 / (1.2326 + 1.3757 * this->omega[k]);
        double beta = 1e-16 / (0.9051 + 1.5410 * this->omega[k]);
        double cij = this->cst_a[k] * 1e-6 * std::pow(this->cst_b[k]*1e-3, 2./3.) * (alpha * (1. - tr) + beta); // [J*m^5/mol^2]
        // kappa_species[k] = 1.3989e-14; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
                                       // kappa_species[k] = 1.3989e-20; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
        kappa_species[k] = cij * 1e6; // [J*m^5/kmol^2]
                                      //
        // kappa_species[k] = 1.3989e-14; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
      } else if (species[k] == "CH4") {
        // kappa_species[k] = 2.0387e-20; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
        kappa_species[k] = 2.0387e-14; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
      } else if (species[k] == "N2") {
        kappa_species[k] = 3.18262e-15;
      }
      this->kappa_spec[k] = kappa_species[k];
      this->c_spec[k] = kappa_species[k] * 1e-6;
    }
    kappa_species[1] = 0.; // ?
    LOOP_k_N(this->n_species) {
      LOOP_l_N(this->n_species) {
        int apos = k * this->n_species + l;
        this->kappa[apos] = std::sqrt(kappa_species[k] * kappa_species[l]);
      }
    }
  } else { // this->n_species == 1
           // Capillary coefficient
    double kappa_species[this->n_species];
    LOOP_k_N(this->n_species) {
      double tr = this->T / this->Tcrit[k];
      double t = 1. - tr;
      if (t <= 0.) {
        kappa_species[k] = 0.;
      } else {
        double temp = std::log(t);
        double nonDimKappa = std::exp(this->k0[k] + this->k1[k] * temp + this->k2[k] * temp * temp);
        int apos = k * this->n_species + k;
        kappa_species[k] = nonDimKappa * this->A_IJ[apos] * std::pow(this->cst_b[k] / this->N_Av, 2./3.);
      }
      double nonDimKappa = (this->Akappa[k] * t + this->Bkappa[k])/100.*0.5*(1.+std::tanh((t+1.)/0.01));
      int apos = k * this->n_species + k;
      kappa_species[k] = nonDimKappa * this->A_IJ[apos] * std::pow(this->cst_b[k], 2./3.);

      // Special cases
      if (species[k] == "H2O") {
        // kappa_species[k] = 1.3989e-14; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
        // kappa_species[k] = 1.3989e-20; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
        // kappa_species[k] = 1.36e-14;
        // kappa_species[k] = 1.634272e-14;
        double alpha = - 1e-16 / (1.2326 + 1.3757 * this->omega[k]);
        double beta = 1e-16 / (0.9051 + 1.5410 * this->omega[k]);
        double cij = this->cst_a[k] * 1e-6 * std::pow(this->cst_b[k]*1e-3, 2./3.) * (alpha * (1. - tr) + beta); // [J*m^5/mol^2]
        // kappa_species[k] = 1.3989e-14; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
                                       // kappa_species[k] = 1.3989e-20; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
        kappa_species[k] = cij * 1e6; // [J*m^5/kmol^2]
      } else if (species[k] == "CH4") {
        // kappa_species[k] = 2.0387e-20; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
        kappa_species[k] = 2.0387e-14; // Miqueue, Journal of Physical Chemistry B, 2011, p. 9620
      }
    }
    kappa_species[1] = 0.; // ?
                           // TODO: pure fluid expression, there is a prefactor beta 
    LOOP_k_N(this->n_species) {
      LOOP_l_N(this->n_species) {
        int apos = k * this->n_species + l;
        this->kappa[apos] = std::sqrt(kappa_species[k] * kappa_species[l]);
      }
    }
  }
}

//----------------------------------------------------------------------------

bool PengRobinson::CheckThreeRoots() {
  double a = this->Am;
  double b = this->Bm;

  double A = a * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
  double B = b * this->P / this->gas_constant_universal / this->T;

  double a0 = -(A * B - B * B - B * B * B);
  double a1 = A - 3 * B * B - 2 * B;
  double a2 = -(1 - B);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; i++)
    if (xZ[i] > B)
      Z.push_back(xZ[i]);

  return Z.size() > 1;
}

void PengRobinson::ComputeMixtureCrit(double &Tcmix, double &Pcmix, double &wcmix, double &Zcmix) {
  Tcmix = 0.;
  Pcmix = 0.;
  wcmix = 0.;
  Zcmix = 0.;
  double Vcmix = 0.;

  LOOP_k_N(this->n_species) {
    Tcmix += this->X[k] * this->Tcrit[k];
    Zcmix += this->X[k] * this->Zcrit[k];
    Vcmix += this->X[k] * this->Vcrit[k];
    wcmix += this->X[k] * this->omega[k];
  }
  Pcmix = Zcmix * this->gas_constant_universal * Tcmix / Vcmix;
}

//----------------------------------------------------------------------------



//----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, const PengRobinson& thermal) {
  os << "--------------------------------------------------------------------------------" << std::endl;
  os << "Number of species: " << thermal.n_species << " ";
  LOOP_k_N(thermal.n_species)
    os << thermal.species[k] << " ";
  os << std::endl;

  os  << "P = " << thermal.P << std::endl
    << "T = " << thermal.T << std::endl
    << "rho = " << thermal.rho << std::endl
    << "Z = " << thermal.Z << std::endl
    << "X = ";
  LOOP_k_N(thermal.n_species)
    os << thermal.X[k] << " ";
  os << std::endl;
  os << "Y = ";
  LOOP_k_N(thermal.n_species)
    os << thermal.Y[k] << " ";
  os << std::endl;
  os  << "MWm = " << thermal.MW_M << std::endl
    << "E = " << thermal.E << std::endl
    << "Partial H = ";
  LOOP_k_N(thermal.n_species)
    os << thermal.hspecies[k] << " ";
  os << std::endl;
  os << "Partial mu = ";
  LOOP_k_N(thermal.n_species)
    os << thermal.muspecies[k] << " ";
  os << std::endl;
  os << "Partial s = ";
  LOOP_k_N(thermal.n_species)
    os << thermal.sspecies[k] << " ";
  os << std::endl;
  os << "Mass Diffusivity = ";
  LOOP_k_N(thermal.n_species * thermal.n_species)
    os << thermal.alphaD_IJ[k] << " ";
  os << std::endl;
  os << "Thermal Diffusivity = ";
  LOOP_k_N(thermal.n_species * thermal.n_species)
    os << thermal.alphaT[k] << " ";
  os << std::endl;
  os << "Binary Diffusivity = ";
  LOOP_k_N(thermal.n_species * thermal.n_species)
    os << thermal.Dm[k] << " ";
  os << std::endl;
  os  << "Kappa = ";
  LOOP_k_N(thermal.n_species * thermal.n_species)
    os << thermal.kappa[k] << " ";
  os << std::endl;
  os  << "mu = " << thermal.mu << std::endl
    << "lambda = " << thermal.lambda << std::endl;
  os << "Lqi = ";
  LOOP_k_N(thermal.n_species-1)
    os << thermal.Lqi[k] << " ";
  os << std::endl;
  os << "kT = ";
  LOOP_k_N(thermal.n_species-1)
    os << thermal.kT[k] << " ";
  os << std::endl;
  os << "Lij = ";
  LOOP_k_N((thermal.n_species-1)*(thermal.n_species-1))
    os << thermal.Lij[k] << " ";
  os << std::endl;
  os << "sos = " << thermal.sos << std::endl;
  os << "gamma = " << thermal.gamma << std::endl;
  os << "--------------------------------------------------------------------------------" << std::endl;
  return os;
}


// To be added..
void PengRobinson::GetRhoLRhoV(const double T_in, const double P_in, double & rhoL, double & rhoV) {
  double Tout;
  double Xvapor [this->n_species];
  double Xliq [this->n_species];
  double Xvapor_guess [this->n_species];
  double Xliq_guess [this->n_species];
  // TPflash
  this->TPFlash(T_in, P_in, Xvapor_guess, Xliq_guess, Xvapor, Xliq, Tout);
  rhoL = this->GetRho_SetMixture_TPX(T_in, P_in, Xliq);
  rhoV = this->GetRho_SetMixture_TPX(T_in, P_in, Xvapor);
}

double PengRobinson::GetRhoFromN(const double* n_in) {
  double rho_i [this->n_species];
  double Y_i [this->n_species];
  double rho_ = 0.0;
  LOOP_k_N(this->n_species) {
    rho_i[k] = n_in[k] * this->GetMW()[k];
    rho_ += rho_i[k];
  }
  return rho_;
}

void PengRobinson::SyncMolarDensities() {
  this->molar_density_total = 0.0;
  LOOP_k_N (this->n_species) {
    double rhoYF_ = this->rho * this->Y[k];
    this->molar_density[k] = rhoYF_ / (this->MW[k]);
    this->molar_density_total += this->molar_density[k];
  }
  this->SyncHelmholtzEnergy();
  Phi = this->helmholtz_energy * this->rho;
  LOOP_k_N (this->n_species) {
    Phi -= this->molar_density[k] * this->muspecies[k] * (this->MW[k]);
  }
}

void PengRobinson::SyncHelmholtzEnergy() {
  double rt = this->gas_constant_universal * this->T;
  double f_ideal = 0.0;
  double f_ex1 = - this->molar_density_total * rt 
    * std::log(1. - this->Bm * this->molar_density_total);
  double f_ex2 = this->Am * this->molar_density_total / (std::sqrt(8) * this->Bm)
    * std::log((1. + (1. - std::sqrt(2)) * this->Bm * this->molar_density_total)
        /(1. + (1. + std::sqrt(2)) * this->Bm * this->molar_density_total));
  LOOP_k_N (this->n_species) {
    f_ideal += this->molar_density[k] * (std::log(molar_density[k]) - 1.0);
  }
  f_ideal *= rt;
  this->helmholtz_energy = f_ideal + f_ex1 + f_ex2;
}

// ----------------------------------------------------------------------------
// @\brief: compute densities, energies, and fugacities of a phase: Algorithm 1
// @param T: temperature
// @param P: pressure
// @param chi: molar fractions of phase
// @param phase: phase type (iLIQPH_TP, iVAPPH_TP)
// ----------------------------------------------------------------------------
void PengRobinson::generateVLstructure(const double T, const double P, const double* chi, const int phase) {
  this->phase = phase;
  // for copying later
  double X_copy [this->num_species];
  LOOP_k_N(this->num_species) X_copy[k] = this->X[k];

  // I have to save the global molar fraction, so that later I can save it, but now I have to replace X with x_TP.
  if (phase == iLIQPH_TP) {
    // set phase mass fraction
    LOOP_k_N(this->num_species) { 
      this->x_TP[k] = chi[k];
    }
    this->SetMolarFraction(chi);
    this->MW_L = this->MW_M; 
    this->gas_constant = this->gas_constant_universal / this->MW_L;
    // compute phase compressibility
    // compute phase density
    this->rho_L = this->GetDensityFromPressureTemperature_TP(P, T, phase);
    this->rhoYF_L = this->rho_L * this->Y[0];
    REPLACEIFWEIRD(this->rhoYF_L);
    REPLACEIFWEIRD(this->rho_L);
    this->Z_L = P * this->MW_L / (this->rho_L * this->gas_constant_universal * T);
    REPLACEIFWEIRD(this->Z_L)
    // compute phase energy
    this->E_L = this->GetEnergyFromTemperature(T, true); // molar (true), J/kmol
    // compute phase fugacity coefficient and fugacity
    this->computeFugacityCoefficient(this->phi_L, this->f_L, this->iLIQPH_TP);
  } else if (phase == iVAPPH_TP) {
    // set phase mass fraction
    LOOP_k_N(this->num_species) { 
      this->y_TP[k] = chi[k];
    }
    this->SetMolarFraction(chi);
    this->MW_V = this->MW_M;
    this->gas_constant = this->gas_constant_universal / this->MW_V;
    // compute phase compressibility
    // compute phase density
    this->rho_V = this->GetDensityFromPressureTemperature_TP(P, T, phase);
    this->rhoYF_V = this->rho_V * this->Y[0];
    REPLACEIFWEIRD(this->rhoYF_V);
    REPLACEIFWEIRD(this->rho_V);
    this->Z_V = P * this->MW_V / (this->rho_V * this->gas_constant_universal * T);
    REPLACEIFWEIRD(this->Z_V)
    // compute phase energy
    this->E_V = this->GetEnergyFromTemperature(T, true); // molar (true), J/kmol
    // compute phase fugacity coefficient and fugacity
    this->computeFugacityCoefficient(this->phi_V, this->f_V, this->iVAPPH_TP);
  } else {
    std::cerr << "Error: RealFluid::generateVLstructure() - "
              << "phase must be either iLIQPH_TP or iVAPPH_TP." << std::endl;
    throw (-1);
  }

  this->SetMolarFraction(X_copy); // restore the global molar fraction
  this->gas_constant = this->gas_constant_universal / this->MW_M;
} // end RealFluid::generateVLstructure
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @brief: compute simple (Am, Bm, Z) temperature and molar derivatives for the real fluid
// check Table 4-10
// ----------------------------------------------------------------------------
void PengRobinson::computeBasicDerivatives(const int phase) {
  this->dAmdT = 0.0; this->d2AmdT2 = 0.0;
  double c2 = 0.; double c3 = 0.; double c5 = 0.;
  double dc1dXk [this->num_species];
  double dc2dXk [this->num_species];
  double dc3dXk [this->num_species];
  double dc5dXk [this->num_species];
  if (this->mixing_rule == "VDW") {
    LOOP_k_N (this->num_species) {
      // Table 10 for dAmdXk
      this->dAmdXk[k] = 0.0;
      // zero-out the derivatives
      dc1dXk[k] = 0.0;
      dc2dXk[k] = 0.0;
      dc3dXk[k] = 0.0;

      double A_i = this->cst_a[k*k];
      double f_i = this->cst_c[k*k];
      double sqrtTci = sqrt(this->Tcrit[k]);
      LOOP_l_N (this->num_species) {
        double X_X = this->X[k] * this->X[l];
        double A_j = this->cst_a[l*l];
        double f_j = this->cst_c[l*l];
        double sqrtTcj = sqrt(this->Tcrit[l]);
        double sqrtAiAj = sqrt(A_i * A_j);
        c2 -= X_X * sqrtAiAj * (f_j * (1. + f_i) / sqrtTcj + f_i * (1. + f_j) / sqrtTci);
        c3 += X_X * sqrtAiAj * f_i * f_j / (sqrtTci * sqrtTcj);
        dc1dXk[k] += 2.0 * this->X[l] * sqrtAiAj * (1. + f_i) * (1. + f_j);
        dc2dXk[k] -= 2.0 * this->X[l] * sqrtAiAj * (f_j * (1. + f_i) / sqrtTcj + f_i * (1. + f_j) / sqrtTci);
        dc3dXk[k] += 2.0 * this->X[l] * sqrtAiAj * f_i * f_j / (sqrtTci * sqrtTcj);
      }
    }
  } else if (this->mixing_rule == "MR3") {
    LOOP_k_N (this->num_species) {
      this->dAmdXk[k] = 0.0;
      dc1dXk[k] = 0.0;
      dc2dXk[k] = 0.0;
      dc3dXk[k] = 0.0;
      dc5dXk[k] = 0.0;
    }
    // Table 10 for dAmdXk
    LOOP_k_N (this->num_species) {
      LOOP_l_N (this->num_species) {
        int i = l;
        int apos_ik = i * this->num_species + k;
        double X_X = this->X[i] * this->X[k];
        double f_ik = this->cst_c[apos_ik]; double Tcrit_ik = this->Tcrit_IJ[apos_ik];
        double prod = 1. + f_ik * (1. - sqrt(this->T/Tcrit_ik));
        double alpha_ik = prod * prod;
        double A_IK_star = this->cst_a[apos_ik];
        // Table 6
        c2 -= 2. * X_X * A_IK_star * f_ik * (1. + f_ik) / sqrt(Tcrit_ik);
        c3 += X_X * A_IK_star * f_ik * f_ik / Tcrit_ik;
        c5 -= 2. * X_X;
        // Table 7
        dc1dXk[k] += 2.0 * this->X[i] * A_IK_star * (1. + f_ik) * (1. + f_ik);
        dc2dXk[k] -= 4.0 * this->X[i] * A_IK_star * f_ik * (1. + f_ik) / sqrt(Tcrit_ik);
        dc3dXk[k] += 2.0 * this->X[i] * A_IK_star * f_ik * f_ik / Tcrit_ik;
        dc5dXk[k] -= 4.0 * this->X[i];
      }
    }
  }
  // save
  this->c2 = c5; this->c3 = c3; 
  LOOP_k_N (this->num_species) {
    this->dc2dXk[k] = dc5dXk[k];
    this->dc3dXk[k] = dc3dXk[k];
  }
  // Now actually compute the things based on Table 9
  double sqrtT = sqrt(this->T);
  this->dAmdT = 0.5 / sqrtT * c2 + c3;
  this->d2AmdT2 = - 0.25 * c2 / (this->T * sqrtT) ;
  LOOP_k_N (this->num_species) {
    this->dAmdXk[k] = dc1dXk[k] + sqrtT * dc2dXk[k] + this->T * dc3dXk[k];
    this->d2AmdTdXk[k] = 0.5 / sqrtT * dc2dXk[k] + dc3dXk[k];
  }

  // Now dZdXk; // check (M.1)
  double Z, m_rho, m_MW;
  if (phase == this->iLIQPH_TP) {
    m_rho = this->rho_L; Z = this->Z_L; m_MW = this->MW_L;
  } else if (phase == this->iVAPPH_TP) {
    m_rho = this->rho_V; Z = this->Z_V; m_MW = this->MW_V;
  } else {
    m_rho = this->rho; m_MW = this->MW_M; Z = this->P  * this->MW_M / (this->rho * this->gas_constant_universal * this->T);
  }

  double Z2 = Z * Z;
  double Z3 = Z2 * Z;
  double R_u = this->gas_constant_universal;
  double RuT_p = R_u * this->T / this->P;
  double Bm_tilde = this->Bm * this->P / R_u / this->T;
  double dBmtildedXk [this->num_species];
  LOOP_k_N (this->num_species) dBmtildedXk[k] = this->dBmdXk[k] / RuT_p; // (M.2.5)
  double da1dXk [this->num_species]; // zeros
  double da2dXk [this->num_species];
  double da3dXk [this->num_species];
  double da4dXk [this->num_species];
  double Ru2T2_p2 = RuT_p * RuT_p;
  double Ru3T3_p3 = Ru2T2_p2 * RuT_p;
  //(2.2.71)
  double a1 = Ru3T3_p3;
  double a2 = Ru3T3_p3 * (1/RuT_p * this->Bm * (this->delta1_p_delta2 - 1.0) - 1.0);
  double a3 = RuT_p * this->Bm * this->Bm * (this->delta1delta2 - this->delta1_p_delta2)
              - this->delta1_p_delta2 * Ru2T2_p2 * this->Bm + this->Am * RuT_p / this->P;
  double a4 = - this->delta1delta2 * RuT_p * this->Bm * this->Bm
              - this->Am * this->Bm / this->P - this->delta1delta2 * this->Bm * this->Bm * this->Bm;

  // temperature derivatives, (M.1.3)
  double da1dT = 3 * Ru3T3_p3 / this->T;
  double da2dT = Ru3T3_p3 / this->T * (2. * Bm_tilde * (this->delta1_p_delta2 - 1.0) - 3.0);
  double da3dT = R_u / this->P * this->Bm * this->Bm * (this->delta1delta2 - this->delta1_p_delta2)
                  - 2. * this->delta1_p_delta2 * Ru2T2_p2 / this->T * this->Bm
                  + this->Am * R_u / this->P / this->P
                  + RuT_p / this->P * this->dAmdT;
  double da4dT = - this->delta1delta2 * R_u / this->P * this->Bm * this->Bm
                 - this->Bm / this->P * this->dAmdT;
  // pressure derivatives, (M.1.3)
  double da1dp = -3 * Ru3T3_p3 / this->P;
  double da2dp = Ru3T3_p3 / this->P * (3. - 2. * Bm_tilde * (this->delta1_p_delta2 - 1.0));
  double da3dp = - RuT_p / this->P * this->Bm * this->Bm * (this->delta1delta2 - this->delta1_p_delta2)
                 + 2. * this->delta1_p_delta2 * Ru2T2_p2 / this->P * this->Bm
                 -2 * this->Am * RuT_p / this->P / this->P;
  double da4dp = this->delta1delta2 * RuT_p * this->Bm * this->Bm / this->P
                 + this->Am * this->Bm / this->P / this->P;

  double denom = 3. * a1 * Z2 + 2. * a2 * Z + a3;
  // M.1.5
  LOOP_k_N (this->num_species) {
    da1dXk[k] = 0.0;
    da2dXk[k] = Ru2T2_p2 * (this->delta1_p_delta2 - 1.0) * this->dBmdXk[k];
    da3dXk[k] = Ru2T2_p2 * (2. * Bm_tilde * (this->delta1delta2 - this->delta1_p_delta2) - this->delta1_p_delta2) * this->dBmdXk[k] + RuT_p / this->P * this->dAmdXk[k];
    da4dXk[k] = - (2. * this->delta1delta2 * RuT_p * this->Bm + this->Am / this->P + 3 * this->delta1delta2 * this->Bm * this->Bm) * this->dBmdXk[k] - this->Bm / this->P * this->dAmdXk[k];

    this->dZdXk[k] = - (Z3 * da1dXk[k] + Z2 * da2dXk[k]
                        + Z * da3dXk[k] + da4dXk[k])
                      / denom; // M.1.4
  }
  this->dZdT = - (Z3 * da1dT + Z2 * da2dT + Z * da3dT + da4dT) / denom; // M.1.1
  this->dZdp = - (Z3 * da1dp + Z2 * da2dp + Z * da3dp + da4dp) / denom; // M.1.2
  // other phasic derivatives (2.2.86, 2.2.87)
  denom = (m_MW + this->delta1 * this->Bm * m_rho) * (m_MW + this->delta2 * this->Bm * m_rho);
  double dpdT = m_rho * R_u / (m_MW - this->Bm * m_rho) - m_rho * m_rho * this->dAmdT / denom;
  double dpdrho = m_MW * R_u * this->T / pow(m_MW - this->Bm * m_rho, 2)
                  - this->Am * m_rho * m_MW * (2*m_MW + this->delta1_p_delta2 * this->Bm * m_rho)/ denom / denom;
  // if (phase == this->iSINGLEPH_TP) throw(-1);
  double drhodXk [this->num_species];
  LOOP_k_N (this->num_species) {
    // K.0.36
    drhodXk[k] = this->P / R_u / this->T * (this->MW[k] * Z - m_MW * this->dZdXk[k]) / Z / Z; REPLACEIFWEIRD(drhodXk[k]);
  }

  // (m.2.1)
  double zdelta1 = Z + this->delta1 * Bm_tilde; double zdelta2 = Z + this->delta2 * Bm_tilde;
  double dAmDiff = this->T * this->dAmdT - this->Am;
  double Rtilde = 1/(this->delta1-this->delta2) / this->Bm;
  double Stilde = dAmDiff;
  double Ttilde = log(zdelta1/zdelta2);
  // double dDepEdT = - dAmDiff * Bm_tilde / this->Bm * (this->dZdT + Z * Bm_tilde/this->T) / zdelta1/ zdelta2
  //                  + this->T * d2AmdT2 * log(zdelta1/zdelta2) / (this->delta1-this->delta2) / this->Bm; REPLACEIFWEIRD(dDepEdT);
  // double dDepEdT = - dAmDiff * Bm_tilde / this->Bm * (this->dZdT + Z /this->T) / zdelta1/ zdelta2
  //                  + this->T * d2AmdT2 * log(zdelta1/zdelta2) / (this->delta1_m_delta2) / this->Bm; REPLACEIFWEIRD(dDepEdT);
  double alpha = 1. / zdelta1 / zdelta2 * (this->delta2 - this->delta1) * Bm_tilde * (this->dZdT + Z / this->T);
  double dDepEdT = 1.0 / (this->delta1_m_delta2 * this->Bm) *
                   (log(zdelta1/zdelta2) * this->T * this->d2AmdT2
                    + dAmDiff * alpha);
  REPLACEIFWEIRD(dDepEdT);
  double dDepEdp = dAmDiff * Bm_tilde * (this->delta1-this->delta2) *(Z / this->P - this->dZdp) / (this->delta1-this->delta2) / this->Bm / zdelta1 / zdelta2; REPLACEIFWEIRD(dDepEdp);

  double dDepEdrho = dAmDiff * m_MW / (this->delta1 * this->Bm * m_rho + m_MW) / (this->delta2 * this->Bm * m_rho + m_MW); REPLACEIFWEIRD(dDepEdrho);
  double dDepEdXk [this->num_species];
  LOOP_k_N (this->num_species) {
    double dRtildedXk = - 1 / this->delta1_m_delta2 * this->dBmdXk[k] / this->Bm / this->Bm; // (M.2.2)
    double dStildedXk = this->T * this->d2AmdTdXk[k] - this->dAmdXk[k];
    double Bmtildediff = Bm_tilde * this->dZdXk[k] - Z * dBmtildedXk[k];
    double dTildedXk = 1.0 * Bmtildediff * (this->delta2 - this->delta1) / (zdelta1 * zdelta2);
    double term1 = Stilde * Ttilde * dRtildedXk;
    double term2 = Rtilde * Ttilde * dStildedXk;
    double term3 = Stilde * Rtilde * dTildedXk;
    dDepEdXk[k] = term1 + term2 + term3; //
    REPLACEIFWEIRD(dDepEdXk[k]);
  }

  // save derivative quantities that is needed to obtain phasic quantities
  if (phase == this->iLIQPH_TP) {
    this->dpdT_L = dpdT;
    this->dpdrho_L = dpdrho;
    this->dZdT_L = this->dZdT;
    LOOP_k_N (this->num_species) {
      this->drhoLdx[k] = drhodXk[k];
      this->dZLdx[k] = this->dZdXk[k];
    }
    this->dDepEdT_L = dDepEdT;
    this->dDepEdp_L = dDepEdp;
    LOOP_k_N (this->num_species) {
      this->dDepEdXk_L[k] = dDepEdXk[k];
    }
    // throw(-1);
  } else if (phase == this->iVAPPH_TP) {
    this->dpdT_V = dpdT;
    this->dpdrho_V = dpdrho;
    this->dZdT_V = this->dZdT;
    LOOP_k_N (this->num_species) {
      this->drhoVdy[k] = drhodXk[k];
      this->dZVdy[k] = this->dZdXk[k];
    }
    this->dDepEdT_V = dDepEdT;
    this->dDepEdp_V = dDepEdp;
    LOOP_k_N (this->num_species) {
      this->dDepEdXk_V[k] = dDepEdXk[k];
    }
    // throw(-1);
  } else {
    this->dPdT = dpdT;
    this->dpdrho = dpdrho;
    LOOP_k_N (this->num_species) this->drhodchi[k] = drhodXk[k];
    // don't do anyting here to avoid replacement of the two-phase derivatives
    // this->dDepEdT = dDepEdT;
    // this->dDepEdp = dDepEdp;
    // LOOP_k_N (this->num_species) {
    //   this->dDepEdXk[k] = dDepEdXk[k];
    // }
  }
} // end of PengRobinson::computeBasicDerivatives
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @brief: compute fugacity coefficients and fugacities
// ----------------------------------------------------------------------------
void PengRobinson::computeFugacityCoefficient(double * phi, double * f, int phase) {
  double lnPhi [this->num_species];

  double Z;
  if (phase == this->iLIQPH_TP) {
    Z = this->Z_L;
  } else if (phase == this->iVAPPH_TP) {
    Z = this->Z_V;
  } else {
    Z = this->P / (this->rho * this->gas_constant * this->T);
  }
  double R_u = this->gas_constant_universal;
  double Bm_tilde = this->Bm * this->P / R_u / this->T;

  LOOP_k_N (this->num_species) {
    int i = k;
    double sum = 0.;
    LOOP_l_N (this->num_species) {
      int j = l;
      int apos = j * this->num_species + i;
      // Table 12
      double Aji = this->cst_a[apos] * std::pow(1.0 + this->cst_c[apos] * (1.0 - std::sqrt(this->T / this->Tcrit_IJ[apos])), 2);
      sum += this->X[j] * Aji;
    }
    // (M.3.18)
    lnPhi[k] = - log(Z - Bm_tilde)
               + this->cst_b[k] / this->Bm * (Z-1.0)
               - log((this->delta1 * Bm_tilde +  Z)/(this->delta2 * Bm_tilde + Z))
               * this->Am / this->Bm / R_u / this->T / (this->delta1 - this->delta2)
               * (2. * sum / this->Am - this->cst_b[k] / this->Bm);
    // printf("lnphi[%d] = %.8e\n", k, lnPhi[k]);
    if (Z > 0.0) {
      // do check
      if (std::isnan(lnPhi[k]) || std::isinf(lnPhi[k])) {
        std::cerr << "Error: PengRobinson::computeFugacityCoefficient: "
                  << "lnPhi[" << k << "] is NaN or Inf." << std::endl;
        printf("k = %d, Z = %.4f, Bm_tilde = %.4f, MW = %.4f, T = %.4f, P = %.4f\n",
               k, Z, Bm_tilde, this->MW_M, this->T, this->P);
        throw (-1);
      }
    }
    phi[k] = exp(lnPhi[k]);
    f[k] = phi[k] * this->P * this->X[k]; // (K.0.1)
  }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @brief: full VLE computation
// ----------------------------------------------------------------------------
int PengRobinson::VLE(const double T_in, const double P_in, const double* X_in,
    double& beta_out, double* x_out, double* y_out) {
  this->T = T_in; this->P = P_in;
  this->SetMolarFraction(X_in);
  this->SetMolecularWeightMixture();
  // Solve single-phase state using Gibbs energy minimization
  this->rho = this->GetDensityFromPressureTemperature_TP(P_in, T_in, this->iSINGLEPH_TP);
  this->computeBasicDerivatives(this->iSINGLEPH_TP); // compute derivatives of the singlephase state
  double Z = P_in / (this->rho * this->gas_constant * T_in);
  // compute fugacity coefficients and fugacities
  this->computeFugacityCoefficient(this->phi_fugacity, this->f_fugacity, this->iSINGLEPH_TP);
  // Estimate K_i using (3.4.1)
  double d [this->num_species]; double K [this->num_species];
  LOOP_k_N (this->num_species) { 
    d[k] = log(this->X[k]) + log(this->phi_fugacity[k]);
    K[k] = this->Pcrit[k] / P_in * exp(5.37 * (1 + omega[k]) * (1. - this->Tcrit[k] / T_in));
  }
  int flag;
  double beta;
  double x [this->num_species], y [this->num_species];
  // RUN SSI for 3 steps: htis will update beta, x, y
  this->SSI(T_in, P_in, X_in, K, beta, x, y, flag, 1e-7, 3);
  if (beta > 0.0 && beta < 1.0) {
    // make sure that phi_L is larger than phi_fugacity
    double sum_phi_L = 0.0; double sum_phi = 0.0; double sum_phi_V = 0.0;
    LOOP_k_N (this->num_species) {
      sum_phi_L = this->X[k] * this->phi_L[k];
      sum_phi = this->X[k] * this->phi_fugacity[k];
      sum_phi_V = this->X[k] * this->phi_V[k];
    }
  }
  // RUN StabilityAnalysis
  int nphase;
  // update K from this stability analysis
  this->StabilityAnalysis(T_in, P_in, X_in, K, d, nphase, 1e-3);
  if (flag == -1 && nphase == 1) {
    // output the SSI because SSI and stability analysis agree that system is single phase
  } else {
    // SSI and stability analysis disagree.
    // Run SSI again with larger nmax and use K estimation from stability analysis
    this->SSI(T_in, P_in, X_in, K, beta, x, y, flag, 1e-7, 1000);
  }
  // output beta, x, y
  beta_out = beta;
  LOOP_k_N (this->num_species) {
    x_out[k] = x[k]; y_out[k] = y[k];
  }
  //evaluate nphase from betea and flag
  this->flag_vle = flag;
  nphase = (flag == -1) ? 1 : 2;
  // printf("Final results:\n");
  // printf("beta_out = %.4f\n", beta_out);
  // printf("nphase = %d\n", nphase);
  // std::cout << "End of VLE computation" << std::endl;
  return nphase;
} // end of PengRobinson::VLE
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @brief: SSI computation: get phase equilibrium
// ----------------------------------------------------------------------------
void PengRobinson::SSI(const double T_in, const double P_in, const double* X_in,
    const double* K_in, double& beta_out, double* x_out, double* y_out, int& flag,
    const double eps, const int max_iter) {
  // check Algorithm 8
  // assume that the temperature and pressure are already set
  double t = T_in; double p = P_in; // zzz...
  double K [this->num_species]; // K values
  LOOP_k_N (this->num_species) K[k] = K_in[k];
  double err = 1e10; // initial error
  int iter = 0; // iteration counter
  while (err > eps) { // Main Loop
    iter++;
    this->niter_ssi = iter; // update this
    if (iter == max_iter) {
      return;
    }
    // Check if single phase
    double sum_XK = 0.;
    double sum_XoK = 0.;

    LOOP_k_N (this->num_species) {
      sum_XK += X_in[k] * K[k];
      sum_XoK += X_in[k] / K[k];
    }

    if ((sum_XK - 1. < 0.0) && (std::fabs(sum_XK - 1.) > 1e-15)) {
      // if ((sum_XK - 1. <= 1e-15)) {
      // if ((sum_XK - 1. < 0.0)) {
      // this mixture is pure liquid
      beta_out = 0.0;
      LOOP_k_N (this->num_species) {
        x_out[k] = X_in[k];
        y_out[k] = 0.0; // no vapor phase
      }
      flag = -1;
      return; // exit
              // } else if (1. - sum_XoK >= -1e-15) {
    } else if ((1. - sum_XoK > 1e-15) && (true)) {
      // hits mixture is pure vapor
      beta_out = 1.0;
      LOOP_k_N (this->num_species) {
        y_out[k] = X_in[k];
        x_out[k] = 0.0; // no liquid phase
      }
      flag = -1;
      return; // exit
    } else {
      flag = 1;
      // this mixture has VLE
      // compute beta from (2.2.37, Rachford-Rice)
      // Returns tuple: F(beta), F'(beta), F''(beta)
      auto rr = [&](double beta) {
        double F   = 0.0;
        double dF  = 0.0;
        double d2F = 0.0;

        for (size_t i = 0; i < this->num_species; ++i) {
          double a = K[i] - 1.0;
          double denom = 1.0 + beta * a;
          // safeguard against division by zero
          if (denom <= 0.0) denom = std::numeric_limits<double>::min();

          double inv = 1.0 / denom;
          double zi  = X_in[i];

          F   += zi * a * inv;
          dF  += -zi * a * a * inv * inv;
          d2F +=  2.0 * zi * a * a * a * inv * inv * inv;
        }
        return std::make_tuple(F, dF, d2F);
      };
      // Max iterations
      const boost::uintmax_t maxit = 500;
      boost::uintmax_t it = maxit;

      // std::pair<double, double> sol = boost::math::tools::toms748_solve(fun, 0., 1., tol, it);
      // beta_out = 0.5 * (sol.first + sol.second);
      // Initial guess (clamped)
      double num = 0.0, den = 0.0;
      for (size_t i = 0; i < this->num_species; ++i) {
        double Ki = std::max(K[i], 1e-30); // guard against div by zero
        num += X_in[i] * (1.0 - 1.0 / Ki);
        den += X_in[i] * ((K[i] - 1.0) / Ki);
      }
      double beta0 = (den != 0.0) ? (num / den) : 0.5;
      beta0 = std::max(0.0, std::min(1.0, beta0)); // clamp to [0, 1]

      // Bracket for halley_iterate
      double min_beta = 0.0;
      double max_beta = 1.0;
      int digits = std::numeric_limits<double>::digits; // full double precision
      double sol = boost::math::tools::halley_iterate(rr, beta0, min_beta, max_beta, digits, it);

      // if (it == maxit) {
      //   // if we reach the maximum number of iterations, we have to return
      //   std::cerr << "Error: PengRobinson::SSI: "
      //             << "Halley's method did not converge after " << maxit << " iterations." << std::endl;
      //   throw (-1);
      // }

      beta_out = sol;
      SAFEGUARD(beta_out, "beta_out");

      // compute x (3.4.2), compute y (3.4.3)
      LOOP_k_N (this->num_species) {
        x_out[k] = X_in[k] / (1. + beta_out * (K[k] - 1.));
        y_out[k] = K[k] * x_out[k];
      }

      // obtain liquid mixture to update fugacity coefficients
      this->generateVLstructure(t, p, x_out, iLIQPH_TP);
      // obtain vapor mixture to update fugacity coefficients
      this->generateVLstructure(t, p, y_out, iVAPPH_TP);

      err = 0.0;
      LOOP_k_N (this->num_species) {
        // compute error
        double err_k = log(K[k]) - log(this->phi_L[k]) + log(this->phi_V[k]);
        // err = std::max(err, abs(err_k));
        err += err_k * err_k; // use squared error
      }
      err = sqrt(err / this->num_species); // normalize error
                                           // update K for next iteration
      LOOP_k_N (this->num_species) {
        K[k] = this->phi_L[k] / this->phi_V[k];
      }
    }
    // printf("SSI iteration %d, err = %.16e\n", iter, err);
    } // end of main loop
    } // end of PengRobinson::SSI
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: run stability analysis based on input K
// @param T: temperature
// @param P: pressure
// @param X: molar fractions of the mixture
// @param K: K values for each species
// @param d: ?
// @param_out nphases: resulting nphases
// @param eps: convergence criterion
// ----------------------------------------------------------------------------
void PengRobinson::StabilityAnalysis(const double T_in, const double P_in, const double* X_in,
    double* K, const double* d_in, int &nphases, const double eps) {
  this->trivial_vle = false; // reset trivial VLE flag
  // check Algorithm 9
  int n = 0;
  double K_0 [this->num_species];
  double chi_V [this->num_species];
  double chi_L [this->num_species];
  double sum_chi_V = 0; double sum_chi_L = 0;
  LOOP_k_N (this->num_species) { 
    K_0[k] = K[k];
    chi_V[k] = K_0[k] * X_in[k];
    chi_L[k] = X_in[k] / K_0[k];
    sum_chi_V += chi_V[k]; sum_chi_L += chi_L[k];
  }
  // normalize chi_V and chi_L
  LOOP_k_N (this->num_species) {
    chi_V[k] /= sum_chi_V;
    chi_L[k] /= sum_chi_L;
  }

  double tpd_L, tpd_V, nphases_L, nphases_V;
  // update chi_V, tpd_V, and nphases_V
  nphases_V = this->UpdateTrialComposition(T_in, P_in, X_in, d_in, iVAPPH_TP, chi_V,
      tpd_V, eps);
  // update chi_L, tpd_L, and nphases_L
  nphases_L = this->UpdateTrialComposition(T_in, P_in, X_in, d_in, iLIQPH_TP, chi_L,
      tpd_L, eps);

  // Now make the final decision based on the values of TPD functions
  if ((nphases_V > 1) && (nphases_L > 1)) {
    // Both trials unstable, pick on wiht lowest TPD function value
    if (tpd_V < tpd_L) {
      // gas-like mixture
      LOOP_k_N (this->num_species) K[k] = chi_V[k] / X_in[k];
    } else {
      // liquid-like mixture
      LOOP_k_N (this->num_species) K[k] = X_in[k] / chi_L[k];
    }
    nphases = 2;
  } else if (nphases_V > 1) {
    // gas-like mixture
    LOOP_k_N (this->num_species) K[k] = chi_V[k] / X_in[k];
    nphases = 2;
  } else if (nphases_L > 1) {
    // liquid-like mixture
    LOOP_k_N (this->num_species) K[k] = X_in[k] / chi_L[k];
    nphases = 2;
  } else {
    // this is a trivial solution; mixture stable as single phase
    trivial_vle = true;
    LOOP_k_N (this->num_species) K[k] = 1.0;
    nphases = 1;
  }
} // end of PengRobinson::StabilityAnalysis
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: update trial composition, to be used in StabilityAnalysis
// @param T: temperature
// @param P: pressure
// @param X: molar fractions of the mixture
// @param d: ?
// @param phase
// @param chi_in: trial composition
// @param_out tpd_phase: TPD function value for the phase
// @param_out chi_out: resulting composition
// @param eps_TPD: convergence criterion
// @param_out nphases: number of phases
// ----------------------------------------------------------------------------
int PengRobinson::UpdateTrialComposition(const double T_in, const double P_in, const double* X_in,
    const double* d, int phase, double* chi, 
    double& tpd_phase,
    const double eps_tpd) {
  // check Algorithm 10
  int n = 0;
  double chi_new [this->num_species];
  // main loop
  double m_phi [this->num_species];
  double err = 1e100; // initial error
  while (err > eps_tpd) {
    n++;
    // if (n > 1000) break;
    err = 0.0;
    // Obtain mixture at phase phase
    this->generateVLstructure(T_in, P_in, chi, phase);
    if (phase == this->iLIQPH_TP) {
      LOOP_k_N (this->num_species) m_phi[k] = this->phi_L[k];
    } else if (phase == this->iVAPPH_TP) {
      LOOP_k_N (this->num_species) m_phi[k] = this->phi_V[k];
    }
    double sum_chi_new = 0.0;
    LOOP_k_N (this->num_species) {
      chi_new[k] = exp(d[k] - log(m_phi[k]));
      sum_chi_new += chi_new[k];
    }

    LOOP_k_N (this->num_species) {
      chi_new[k] /= sum_chi_new; // normalize
      err += (chi_new[k] - chi[k]) * (chi_new[k] - chi[k]);
    }
    err = sqrt(err);
    // LOOP_k_N (this->num_species) {
    //   if (chi_new[k] < 0. || chi_new[k] > 1.0) {
    //     printf("d[%d] = %.4e, logPhi[%d] = %.4e, chi_new[%d] = %.4e, chi[%d] = %.4e\n",
    //         k, d[k], k, log(m_phi[k]), k, chi_new[k], k, chi[k]);
    //
    //     std::cout << "Weird chi_new values in UpdateTrialComposition: "
    //       << "chi_new[" << k << "] = " << chi_new[k] << std::endl;
    //     exit(-1);
    //   }
    // }

    // update
    LOOP_k_N (this->num_species) {
      chi[k] = chi_new[k];
    }
  }
  // Now compute the TPD function (3.4.4)
  tpd_phase = 1.0;
  assert(phase == this->iLIQPH_TP || phase == this->iVAPPH_TP);
  this->generateVLstructure(T_in, P_in, chi_new, phase);
  LOOP_k_N (this->num_species) {
    tpd_phase += chi[k] * (log(chi[k]) + log(m_phi[k]) - d[k] - 1.0);
  }
  // Determine if this is a trivial soluiton
  bool trivial = false;
  double err_trivial = 0;
  LOOP_k_N (this->num_species) err_trivial += (chi_new[k] - X_in[k]) * (chi_new[k] - X_in[k]);
  if (err_trivial < 1e-5) {
    trivial = true;
  } 

  int nphases;
  if ((tpd_phase < 0.) && !trivial) {
    // VLE exists
    nphases = 2;
  } else {
    // this is a trivial solution
    nphases = 1;
  }
  return nphases;
} // end of PengRobinson::UpdateTrialComposition
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: solve ERHO problem, will be used in determining T,P from conservatives
// @param in e_in: specific internal energy
// @param in rho_in: density
// @param in X_in: mixture molar fractions
// @param out T_out: resulting temperature
// @param out P_out: resulting pressure
// @param out beta_out: resulting vapor volume fraction
// @param out x_out: resulting liquid molar fractions
// @param out y_out: resulting vapor molar fractions
// ----------------------------------------------------------------------------
void PengRobinson::solveERHO(const double e_in, const double rho_in, const double* X_in,
    double& T_out, double& P_out, double& beta_out, double* x_out, double* y_out,
    const double eps_ERHO, const double deltaT_max, const double deltaP_max,
    const double Nmax,
    const double T_0, const double P_0, const double lambda_0, const bool accelerate) {
  int n = 0;
  // e_in is per unit mass, [J/kg]
  // Initialize temperature and pressure guesses;
  double T_guess = T_0; // initial guess
  double P_guess = P_0; // initial guess
  double err = 1e10; // initial error

  // Acceleration parameters
  bool step_accelerate = false; // always start with no acceleration
  int n_acc = 0; // number of consecutive accelerated so far
  const bool N_acc = 5; // number of consecutive iterations to accelerate
  const double eps_acc = 1e-7; // error threshold for acceleration
                               //
  // Main loop
  while (err > eps_ERHO && n < Nmax) {
    n++; // increment iteration counter
    if (!step_accelerate) {
      n_acc = 0; // reset acceleration counter
      // Solve the Tp problem to get beta, x, y (updated beta, x_TP, y_TP, E, rho)
      this->solveTP(T_guess, P_guess, X_in);
    }
    // Get all mixture derivatives
    this->computeAdditionalDerivatives(T_guess, P_guess, X_in, this->betaV, this->x_TP, this->y_TP);
    this->SetMolarFraction(X_in); // make sure we have the correct MW_M...
    double F1 = (this->E_mass - e_in) * this->MW_M; // remember that these are in J/kg, so we convert it to J/mol
    double F2 = this->rho - rho_in;
    double dF1dT = this->dEdT; // J/kmol/K
    double dF1dp = this->dEdp; // J/kmol/Pa
    double dF2dT = this->drhodT;
    double dF2dp = this->drhodp;
    // Assemble Jacobian and compute variations, (3.4.11) and the inverse
    Eigen::MatrixXd scriptJ(2, 2); // Jacobian matrix
    Eigen::VectorXd F(2); // Objective function vector
    scriptJ(0, 0) = dF1dT; scriptJ(0, 1) = dF1dp;
    scriptJ(1, 0) = dF2dT; scriptJ(1, 1) = dF2dp;
    F(0) = F1; F(1) = F2;
    Eigen::MatrixXd scriptJ_inv = scriptJ.inverse(); // Inverse of Jacobian
    Eigen::VectorXd delta = -scriptJ_inv * F; // (3.4.11), look at line 19 of the Algorithm 13
    double DeltaT = delta(0); // temperature change
    double DeltaP = delta(1); // pressure change
    // Declare error
    err = F1*F1 + F2*F2; // or, fabs(DeltaT) + fabs(DeltaP);
    // err = fabs(DeltaT) + fabs(DeltaP);
    // Update Unknowns
    double T_new = T_guess + DeltaT;
    double P_new = P_guess + DeltaP;
    // Optional: temperature/pressure update can be done using line-search method
    double lambda = std::min(lambda_0, std::min(fabs(deltaT_max/DeltaT), fabs(deltaP_max/DeltaP)));
    T_new = T_guess + lambda * DeltaT;
    P_new = P_guess + lambda * DeltaP;
    // update guesses
    T_guess = T_new;
    P_guess = P_new;
    // TODO: evaluate if we accelerate the step
    step_accelerate = (n_acc <= N_acc) && (err < eps_acc);
    step_accelerate = step_accelerate && accelerate; // only accelerate if global accelerate is true
    if (step_accelerate) {
      n_acc++; // increment acceleration counter
      // Approximate beta, x, y using VLE derivatives obtained from previous iterations
      double sum_x = 0.0, sum_y = 0.0;
      LOOP_k_N (this->num_species) {
        this->x_TP[k] += this->dxdT[k] * DeltaT + this->dxdp[k] * DeltaP;
        this->y_TP[k] += this->dydT[k] * DeltaT + this->dydp[k] * DeltaP;
        sum_x += this->x_TP[k]; sum_y += this->y_TP[k];
      }
      // Normalize x and y
      LOOP_k_N (this->num_species) {
        this->x_TP[k] /= sum_x; // normalize liquid molar fractions
        this->y_TP[k] /= sum_y; // normalize vapor molar fractions
        REPLACEIFWEIRD(this->x_TP[k]);
        REPLACEIFWEIRD(this->y_TP[k]);
      }
      double K [this->num_species]; // K values
                                    // Compute K values
      LOOP_k_N (this->num_species) {
        K[k] = this->y_TP[k] / this->x_TP[k];
        REPLACEIFWEIRD(K[k]);
      }
      // compute beta using (2.2.37, Rachford-Rice)
      this->betaV = this->solveRachfordRice(X_in, this->x_TP, this->y_TP);
    } // end of if step_accelerate
  } // end of solveERHO main loop
  T_out = T_guess; // resulting temperature
  P_out = P_guess; // resulting pressure
  beta_out = this->betaV;
  LOOP_k_N (this->num_species) {
    x_out[k] = this->x_TP[k]; // resulting liquid molar fractions
    y_out[k] = this->y_TP[k]; // resulting vapor molar fractions
  }
  this->SyncZPhi();
  // printf("takes %d iterations to solve ERHO, err = %.4e\n", n, err);
} // end of PengRobinson::solveERHO
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: solve ERHO problem, will be used in determining T,P from conservatives, using mass fractions
// @param in e_in: specific internal energy
// @param in rho_in: density
// @param in Y_in: mixture mass fractions
// @param out T_out: resulting temperature
// @param out P_out: resulting pressure
// @param out beta_out: resulting vapor volume fraction
// @param out x_out: resulting liquid molar fractions
// @param out y_out: resulting vapor molar fractions
// ----------------------------------------------------------------------------
void PengRobinson::solveERHOY(const double e_in, const double rho_in, const double* Y_in,
    double& T_out, double& P_out, double& beta_out, double* x_out, double* y_out,
    const double eps_ERHO, const double deltaT_max, const double deltaP_max,
    const double Nmax,
    const double T_0, const double P_0, const double lambda_0, const bool accelerate) {
  this->SetMassFractionFromY(Y_in);
  this->SetMolecularWeightMixture();
  this->SetMolarFractionFromY(); // this updates X
  double X_in [this->num_species];
  LOOP_k_N(this->num_species) X_in[k] = this->X[k]; // copy X to X_in
  // Now call solveERHO
  this->solveERHO(e_in, rho_in, X_in, T_out, P_out, beta_out, x_out, y_out,
      eps_ERHO, deltaT_max, deltaP_max, Nmax, T_0, P_0, lambda_0, accelerate);
} // end of PengRobinson::solveERHOY
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: solve PRHO problem, will be used in Double Flux model as we update from pressure
// @param in P_in: pressure
// @param in rho_in: density
// @param in X_in: mixture molar fractions
// @param out T_out: resulting temperature
// @param out beta_out: resulting vapor volume fraction
// @param out x_out: resulting liquid molar fractions
// @param out y_out: resulting vapor molar fractions
// ----------------------------------------------------------------------------
void PengRobinson::solvePRHO(const double P_in, const double rho_in, const double* X_in,
    double& T_out, double& beta_out, double* x_out, double* y_out,
    const double eps_PRHO, const double deltaT_max, const double N_max,
    const double T_0, const double lambda_0) {
  int n = 0;
  // Initialize temperature guess and known value of density
  double T_guess = T_0; // initial guess
  this->rho = rho_in;
  this->P = P_in;
  double err = 1e10;
  double T_new;
  // main loop
  while ((err > eps_PRHO) && (n <= N_max)) {
    n++;
    // Solve TP problem to get beta, x, y (updated beta, x_TP, y_TP)
    this->solveTP(T_guess, P_in, X_in);
    // Get all mixture derivatives
    this->computeAdditionalDerivatives(T_guess, P_in, X_in, this->betaV, this->x_TP, this->y_TP);
    // printf("Iteration %d, T_guess = %.4f, sos = %.4f, rho = %.4f, P = %.4f\n",
    //         n, T_guess, this->sos, this->rho, this->P);
    // Update objective function and its derivatives
    double F = this->rho - rho_in;
    double dFdT = this->drhodT; // computed in computeAdditionalDerivatives
    err = F;
    // Update unknown; T_new from (3.4.10) and DeltaT = T_new - T_guess
    T_new = T_guess - F / dFdT; // (3.4.10)
    double dT = T_new - T_guess; // temperature change
    // Optional: temperature update using line-search method
    double lambda = std::min(lambda_0, fabs(deltaT_max / dT));
    // Update temperature using modified version of (3.4.10)
    T_new = T_guess - lambda * F / dFdT;
    // printf("Iteration %d, T_guess = %.4f, T_new = %.4f, F = %.4e, dFdT = %.4e, err = %.4e, lambda= %.4e, dT = %.4e\n",
    //        n, T_guess, T_new, F, dFdT, err, lambda, dT);
    T_guess = T_new; // update temperature guess
    if (n == N_max) {
      printf("Warning: PengRobinson::solvePRHO: reached maximum number of iterations %d, err = %.4e\n",
          N_max, err);
      printf("Last T_guess = %.4f, P_in = %.4f, rho_in = %.4f\n", T_guess, P_in, rho_in);
    }
  } // end of main loop
  T_out = T_guess; // resulting temperature
  beta_out = this->betaV;
  LOOP_k_N (this->num_species) {
    x_out[k] = this->x_TP[k]; // resulting liquid molar fractions
    y_out[k] = this->y_TP[k]; // resulting vapor molar fractions
  }
  // update sos, cp, Z
  this->computeAdditionalDerivatives(T_guess, P_in, X_in, this->betaV, this->x_TP, this->y_TP);
  this->SyncZPhi();
} //end of PengRobinson::solvePRHO
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: solve PRHO problem, will be used in Double Flux model as we update from pressure
// @param in P_in: pressure
// @param in rho_in: density
// @param in Y_in: mixture mass fractions
// @param out T_out: resulting temperature
// @param out beta_out: resulting vapor volume fraction
// @param out x_out: resulting liquid molar fractions
// @param out y_out: resulting vapor molar fractions
// ----------------------------------------------------------------------------
void PengRobinson::solvePRHOY(const double P_in, const double rho_in, const double* Y_in,
    double& T_out, double& beta_out, double* x_out, double* y_out,
    const double eps_PRHO, const double deltaT_max, const double N_max,
    const double T_0, const double lambda_0) {
  this->SetMassFractionFromY(Y_in);
  this->SetMolecularWeightMixture();
  this->SetMolarFractionFromY(); // this updates X
  double X_in [this->num_species];
  LOOP_k_N(this->num_species) X_in[k] = this->X[k]; // copy X to X_in
  // Now call solvePRHO
  this->solvePRHO(P_in, rho_in, X_in, T_out, beta_out, x_out, y_out,
      eps_PRHO, deltaT_max, N_max, T_0, lambda_0);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// @\brief: solve VLE and set energies
// @param in T_in: temperature
// @param in P_in: pressure
// @param in X_in: mixture molar fractions
// ----------------------------------------------------------------------------
void PengRobinson::solveTP(const double T_in, const double P_in, const double* X_in, const bool calc_derivatives) {
  // Solve VLE problem
  int nphase = this->VLE(T_in, P_in, X_in, this->betaV, this->x_TP, this->y_TP);
  // Obtain liquid mixture
  this->generateVLstructure(T_in, P_in, this->x_TP, iLIQPH_TP);
  // Obtain vapor mixture
  this->generateVLstructure(T_in, P_in, this->y_TP, iVAPPH_TP);
  // Compute mixture density (2.2.34)
  double M_over_rho;
  if (nphase > 1) {
    // both MW_L and MW_V are defined
    M_over_rho = this->betaV * this->MW_V / this->rho_V +
      (1.0 - this->betaV) * this->MW_L / this->rho_L;
  } else {
    if (this->betaV > 0.5) { // pure vapor phase
      M_over_rho = this->betaV * this->MW_V / this->rho_V;
    } else { // pure liquid phase
      M_over_rho = (1.0 - this->betaV) * this->MW_L / this->rho_L;
    }
  }

  if (M_over_rho <= 0.0 || std::isnan(M_over_rho) || std::isinf(M_over_rho)) {
    printf("rho_L = %.4f, rho_V = %.4f, MW_L = %.4f, MW_V = %.4f, betaV = %.4f\n",
           this->rho_L, this->rho_V, this->MW_L, this->MW_V, this->betaV);
    std::cerr << "Error: PengRobinson::solveTP: M_over_rho is non-positive: "
              << M_over_rho << std::endl;
    // assert(false);
  }

  // recompute M with X_in to make sure this is correct
  this->SetMolarFraction(X_in); // update MW_M
  this->gas_constant = this->gas_constant_universal / this->MW_M; // update gas constant
  this->rho = this->MW_M / M_over_rho; // mixture density
  this->Am = this->GetAm(T_in); this->Bm = this->GetBm();
  this->computeBasicDerivatives(this->iSINGLEPH_TP); // compute derivatives of the singlephase state
  // Compute mixture internal energy (2.2.35) //
  double E_volume;
  if (nphase > 1) {
    E_volume = this->betaV * this->E_V + (1.0 - this->betaV) * this->E_L;
  } else {
    if (this->betaV > 0.5) { // pure vapor phase
      E_volume = this->betaV * this->E_V;
    } else { // pure liquid phase
      E_volume = (1.0 - this->betaV) * this->E_L;
    }
  }

  this->E = E_volume; // J/mol
  this->E_mass = E_volume / this->MW_M; // convert to per kg
                                        //
  // Compute mixture additional energies
  double h = this->E_mass + P_in / this->rho; // J/kg
  // TODO: entropy and gibbs?
  if (calc_derivatives)
    this->computeAdditionalDerivatives(T_in, P_in, X_in, this->betaV, this->x_TP, this->y_TP);
  this->SyncZPhi(); // sync Z_phi
} // end of PengRobinson::solveTP
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: solve VLE and set energies
// @param in T_in: temperature
// @param in P_in: pressure
// @param in Y_in: mixture mass fractions
// ----------------------------------------------------------------------------
void PengRobinson::solveTPY(const double T_in, const double P_in, const double* Y_in, const bool calc_derivatives) {
  this->SetMassFractionFromY(Y_in);
  this->SetMolecularWeightMixture();
  this->SetMolarFractionFromY(); // this updates X
  double X_in [this->num_species];
  LOOP_k_N(this->num_species) X_in[k] = this->X[k]; // copy X to X_in
  this->solveTP(T_in, P_in, X_in, calc_derivatives);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: compute additional derivatives for ERHO or PRHO problems
// @param in T_in: temperature
// @param in P_in: pressure
// @param in X_in: mixture molar fractions
// @param in beta_in: vapor volume fraction
// @param in x_in: liquid molar fractions
// @param in y_in: vapor molar fractions
void PengRobinson::computeAdditionalDerivatives(const double T_in, const double P_in, const double* X_in,
      const double beta_in, const double* x_in, const double* y_in) {
    // Algorithm 4
  // Get basic derivatives
  this->computeTPvariations(T_in, P_in, X_in, beta_in, x_in, y_in);
  // phase derivatives
  this->computeAdditionalPhaseDerivatives(T_in, P_in, X_in, this->iLIQPH_TP, x_in);
  this->computeAdditionalPhaseDerivatives(T_in, P_in, X_in, this->iVAPPH_TP, y_in);
  // compute dEdT from (K.0.42)
  this->dEdT = this->dBetadT * (this->E_V - this->E_L)
              + this->betaV * this->dEdT_V + (1.0 - this->betaV) * this->dEdT_L; // J/mol/K
  // compute dEdp from (K.0.43)
  this->dEdp = this->dBetadp * (this->E_V - this->E_L)
              + this->betaV * this->dEdp_V + (1.0 - this->betaV) * this->dEdp_L;
  // compute drhodT from (2.2.49)
  double M_over_rho_V = (this->betaV != 0.0) ? this->MW_V / this->rho_V : 0.0;
  double M_over_rho_L = (this->betaV != 1.0) ? this->MW_L / this->rho_L : 0.0;
  double parant = this->dBetadT * (M_over_rho_V - M_over_rho_L);

  // evaluate single-phased ness to avoid NaN or Inf
  bool add_vapor = (!std::isnan(this->rho_V) && !std::isinf(this->rho_V) && this->rho_V > 0.0 &&
                    !std::isnan(this->MW_V) && !std::isinf(this->MW_V) && this->MW_V > 0.0);
  bool add_liquid = (!std::isnan(this->rho_L) && !std::isinf(this->rho_L) && this->rho_L > 0.0 &&
      !std::isnan(this->MW_L) && !std::isinf(this->MW_L) && this->MW_L > 0.0);

  // prev: 09242025
  // if (this->betaV != 0.0) parant += this->betaV / this->rho_V / this->rho_V * (this->rho_V * this->dMWVdT - this->MW_V * this->drhodT_V);
  // if (this->betaV != 1.0) parant += (1. - this->betaV) / this->rho_L / this->rho_L * (this->rho_L * this->dMWLdT - this->MW_L * this->drhodT_L);
  if (add_vapor) parant += this->betaV / this->rho_V / this->rho_V * (this->rho_V * this->dMWVdT - this->MW_V * this->drhodT_V);
  if (add_liquid) parant += (1. - this->betaV) / this->rho_L / this->rho_L * (this->rho_L * this->dMWLdT - this->MW_L * this->drhodT_L);
  this->drhodT = - this->rho * this->rho / this->MW_M * parant;
  if (std::isnan(this->drhodT) || std::isinf(this->drhodT)) {
    // replace with either liquid or vapor derivative
    this->drhodT = (this->betaV > 0.5) ? this->drhodT_V : this->drhodT_L;
  }
  // if still NaN or Inf, throw error
  if (std::isnan(this->drhodT) || std::isinf(this->drhodT)) {
    printf("T_in = %.4f, P_in = %.4f, X_in[0] = %.4e, parant = %g\n",
           T_in, P_in, X_in[0], parant);
    std::cerr << "Error: PengRobinson::computeAdditionalDerivatives: "
              << "drhodT is NaN or Inf: " << this->drhodT << std::endl;
    throw(-1);
  }
  // compute drhodp from (2.2.50)
  parant = this->dBetadp * (M_over_rho_V - M_over_rho_L);
  if (add_vapor) parant += this->betaV / this->rho_V / this->rho_V * (this->rho_V * this->dMWVdp - this->MW_V * this->drhodp_V); 
  if (add_liquid) parant += (1. - this->betaV) / this->rho_L / this->rho_L * (this->rho_L * this->dMWLdp - this->MW_L * this->drhodp_L);
  this->drhodp = - this->rho * this->rho / this->MW_M * parant;
  if (std::isnan(this->drhodp) || std::isinf(this->drhodp)) {
    // replace with either liquid or vapor derivative
    this->drhodp = (this->betaV > 0.5) ? this->drhodp_V : this->drhodp_L;
  }
  // if still NaN or Inf, throw error
  if (std::isnan(this->drhodp) || std::isinf(this->drhodp)) {
    printf("T_in = %.4f, P_in = %.4f, X_in[0] = %.4e, parant = %g\n",
           T_in, P_in, X_in[0], parant);
    std::cerr << "Error: PengRobinson::computeAdditionalDerivatives: "
              << "drhodp is NaN or Inf: " << this->drhodp << std::endl;
    throw(-1);
  }

  this->dpdrho = 1.0 / this->drhodp; // compute dpdrho from drhodp
  // compute dpdT from (K.0.32)
  this->dpdT = - this->drhodT / this->drhodp;
  // compute dZdT from (K.0.40)
  this->dZdT = this->P * this->MW_M / this->gas_constant_universal * (this->rho - this->T * this->drhodT) / pow(this->rho * this->T, 2.);
  // compute dZdp from (K.0.41)
  this->dZdp = this->MW_M / this->gas_constant_universal / this->T * (this->rho - this->P * this->drhodp);
  // compute alpha_p kappa_T (2.2.21)
  // this->alpha_p = 1/this->rho * this->dpdT / this->dpdrho;
  this->alpha_p = - this->drhodT / this->rho;
  this->kappa_T = 1.0 /this->rho / this->dpdrho;
  SAFEGUARD(this->alpha_p, "alpha_p");
  SAFEGUARD(this->kappa_T, "kappa_T");
  // compute C_p from (2.2.38)
  this->cp = - this->P * this->MW_M / (this->rho * this->rho) * this->drhodT; // this is in J/mol/K
  if (add_liquid) this->cp += (1.0 - this->betaV) * this->dEdT_L;
  if (add_vapor)  this->cp += this->betaV * this->dEdT_V;
  if (std::isnan(this->cp) || std::isinf(this->cp)) {
    // replace with either liquid or vapor derivative
    this->cp = (this->betaV > 0.5) ? this->dEdT_V : this->dEdT_L;
    this->cp -= this->P * this->MW_M / (this->rho * this->rho) * this->drhodT;
  }
  SAFEGUARD(this->cp, "cp");
  // compute C_v from (2.2.25)
  this->cv = this->cp - this->T * this->MW_M / this->rho * this->alpha_p * this->alpha_p / this->kappa_T;
  // compute gamma = cp/cv
  this->gamma = this->cp / this->cv;
  // compute kappa_S (2.2.22)
  this->kappa_S = this->kappa_T - this->T * this->alpha_p * this->alpha_p / this->rho / (this->cp / this->MW_M);
  // compute sos from (2.2.23)
  double sos2 = 1 / (this->kappa_S * this->rho);
  if (sos2 < 0.0) {
    // use (3.4.16)
    // compute sos2_L and sos2_V
    double sos2_L = this->gamma * this->dpdrho_L;
    double sos2_V = this->gamma * this->dpdrho_V;
    sos2 = 1. / (this->betaV * this->betaV * this->gamma / sos2_V
                  + (1.0 - this->betaV) * (1.0 - this->betaV) * this->gamma / sos2_L
                  + this->betaV * (1. - this->betaV) * this->rho_L * (1/this->P + this->kappa_T));

    if (sos2 < 0.0) { // use Wood et al (3.4.17)
      sos2 = 1. / this->rho / (this->betaV / this->rho_V / sos2_V + (1.-this->betaV) / this->rho_L / sos2_L);
    }
  }
  if (sos2 < 0.0) {
    std::cerr << "Error: PengRobinson::computeAdditionalDerivatives: "
              << "sos2 is negative: " << sos2 << std::endl;
    exit(-1);
  }
  this->sos = sqrt(sos2); // speed of sound
} // end of PengRobinson::computeAdditionalDerivatives
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\brief: compute variations in T,P,X, x, y
// @param in T_in: temperature
// @param in P_in: pressure
// @param in X_in: mixture molar fractions
// @param in beta_in: vapor volume fraction
// @param in x_in: liquid molar fractions
// @param in y_in: vapor molar fractions
void PengRobinson::computeTPvariations(const double T_in, const double P_in, const double* X_in,
      const double beta_in, const double* x_in, const double* y_in) {
  // Algorithm 3
  // Obtain liquid and vapor mixtures
  this->generateVLstructure(T_in, P_in, x_in, iLIQPH_TP);
  this->computeFugacityDerivatives(x_in, iLIQPH_TP); // dlogPhiL/dy, dlogPhiL/dT, dlogPhiL/dp

  this->generateVLstructure(T_in, P_in, y_in, iVAPPH_TP);
  this->computeFugacityDerivatives(y_in, iVAPPH_TP); // dlogPhiV/dx, dlogPhiV/dT, dlogPhiV/dp
  // compute scriptA, b_T and b_P (K.0.23-25)
  Eigen::MatrixXd scriptA; scriptA.resize(this->num_species, this->num_species);
  Eigen::VectorXd b_T, b_p; b_T.resize(this->num_species); b_p.resize(this->num_species);
  LOOP_l_N (this->num_species) {
    int i = l; // to help consistency with reference
    b_T(i) = beta_in * (1. - beta_in) * (this->d_logPhiV_dT[i] - this->d_logPhiL_dT[i]);
    b_p(i) = beta_in * (1. - beta_in) * (this->d_logPhiV_dp[i] - this->d_logPhiL_dp[i]);
    LOOP_k_N (this->num_species) {
      int apos = i * this->num_species + k;
      double delta_ki = (k == l) ? 1.0 : 0.0; // Kronecker delta
      double sumj = 0;
      for (int j = 0; j < this->num_species; j++) {
        int apos_ij = i * this->num_species + j;
        double delta_jk = (j == k) ? 1.0 : 0.0; // Kronecker delta
        sumj += this->y_TP[j] / beta_in * this->d_logPhiV_dy[apos_ij] * (1 - delta_jk/this->y_TP[j])
                + this->x_TP[j] / (1 - beta_in) * this->d_logPhiL_dx[apos_ij] * (1 - delta_jk/this->x_TP[j]);
      }
      scriptA(i, k) = (1 - delta_ki * this->X[k] / this->y_TP[k] / this->x_TP[k])
                      + beta_in * (1 - beta_in) * sumj;
      // SAFEGUARD(scriptA(i, k), "scriptA[" + std::to_string(i) + ", " + std::to_string(k) + "]");
    }
  }
  // solve scriptA chi = b
  Eigen::MatrixXd scriptA_inv = scriptA.inverse();
  Eigen::VectorXd chi_T = scriptA_inv * b_T; // chi = scriptA^-1 * b_T
  Eigen::VectorXd chi_P = scriptA_inv * b_p; // chi = scriptA^-1 * b_P
  if (beta_in == 1.0 || beta_in == 0.0) {
    chi_T.setZero(); // zero-out chi_T and chi_P if beta_in is 0 or 1
    chi_P.setZero(); // zero-out chi_T and chi_P if beta_in is 0 or 1
  }
  // compute derivatives (k.0.7. K.0.8, K.0.16-18)
  // Zero-out
  this->dBetadT = 0.0; this->dBetadp = 0.0;
  LOOP_l_N (this->num_species) {
    this->dxdp[l] = 0.0;
    this->dydp[l] = 0.0;
    this->dxdT[l] = 0.0;
    this->dydT[l] = 0.0;
    this->dKdp[l] = 0.0;
    this->dKdT[l] = 0.0;
  }

  LOOP_k_N (this->num_species) {
    this->dBetadT += chi_T(k);
    this->dBetadp += chi_P(k);
    LOOP_l_N (this->num_species) {
      int i = l;
      int apos = i * this->num_species + k;
      double delta_ik = (k == l) ? 1.0 : 0.0; // Kronecker delta
      this->dxdp[i] += this->x_TP[i] / (1. - beta_in) * (1. - delta_ik / this->x_TP[k]) * chi_P(k); REPLACEIFWEIRD(this->dxdp[i]);
      this->dydp[i] -= this->y_TP[i] / beta_in * (1. - delta_ik / this->y_TP[k]) * chi_P(k); REPLACEIFWEIRD(this->dydp[i]);
      this->dxdT[i] += this->x_TP[i] / (1. - beta_in) * (1. - delta_ik / this->x_TP[k]) * chi_T(k); REPLACEIFWEIRD(this->dxdT[i]);
      this->dydT[i] -= this->y_TP[i] / beta_in * (1. - delta_ik / this->y_TP[k]) * chi_T(k); REPLACEIFWEIRD(this->dydT[i]);
      this->dKdp[i] += 1 / (beta_in*(1-beta_in)*this->x_TP[i]) * 
                       chi_P(k) * (this->X[i] / this->x_TP[i] * delta_ik - this->y_TP[k]); REPLACEIFWEIRD(this->dKdp[i]);
      this->dKdT[i] += 1 / (beta_in*(1-beta_in)*this->x_TP[i]) * 
                       chi_T(k) * (this->X[i] / this->x_TP[i] * delta_ik - this->y_TP[k]); REPLACEIFWEIRD(this->dKdT[i]);

      SAFEGUARD(this->dxdT[i], "dxdT[" + std::to_string(i) + "]");
      SAFEGUARD(this->dydT[i], "dydT[" + std::to_string(i) + "]");
      SAFEGUARD(this->dxdp[i], "dxdp[" + std::to_string(i) + "]");
      SAFEGUARD(this->dydp[i], "dydp[" + std::to_string(i) + "]");
    }
  }

  // Derivatives of MW_V and MW_L (K.0.28) - (K.0.31)
  this->dMWLdp = 0.0; this->dMWLdT = 0.0;
  this->dMWVdp = 0.0; this->dMWVdT = 0.0;
  // (2.2.40)
  double v_species [this->num_species];
  LOOP_k_N (this->num_species) v_species[k] = beta_in * y_in[k];

  double betaSq = beta_in * beta_in;
  double betaPSq = (1.0 - beta_in) * (1.0 - beta_in);
  // LOOP_k_N (this->num_species) {
  //   this->dMWVdT += (chi_T(k) - v_species[k] * this->dBetadT) * this->MW[k];
  //   this->dMWVdp += (chi_P(k) - v_species[k] * this->dBetadp) * this->MW[k];
  //   this->dMWLdT += ((this->X[k] - v_species[k]) * this->dBetadT - chi_T(k) * (1-beta_in)) * this->MW[k];
  //   this->dMWLdp += ((this->X[k] - v_species[k]) * this->dBetadp - chi_P(k) * (1-beta_in)) * this->MW[k];
  // }
  // this->dMWVdT /= betaSq;
  // this->dMWVdp /= betaSq;
  // this->dMWLdT /= betaPSq;
  // this->dMWLdp /= betaPSq;
  LOOP_k_N (this->num_species) {
    this->dMWVdT += this->dydT[k] * this->MW[k];
    this->dMWVdp += this->dydp[k] * this->MW[k];
    this->dMWLdT += this->dxdT[k] * this->MW[k];
    this->dMWLdp += this->dxdp[k] * this->MW[k];
  }
  // throw(-1);

} // end of PengRobinson::computeTPvariations
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\paramin chi_in: input molar fractions
// @param phase: phase type (iLIQPH_TP or iVAPPH_TP) to set fugacity derivatives
// @\brief: compute fugacity derivatives, d_logPhi*
// ----------------------------------------------------------------------------
void PengRobinson::computeFugacityDerivatives(const double* chi_in, const int phase) {
  // Copy molar fractions so that we can restore them later
  double X_copy [this->num_species];
  LOOP_k_N (this->num_species) X_copy[k] = this->X[k];

  // variables to be set altered; generalize
  double* m_d_logPhi_dXk;
  double* m_d_logPhi_dT;
  double* m_d_logPhi_dp;
  // assume T, P were already set, such that
  this->SetMolarFraction(chi_in); // updates X, Y, MW_M
  this->gas_constant = this->gas_constant_universal / this->MW_M; // update gas constant
  this->Am = this->GetAm(this->T); // make sure Am Bm are computed
  this->Bm = this->GetBm();

  // determine which quantities to update based on input phase
  double m_rho, Z, m_MW;
  m_MW = this->MW_M; // automatically updated by setting chi_in
  if (phase == this->iLIQPH_TP) {
    m_rho = this->rho_L; Z = this->Z_L;
    m_d_logPhi_dXk = this->d_logPhiL_dx;
    m_d_logPhi_dT = this->d_logPhiL_dT;
    m_d_logPhi_dp = this->d_logPhiL_dp;
  } else if (phase == this->iVAPPH_TP) {
    m_rho = this->rho_V; Z = this->Z_V;
    m_d_logPhi_dXk = this->d_logPhiV_dy;
    m_d_logPhi_dT = this->d_logPhiV_dT;
    m_d_logPhi_dp = this->d_logPhiV_dp;
  } else {
    // Fugacity derivatives only done for liquid and vapor phases
    std::cerr << "Error: PengRobinson::computeFugacityDerivatives: "
              << "Unknown phase type: " << phase << std::endl;
    throw (-1);
  }

  this->computeBasicDerivatives(phase); // d*dT, d*dp, d*dXk. Based on input chi_in

  // check out (Sec. M.3.2); Good Luck ;)
  // Start computing!
  // d_logPhiV_dy, d_logPhiL_dx
  double term1, term2, term3, term4, term5;
  double scriptU, scriptX, scriptY, scriptW, scriptJ;
  double scriptD, scriptE;
  double scriptK [this->num_species];
  double scriptT [this->num_species];
  double scriptV [this->num_species];
  double scriptJDeltai [this->num_species];
  double dscriptUdXk [this->num_species];
  double dscriptKdXk [this->num_species * this->num_species];
  double dscriptVdXk [this->num_species * this->num_species];
  double dscriptTdXk [this->num_species * this->num_species];
  double dscriptWdXk [this->num_species];
  double dscriptXdXk [this->num_species];
  double dscriptYdXk [this->num_species];
  double dscriptDdXk [this->num_species];
  double dscriptEdXk [this->num_species];
  double dscriptJDeltai_dXk [this->num_species * this->num_species];

  double Aij [this->num_species * this->num_species];
  double dAijdT [this->num_species * this->num_species];

  double R_u = this->gas_constant_universal;
  double RuT_p = this->gas_constant_universal * this->T / this->P;
  double Bm_tilde = this->Bm / RuT_p;
  double dBmtildedXk [this->num_species];
  LOOP_k_N (this->num_species) dBmtildedXk[k] = this->dBmdXk[k] / RuT_p;

  // (M.3.23)
  scriptU = this->Am / this->Bm / this->gas_constant_universal / this->T / (this->delta2 - this->delta1);
  scriptX = Z + this->delta2 * Bm_tilde; scriptY = Z + this->delta1 * Bm_tilde;
  scriptW = scriptX/scriptY;
  scriptD = this->Bm / (delta2*Bm_tilde + Z) / (delta1*Bm_tilde + Z);

  double S1 [this->num_species]; double S2 [this->num_species];
  LOOP_k_N (this->num_species) {S1[k] = 0.0; S2[k] = 0.0;}

  LOOP_k_N (this->num_species) {
    int i = k;
    LOOP_l_N (this->num_species) {
      int j = l;
      double delta_ij = (i == j) ? 1.0 : 0.0; // Kronecker delta S1[i] += this->delta1 * (delta_ij - this->X[j]);
      S2[i] += this->delta2 * (delta_ij - this->X[j]);
      int apos = i * this->num_species + j;
      // Determine Aij based on mixing rule
      if (this->mixing_rule == "VDW") {
        double f_i = this->cst_c[i*i]; double Tcrit_i = this->Tcrit[i];
        double alpha_i = 1. + f_i * (1. - sqrt(this->T/Tcrit_i)); alpha_i = alpha_i * alpha_i;
        double dalphaidT = - f_i / sqrt(Tcrit_i * this->T) * sqrt(alpha_i);

        double f_j = this->cst_c[j*j]; double Tcrit_j = this->Tcrit[j];
        double alpha_j = 1. + f_j * (1. - sqrt(this->T/Tcrit_j)); alpha_j = alpha_j * alpha_j;
        double dalphajdT = - f_j / sqrt(Tcrit_j * this->T) * sqrt(alpha_j);

        Aij[apos] = sqrt(this->cst_a[i] * this->cst_a[j] * alpha_i * alpha_j);
        dAijdT[apos] = 0.5 * sqrt(this->cst_a[i] * this->cst_a[j] / alpha_i / alpha_j) 
                      * (alpha_j * dalphaidT + alpha_i * dalphajdT);
      } else if (this->mixing_rule == "MR3") {
        double f_ij = this->cst_c[apos]; double Tcrit_ij = this->Tcrit_IJ[apos];
        double prod = 1. + f_ij * (1. - sqrt(this->T/Tcrit_ij));
        double alpha_ij = prod * prod;
        double A_IJ = this->cst_a[apos] * alpha_ij;
        double dalphadT_IJ = - f_ij / sqrt(Tcrit_ij * this->T) * prod;
        Aij[apos] = A_IJ;
        dAijdT[apos] = this->cst_a[apos] * dalphadT_IJ;
      }
    }
  }

  double denom = pow((delta2*Bm_tilde+Z)*(delta1*Bm_tilde+Z),2);
  // real computation
  // d_logPhi*_dx/y
  LOOP_l_N (this->num_species) {
    int i = l; // to make it easy..
    double S2i = S2[i];
    double S1i = S1[i];

    double sum_XjAji = 0.0;
    double sum_XjdAjidT = 0.0;
    for (int j = 0; j < this->num_species; j++) {
      int apos_j = j * this->num_species + i;
      double Aji = Aij[apos_j];
      double dAji_dT = dAijdT[apos_j];

      sum_XjAji += this->X[j] * Aji;
      sum_XjdAjidT += this->X[j] * dAji_dT;
    }
    scriptE = Bm_tilde * (S2i * delta1 - S1i * delta2) + Z * (S2i - S1i);
    scriptV[i] = 2 * sum_XjAji / this->Am - this->cst_b[i] / this->Bm;
    scriptT[i] = scriptU * scriptV[i];
    scriptK[i] = (S2i - S1i) / (this->delta2 - this->delta1);
    scriptJ = pow(RuT_p, -2.0) / m_rho / (delta2 * Bm_tilde + Z) / (delta1 * Bm_tilde + Z);
    scriptJDeltai[i] = 1/RuT_p * scriptD * scriptE;
    // l = i
    LOOP_k_N (this->num_species) {
      // k  = k
      int apos = l * this->num_species + k;
      double delta_ij = (k == l) ? 1.0 : 0.0;
      double Bik = 0.0; // (G.0.5)
      dscriptUdXk[k] = 1/(R_u * this->T * (this->delta2 - this->delta1))
        * (this->dAmdXk[k] / this->Bm - this->Am * this->dBmdXk[k] / (this->Bm * this->Bm))
        - this->Am / this->Bm / this->gas_constant_universal / this->T * (this->delta2 - this->delta1);
      dscriptXdXk[k] = this->dZdXk[k] + this->delta2 / RuT_p * this->dBmdXk[k] + Bm_tilde * this->delta2;
      dscriptYdXk[k] = this->dZdXk[k] + this->delta1 / RuT_p * this->dBmdXk[k] + Bm_tilde * this->delta1;
      dscriptWdXk[k] = dscriptXdXk[k] / scriptY - scriptX / (scriptY * scriptY) * dscriptYdXk[k];

      dscriptDdXk[k] = this->dBmdXk[k] / (delta2*Bm_tilde + Z) / (delta1*Bm_tilde + Z) - this->Bm * Bm_tilde * Bm_tilde * (2.*this->delta1delta2) / denom
        - this->Bm * dBmtildedXk[k] * (2 *this->delta1delta2 * Bm_tilde + Z) / denom
        - this->Bm * this->dZdXk[k] * (this->delta1_p_delta2 * Bm_tilde + 2*Z) / denom
        - this->Bm * Bm_tilde * Z * (this->delta1_p_delta2) / denom;
      dscriptEdXk[k] = dBmtildedXk[k] * (S2i * delta1 - S1i * delta2) 
                       + Bm_tilde * (delta1*(S2i+delta2) - delta2 * (S1i+delta1))
                       + this->dZdXk[k] * (S2i - S1i) + Z * (delta1 - delta2);

      dscriptKdXk[apos] = -(delta2-delta1) * ((delta2-delta1) + S2i-S1i) / pow((this->delta2 - this->delta1), 2);

      double sum_XjAjidAmdXk = 0.0;
      for (int j = 0; j < this->num_species; j++) {
        int apos_j = j * this->num_species + i;
        double Aji = Aij[apos_j];
        sum_XjAjidAmdXk += this->X[j] * Aji * this->dAmdXk[k];
      }

      dscriptVdXk[apos] = - (Bik * this->Bm - this->cst_b[i] * this->dBmdXk[k]) / (this->Bm * this->Bm)
        + (2 * this->Am * Aij[k*this->num_species+i] - 2 * sum_XjAjidAmdXk) / this->Am / this->Am;

      dscriptTdXk[apos] = dscriptUdXk[k] * scriptV[i] + scriptU * dscriptVdXk[apos];
      dscriptJDeltai_dXk[apos] = 1/RuT_p * (scriptD * dscriptEdXk[k] + scriptE * dscriptDdXk[k]);

      // computation of all the terms
      term1 = this->cst_b[i] / this->Bm * this->dZdXk[k]
        + (Z-1.0) * (Bik * this->Bm - this->cst_b[i] * this->dBmdXk[k]) / (this->Bm * this->Bm); REPLACEIFWEIRD(term1);
      term2 = 1. / (Z - Bm_tilde) * (this->dZdXk[k] - dBmtildedXk[k]); REPLACEIFWEIRD(term2);
      term3 = dscriptTdXk[apos] * log(scriptW) + scriptT[i] / scriptW * dscriptWdXk[k]; REPLACEIFWEIRD(term3);
      term4 = scriptK[i] * scriptU/scriptW * dscriptWdXk[k] + scriptK[i] * dscriptUdXk[k] * log(scriptW)
              + scriptU * log(scriptW) * dscriptKdXk[apos]; REPLACEIFWEIRD(term4);
      term5 = scriptJDeltai[i] * dscriptUdXk[k] + scriptU * dscriptJDeltai_dXk[apos]; REPLACEIFWEIRD(term5);
      // set d_logPhi*_dXk
      m_d_logPhi_dXk[apos] = term1 - term2 - term3 + term4 - term5;
      SAFEGUARD(m_d_logPhi_dXk[apos], "d_logPhi_dXk[" + std::to_string(apos) + "]");
    }
    // d_logPhi*_dT
    double Deltai = this->Bm * m_MW / Z * (Bm_tilde * (S2i*delta1-S1i*delta2) + Z * (S2i - S1i)); // (M.3.13)
    double T2 = this->T * this->T;
    double T3 = T2 * this->T;
    double dscriptUdT = 1/(this->Bm * R_u * (this->delta2-this->delta1)) / this->T / this->T * (this->dAmdT * this->T - this->Am);
    double dscriptVdT = (2. * this->Am * sum_XjdAjidT - this->dAmdT * 2 * sum_XjAji) / this->Am / this->Am;
    double term3Bracket = Bm_tilde * (this->delta1 - this->delta2) * (this->dZdT + Z / this->T);

    double dDeltaidT = - this->Bm * this->Bm * m_rho / Z / this->T * (S2i * this->delta1 - S1i * this->delta2) * (this->dZdT + Z);
    double dBmtildedT = - Bm_tilde / this->T;
    double dscriptJdT = this->P / (R_u * T3 * denom) * (Bm_tilde * (Z*this->delta1_p_delta2 + this->delta1delta2) + Z * Z) * (this->dZdT * this->T - Z)
                      - this->P * Z / (R_u / T2 * denom) * (dBmtildedT) * (2*this->delta1delta2 * Bm_tilde + Z * this->delta1_p_delta2)
                      - this->P * Z / (R_u / T * denom) * this->dZdT * (Bm_tilde * this->delta1_p_delta2 + 2 * Z);

    term1 = this->cst_b[i] / this->Bm * this->dZdT; REPLACEIFWEIRD(term1);
    term2 = 1/(Z - Bm_tilde) * (this->dZdT + Bm_tilde / this->T); REPLACEIFWEIRD(term2);
    double term3a = (dscriptUdT * scriptV[i] + scriptU * dscriptVdT) * log(scriptX/scriptY);
    double term3b = scriptU * scriptV[i] / scriptX / scriptY * term3Bracket;

    term3 = (dscriptUdT * scriptV[i] + scriptU * dscriptVdT) * log(scriptX/scriptY)
      + scriptU * scriptV[i] / scriptX / scriptY * term3Bracket; REPLACEIFWEIRD(term3);
    term4 = scriptK[i] *
            ( dscriptUdT * log(scriptX/scriptY)
              + scriptU / (scriptX * scriptY) * term3Bracket); REPLACEIFWEIRD(term4);
    term5 = scriptU * scriptJ * dDeltaidT
            + scriptU * Deltai * dscriptJdT
            + scriptJ * Deltai * dscriptUdT; REPLACEIFWEIRD(term5);
    m_d_logPhi_dT[i] = term1 - term2 - term3 + term4 - term5;
    SAFEGUARD(m_d_logPhi_dT[i], "d_logPhi_dT[" + std::to_string(i) + "]");
    // printf("scriptV[%d] = %e, dscriptVdT = %e\n",
    //        i, scriptV[i], dscriptVdT);
    // printf("term1 = %e, term2 = %e, term3 = %e, term4 = %e, term5 = %e\n", term1, term2, term3, term4, term5);
    // printf("d_logPhi_dT[%d] = %e\n", i, m_d_logPhi_dT[i]);
    // if (phase == this->iLIQPH_TP && i == 1) throw(-1); // debug

    // d_logPhi*_dp
    double dscriptXdp = this->dZdp + this->delta2 * this->Bm / R_u / this->T;
    double dscriptYdp = this->dZdp + this->delta1 * this->Bm / R_u / this->T;
    double dDeltaidp = this->Bm * this->Bm * m_MW / Z / R_u / this->T * (S2i * this->delta1 - S1i * this->delta2) * (1.0 - this->P / Z * this->dZdp);
    double dBmtildedp = Bm_tilde / this->P;
    double dscriptJdp = (Z + this->P * this->dZdp) / R_u / m_MW / this->T / sqrt(denom)
                       - this->P * Z / (R_u * m_MW * this->T * denom) * ((this->delta2 * dBmtildedp + this->dZdp) * (this->delta1 * Bm_tilde + Z)
                       - (this->delta2 * Bm_tilde + Z) * (this->delta1 * dBmtildedp + this->dZdp));

    term1 = this->cst_b[i] / this->Bm * this->dZdp; REPLACEIFWEIRD(term1);
    term2 = 1. / (Z - Bm_tilde) * (this->dZdp - Bm_tilde / this->P); REPLACEIFWEIRD(term2);
    term3 = scriptU * scriptV[i] / scriptX / scriptY * (dscriptXdp * scriptY - scriptX * dscriptYdp); REPLACEIFWEIRD(term3);
    term4 = scriptU / scriptX / scriptY * (dscriptXdp * scriptY - scriptX * dscriptYdp); REPLACEIFWEIRD(term4);
    term5 = scriptU * scriptJ * dDeltaidp + scriptU * Deltai * dscriptJdp; REPLACEIFWEIRD(term5);
    m_d_logPhi_dp[i] = term1 - term2 - term3 + term4 - term5;
    SAFEGUARD(m_d_logPhi_dp[i], "d_logPhi_dp[" + std::to_string(i) + "]");
  } // end of loop i

  // Restore molar fractions
  this->SetMolarFraction(X_copy);
  this->gas_constant = this->gas_constant_universal / this->MW_M; // update gas constant
} // end of PengRobinson::computeFugacityDerivatives
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @\paramin T,P,X,phase,chi_in
// @\brief: compute phase derivatives
// ----------------------------------------------------------------------------
void PengRobinson::computeAdditionalPhaseDerivatives
(const double T_in, const double P_in, const double* X_in, const int phase, const double* chi_in) {
  // copy current molar fractions to set later
  double X_copy [this->num_species];
  LOOP_k_N (this->num_species) {
    X_copy[k] = this->X[k];
    assert(X_copy[k] == X_in[k]);
  }
  this->SetMolarFraction(chi_in); // set molar fractions
  this->gas_constant = this->gas_constant_universal / this->MW_M; // update gas constant
  
  // Algorithm 5
  double m_dpdrho, m_dpdT, m_dZdT, m_dZdp;
  double *m_drhodchi;
  double *m_dZdchi;
  double *m_dchidT, *m_dchidp;
  double m_Z, m_rho, m_MW;
  double dDelEdelT, dDelEdelP;
  double *dDelEdelChi;
  if (phase == this->iLIQPH_TP) {
    m_rho = this->rho_L; m_Z = this->Z_L; m_MW = this->MW_L;
    m_dpdrho = this->dpdrho_L;
    m_dpdT = this->dpdT_L;
    m_drhodchi = this->drhoLdx;
    m_dZdT = this->dZdT_L;
    m_dZdp = this->dZdp_L;
    m_dZdchi = this->dZLdx;
    dDelEdelT = this->dDepEdT_L;
    dDelEdelP = this->dDepEdp_L;
    dDelEdelChi = this->dDepEdXk_L;
    m_dchidT = this->dxdT; m_dchidp = this->dxdp;
  } else if (phase == this->iVAPPH_TP) {
    m_rho = this->rho_V; m_Z = this->Z_V; m_MW = this->MW_V;
    m_dpdrho = this->dpdrho_V;
    m_dpdT = this->dpdT_V;
    m_drhodchi = this->drhoVdy;
    m_dZdT = this->dZdT_V;
    m_dZdp = this->dZdp_V;
    m_dZdchi = this->dZVdy;
    dDelEdelT = this->dDepEdT_V;
    dDelEdelP = this->dDepEdp_V;
    dDelEdelChi = this->dDepEdXk_V;
    m_dchidT = this->dydT; m_dchidp = this->dydp;
  } else {
    std::cerr << "Error: PengRobinson::computeAdditionalPhaseDerivatives: "
              << "Unknown phase type: " << phase << std::endl;
    throw (-1);
  }

  if (std::isnan(m_MW) || std::isinf(m_MW)) m_MW = 0.0;

  // (K.0.37, 34, 35)
  // Compute dZdchi, dzdT, dZdp (K.0.38,39)
  double drhoEtadp = 1. / m_dpdrho;
  double drhoEtadT = - m_dpdT / m_dpdrho;
  LOOP_k_N (this->num_species) {
    drhoEtadT += m_drhodchi[k] * m_dchidT[k];
    drhoEtadp += m_drhodchi[k] * m_dchidp[k];
    m_dZdT += m_dZdchi[k] * m_dchidT[k];
    m_dZdp += m_dZdchi[k] * m_dchidp[k];
  }
  if (phase == this->iLIQPH_TP) {
    this->dZdT_L = m_dZdT; // set derivatives for liquid phase
    this->dZdp_L = m_dZdp; // set derivatives for liquid phase
    this->drhodT_L = drhoEtadT; // set derivatives for liquid phase
    this->drhodp_L = drhoEtadp; // set derivatives for liquid phase
  } else if (phase == this->iVAPPH_TP) {
    this->dZdT_V = m_dZdT; // set derivatives for liquid phase
    this->dZdp_V = m_dZdp; // set derivatives for liquid phase
    this->drhodT_V = drhoEtadT; // set derivatives for liquid phase
    this->drhodp_V = drhoEtadp; // set derivatives for liquid phase
  }
  // Ideal gas energy and Cp
  double Cv_ig = this->GetIdealGasCp(T_in, true) - this->gas_constant_universal; // J/kmol-K
  double E_eta_ig_v [this->num_species];
  LOOP_k_N (this->num_species) {
    E_eta_ig_v[k] = this->hspecies_ig[k] * T_in * this->gas_constant_universal;
  }
  // energy departure gradients (K.0.46, K.0.47, 48, 44, 45)
  double dE_eta_dChi [this->num_species];

  LOOP_k_N(this->num_species) {
    dE_eta_dChi[k] = E_eta_ig_v[k] + (dDelEdelChi[k]); // (K.0.46)
    REPLACEIFWEIRD(dE_eta_dChi[k]);
  }

  double dE_eta_dT = Cv_ig + dDelEdelT; // (K.0.47)
  double dE_eta_dp = dDelEdelP; // (K.0.48)

  double surplus_T = 0.0;
  double surplus_p = 0.0;
  // pressure and temperature gradients added with dchid*
  LOOP_k_N (this->num_species) {
    surplus_T += dE_eta_dChi[k] * m_dchidT[k]; // (K.0.44)
    surplus_p += dE_eta_dChi[k] * m_dchidp[k]; // (K.0.45)
  }
  dE_eta_dT += surplus_T; // (K.0.44)
  dE_eta_dp += surplus_p; // (K.0.45)

  SAFEGUARD(dE_eta_dT, "dE_eta_dT");
  SAFEGUARD(dE_eta_dp, "dE_eta_dp");
  // if (dE_eta_dT < 0.0) {
  //   std::cerr << "Error: PengRobinson::computeAdditionalPhaseDerivatives: "
  //             << "dE_eta_dT < 0.0: " << dE_eta_dT << std::endl;
  //   throw(-1);
  // }

  if (phase == this->iLIQPH_TP) {
    this->dEdT_L = dE_eta_dT; // set derivatives for liquid phase
    this->dEdp_L = dE_eta_dp; // set derivatives for liquid phase
    LOOP_k_N (this->num_species) {
      this->dEdchi_L[k] = dE_eta_dChi[k]; // set derivatives for liquid phase
    }
  } else if (phase == this->iVAPPH_TP) {
    this->dEdT_V = dE_eta_dT; // set derivatives for liquid phase
    this->dEdp_V = dE_eta_dp; // set derivatives for liquid phase
    LOOP_k_N (this->num_species) {
      this->dEdchi_V[k] = dE_eta_dChi[k]; // set derivatives for liquid phase
    }
  }
  // restore molar fractions
  this->SetMolarFraction(X_copy); // restore molar fractions
  this->gas_constant = this->gas_constant_universal / this->MW_M; //
} // end of PengRobinson::computeAdditionalPhaseDerivatives

double PengRobinson::GetDensityFromPressureTemperature_TP(const double p_in,
                                                    const double T_in, const int phase) {

  double Z;
  double m_MW;
  if (phase == this->iLIQPH_TP) {
    m_MW = this->MW_L;
  } else if (phase == this->iVAPPH_TP) {
    m_MW = this->MW_V;
  } else {
    m_MW = this->MW_M; // mixture
  }
  double a = this->GetAm(T_in); // make sure Am Bm are computed
  double b = this->GetBm();
  this->Am = a; this->Bm = b;
  double Bm_tilde = this->Bm * this->P / (this->gas_constant_universal) / this->T;
  double Ru = this->gas_constant_universal;

  // check notes, or Tudisco thesis (2021) (2.2.71)
  double a1 = pow(Ru * T_in / p_in, 3);
  double a2 = a1 * (p_in / Ru / T_in * b * (this->delta1_p_delta2 - 1.) - 1.);
  double a3 = Ru * T * b * b / p_in * (this->delta1delta2 - delta1_p_delta2)
    - this->delta1_p_delta2 * b * pow(Ru * T_in / p_in, 2)
    + a * Ru * T_in / p_in / p_in;
  double a4 = - this->delta1delta2 * Ru * T_in * b * b / p_in
    - a * b / p_in - this->delta1delta2 * b * b * b;

  double A = a2 / a1; double B = a3 / a1; double C = a4 / a1;
  std::vector<double> xZ(3);
  std::vector<double> Z_physical;
  int n = SolveP3(&xZ[0], A, B, C);
  if (n == 1) { // one root
    Z = xZ[0];
  } else {
    std::sort(xZ.begin(), xZ.end()); // this sorting is important to get the correct roots
    for (int i = 0; i < n; i++) {
      if (xZ[i] > Bm_tilde) {
        Z_physical.push_back(xZ[i]);
      }
    }

    double Z_L, Z_V;
    if (Z_physical.size() == 0) {
      Z_L = 0.0; // liquid phase
      Z_V = 0.0; // liquid phase
    } else {
      Z_L = Z_physical[0]; // liquid phase
      Z_V = Z_physical[Z_physical.size()-1]; // liquid phase
    }
    if (phase == this->iLIQPH_TP) {
      // Z = Z2;
      Z = Z_L;
    } else if (phase == this->iVAPPH_TP) {
      // Z = Z1;
      Z = Z_V;
    } else {
      double Z2= Z_L; // liquid phase
      double Z1= Z_V; // vapor phase
                      // two-phase, In eed to compute the Gibbs free energy and choose what solution to return
      double deltaG_mix = 0.0; // scaled by Ru * T
                               // Check the Z1 and Z2 are greater than Bm_tilde
                               // check (M.3.18)
      double min_lnPhi = 1e99;
      int min_k = -1;

      // This is the original implementation based on Tudisco
      // LOOP_k_N(this->num_species) {
      //   double sum = 0.0;
      //   LOOP_l_N(this->num_species) {
      //     int apos = l * num_species + k;
      //     sum += this->X[l] * this->cst_a[apos];
      //   }
      //   double lnPhi_k_Z2 =  - log(Z2 - Bm_tilde)
      //     + this->cst_b[k] / this->Bm * (Z2-1)
      //     - log((this->delta2 * Bm_tilde +  Z2)/(this->delta1 * Bm_tilde + Z2))
      //     * this->Am / this->Bm / Ru / this->T / (this->delta2 - this->delta1)
      //     * (2. * sum / this->Am - this->cst_b[k] / this->Bm);
      //   double lnPhi_k_Z1 =  - log(Z1 - Bm_tilde)
      //     + this->cst_b[k] / this->Bm * (Z1-1)
      //     - log((this->delta2 * Bm_tilde +  Z1)/(this->delta1 * Bm_tilde + Z1))
      //     * this->Am / this->Bm / Ru / this->T / (this->delta2 - this->delta1)
      //     * (2. * sum / this->Am - this->cst_b[k] / this->Bm);
      //   deltaG_mix += this->X[k] * (lnPhi_k_Z2 - lnPhi_k_Z1);
      // }

      double lnPhi_k_Z2 = - log(Z2 - Bm_tilde) 
        - this->Am / this->Bm / Ru / T_in / (this->delta1_m_delta2)
        * log((Z2+delta1*Bm_tilde)/(Z2+delta2*Bm_tilde))
        + Z2-1.0;

      double lnPhi_k_Z1 = - log(Z1 - Bm_tilde) 
        - this->Am / this->Bm / Ru / T_in / (this->delta1_m_delta2)
        * log((Z1+delta1*Bm_tilde)/(Z1+delta2*Bm_tilde))
        + Z1-1.0;
      deltaG_mix = lnPhi_k_Z2 - lnPhi_k_Z1;
      assert(std::isnan(deltaG_mix) == false);
      assert(std::isinf(deltaG_mix) == false);
      if (deltaG_mix > 0.0) {
        // return vapor phase
        Z = Z1;
      } else {
        // return liquid phase
        Z = Z2;
      }
    }
  }
  double rho_ = p_in * m_MW / Z / Ru / T_in;
  return rho_;
}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
double PengRobinson::GetAm(double T_in, const bool compute_gradients) {
  double Am = 0.0;
  if (this->mixing_rule == "MR3") {
    LOOP_k_N(this->num_species) {
      LOOP_l_N(this->num_species) { // => theoretically faster
        int apos = k * num_species + l;
        double X_X = this->X[l] * this->X[k];
        double A_IJ = this->cst_a[apos] * pow(
            1.0 + this->cst_c[apos] * (1.0 - sqrt(T_in / this->Tcrit_IJ[apos])),
            2);
        Am += X_X * A_IJ;
      }
    }
  } else if (this->mixing_rule == "VDW") {
    // fill alphas
    LOOP_k_N(this->num_species) {
      int idk = k * k; // diagonal index
      double A_i = this->cst_a[idk]; double f_i = this->cst_c[idk];
      double alpha_i = pow(1 + f_i * (1 - sqrt(T_in / this->Tcrit[k])), 2.0);
      LOOP_l_N(this->num_species) {
        int idl = l * l; // diagonal index
        double A_j = this->cst_a[idl]; double f_j = this->cst_c[idl];
        double alpha_j = pow(1 + f_j * (1 - sqrt(T_in / this->Tcrit[l])), 2.0);
        Am += this->X[k] * this->X[l] * sqrt(A_i * A_j * alpha_i * alpha_j);
      }
    }
  } else {
    std::cout << "Error: PengRobinson::GetAm() - unknown mixture rule: "
              << this->mixing_rule << std::endl;
  }

  // if (compute_gradients) this->computeBasicDerivatives();

  return Am;
}

//----------------------------------------------------------------------------

double PengRobinson::GetBm() {
  double Bm = 0.0;
  LOOP_k_N(this->num_species)
    Bm += X[k] * this->cst_b[k];

  return Bm;
}

void PengRobinson::SetRealFluidConstants() {
  //! \brief Assign arrays of critical point data
  LOOP_k_N(this->num_species) {
    LOOP_l_N(this->num_species) {
      double tmp_k = 0.0; // zero binary interation parameter

      int apos = k * num_species + l;

      this->Tcrit_IJ[apos] =
          sqrt(this->Tcrit[l] * this->Tcrit[k]) * (1.0 - tmp_k);
      this->Vcrit_IJ[apos] =
          pow(pow(this->Vcrit[l], 1.0 / 3.0) + pow(this->Vcrit[k], 1.0 / 3.0),
              3.0) / 8.0;
      this->Zcrit_IJ[apos] = 0.5 * (this->Zcrit[l] + this->Zcrit[k]);
      if (this->mixing_rule == "MR3" ) {
        this->Pcrit_IJ[apos] = this->Zcrit_IJ[apos] * this->gas_constant_universal
          * this->Tcrit_IJ[apos] / this->Vcrit_IJ[apos];
      } else if (this->mixing_rule == "VDW") {
        this->Pcrit_IJ[apos] = sqrt(this->Pcrit[l] * this->Pcrit[k]);
      } else {
        std::cout << "Error: PengRobinson::SetRealFluidConstants() - unknown mixture rule: "
                  << this->mixing_rule << std::endl;
      }
      this->omega_IJ[apos] = 0.5 * (this->omega[l] + this->omega[k]);

    }
  }

  //! \brief Assign constants used for computing real fluid effects
  LOOP_k_N(this->num_species) {
    this->cst_b[k] = 0.077796 * this->gas_constant_universal * this->Tcrit[k]
        / this->Pcrit[k];
    this->dBmdXk[k] = this->cst_b[k];

    LOOP_l_N(this->num_species) {
      int apos = k * num_species + l;
      this->cst_a[apos] = 0.457236
          * pow(this->gas_constant_universal * this->Tcrit_IJ[apos], 2.0)
          / this->Pcrit_IJ[apos];
      this->cst_c[apos] = 0.37464 + 1.54226 * this->omega_IJ[apos]
          - 0.26992 * pow(this->omega_IJ[apos], 2);

    }
  }
}
double PengRobinson::GetIdealGasEnthalpy(double T_in, const bool molar) {
//  const double T_min = 300.0;
//  if (T_in < T_min) T_in = T_min;

  double T1 = T_in;
  double T2 = T_in * T_in;
  double T3 = T2 * T_in;
  double T4 = T3 * T_in;

  this->h_ig = 0.0;
  LOOP_l_N(this->num_species) {
    if (T_in < 1000.0) {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1 / 2.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 3.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 4.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 5.0;
     this->hspecies_ig[l] +=
         this->nasa_poly_coeff[l * this->NasaCoef + 5] / T1;
    } else {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->num_species];
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->num_species] * T1 / 2.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->num_species] * T2 / 3.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->num_species] * T3 / 4.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->num_species] * T4 / 5.0;
     this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 5
         + this->NasaCoef * this->num_species] / T1;
    }
    if (!this->m_includeChem)
      this->hspecies_ig[l] -= this->nasa_poly_coeff[l * this->NasaCoef + 5] / T1;
    this->h_ig += this->X[l] * (hspecies_ig[l]);
  }
  this->h_ig = h_ig * (T_in * this->gas_constant_universal);
  if (!molar) this->h_ig /= this->MW_M; // convert to J/kg
  return this->h_ig;
}

double PengRobinson::GetIdealGasCp(double T_in, const bool molar) {
//  const double T_min = 300.0;
//  if (T_in < T_min) T_in = T_min;

  double T1 = T_in;
  double T2 = T_in * T_in;
  double T3 = T2 * T_in;
  double T4 = T3 * T_in;

  double cp_ig = 0.0;

  LOOP_l_N(this->num_species) {
    if (T_in <= 1000.0) {
      this->cpspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4;
    } else {
      this->cpspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->num_species];
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->num_species] * T1;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->num_species] * T2;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->num_species] * T3;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->num_species] * T4;
    }
    double tmp = cpspecies_ig[l] * this->gas_constant_universal;
    cp_ig += this->X[l] * tmp; //cp_kmol
  }
  if (std::isnan(cp_ig) || std::isinf(cp_ig)) {
    std::cerr << "Error: PengRobinson::GetIdealGasCp: cp_ig is NaN or Inf, T = "
              << T_in << ", P = " << this->P << std::endl;
    throw (-1);
  }
  if (!molar) cp_ig /= this->MW_M; // convert to J/kg-K
  return cp_ig;
}

double PengRobinson::GetEnergyFromTemperature(const double T_in, const bool molar) {
  //! brief Computes the internal energy of real per unit mass [J/kg]
  /* 
  // ORIGINAL CHARLESX
  this->SetRealFluidThermodynamics(this->V, T_in);
  double dep = (this->Am - T_in * this->dAmdT) * this->K1 / this->MW_M;
  */
  if (this->phase == this->iLIQPH_TP && this->MW_L == 0.0) return 0.0;
  if (this->phase == this->iVAPPH_TP && this->MW_V == 0.0) return 0.0;
  // compute dAmdT and d2AmdT2

  double R_u = this->gas_constant_universal;
  double m_rho, m_Z, m_MW;
  if (this->phase == this->iLIQPH_TP) {
    m_rho = this->rho_L; m_Z = this->Z_L; m_MW = this->MW_L;
  } else if (this->phase == this->iVAPPH_TP) {
    m_rho = this->rho_V; m_Z = this->Z_V; m_MW = this->MW_V;
  } else {
    m_rho = this->rho; m_Z = this->P * this->MW_M / this->rho / this->gas_constant_universal / T_in; // mixture
    m_MW = this->MW_M; // mixture
  }

  this->computeBasicDerivatives(this->phase); // to get departure function
  double Bm_tilde = this->Bm * this->P / R_u / this->T;
  // (2.2.79)
  double dep = 1. / (this->delta1 - this->delta2) / this->Bm
              * (this->dAmdT * T_in - this->Am)
              * log((this->delta1 * this->Bm * m_rho + m_MW) / (this->delta2 * this->Bm * m_rho + m_MW));
  // // (2.2.82)
  // double dep2 = 1. / (this->delta1 - this->delta2) / this->Bm
  //             * (this->dAmdT * T_in - this->Am)
  //             * log((this->delta1 * Bm_tilde + m_Z) / (this->delta2 * Bm_tilde + m_Z));
  // if (std::fabs(dep2 - dep) > 1e-6) {
  //   std::cerr << "Error: PengRobinson::GetEnergyFromTemperature: dep2 != dep, dep2 = "
  //             << dep2 << ", dep = " << dep << std::endl;
  //   throw(-1);
  // }

  REPLACEIFWEIRD(dep);
  if (this->phase == this->iLIQPH_TP) {
    this->depE_L = dep; // set departure function for liquid phase
  } else if (this->phase == this->iVAPPH_TP) {
    this->depE_V = dep; // set departure function for vapor phase
  }

  double E_ig_molar = this->GetIdealGasEnthalpy(T_in, true) - this->gas_constant_universal * T_in;
  double E_molar = E_ig_molar + dep;
  double res = E_molar;


  if (!molar) res /= this->MW_M; // convert to J/kg

  return res;

  //energy per unit mass
  // return this->GetIdealGasEnthalpy(T_in) - this->gas_constant * T_in + dep;
  //return this->GetRealFluidEnthalpy(T_in,this->P);
}

//----------------------------------------------------------------------------
void PengRobinson::SetMixture_TPRYF(const double T, const double P, const double RYF, const double YFGuess, const bool flash) {
  assert(this->num_species == 2); // this function is only for two species
  // BK
  // set rho, p and RYF
  // iterate over YF (input guess YF) to get the correct pressure
  this->T = T;
  this->P = P;
  double P_target = P;
  double Y_ [this->num_species];

  auto resFun = [&](double m_YF) {
    // if (std::isnan(m_YF) || std::isinf(m_YF) || m_YF < 0. || m_YF > 1.) {
    //   printf("error in m_YF = %f\n", m_YF);
    //   printf("T = %f, P = %f, RYF = %f\n", this->T, this->P, RYF);
    //   std::cout << std::endl;
    //   exit(1);
    // }
    if (m_YF > 1.0) m_YF = 1.0;
    if (m_YF < 0.0) m_YF = 0.0;
    if (flash) { // flash, this is typically used for two-phase RIM calculation
      double Yin [this->num_species];
      Yin[0] = m_YF; // set first species
      Yin[1] = 1.0 - m_YF; // set second species

      this->solveTPY(this->T, this->P, Yin);
      double rhoYF = this->rho * m_YF;
      double res = (rhoYF - RYF) / RYF; // relative error
      return res;
    } else { // no flash, we use pressure as the target
      Y_[0] = m_YF; // set first species
      Y_[1] = 1.0 - m_YF; // set second species
                          //
      double rho_guess = RYF / m_YF;
      this->SetMixture_TRY(this->T, rho_guess, Y_);
      double res = (P_target - this->P) / P_target; // relative error
      return res;
    }
  };

  double step = 1e-4;
  // initial guess
  // bool isRising = (resFun(YFGuess*(1+step))-resFun(YFGuess)) > 0;
  bool isRising = true;
  // tolerance
  boost::math::tools::eps_tolerance<double> tol(24);
  const double tolerance = 1e-8; // 1e-2 is a good value for RYF, but can be adjusted if needed
  auto residual_termination = [&](double a, double b) {
    double fa = std::fabs(resFun(a));
    double fb = std::fabs(resFun(b));
    bool res = (fa < tolerance) || (fb < tolerance);
    return res;
  };
  // Max iterations
  const boost::uintmax_t maxit = 500;
  boost::uintmax_t it = maxit;
  double factor = 1.10;

  std::pair<double, double> sol = boost::math::tools::bracket_and_solve_root(resFun, YFGuess, factor, isRising, residual_termination, it);


  if (it >= maxit) {
    std::cerr << "Error: SetMixture_TPRYF: Max iterations reached. T = " << T << ", P = " << P << ", RYF = " << RYF << std::endl;
    throw std::runtime_error("Max iterations reached in SetMixture_TPRYF");
  }

  double YF_sol = 0.5 * (sol.first + sol.second);
  if (std::fabs(resFun(YF_sol)) > tolerance || YF_sol > 1.0 || YF_sol < 0.0) { // sometimes this needs to be adjusted.
    YF_sol = sol.first;
    double res_first = resFun(sol.first);
    double res_second = resFun(sol.second);
    if (std::fabs(res_first) < std::fabs(res_second)) {
      YF_sol = sol.first;
    } else {
      YF_sol = sol.second;
      if (YF_sol > 1.0 || YF_sol < 0.0) {
        YF_sol = sol.first;
      }
    }
  }

  Y_[0] = YF_sol; // set first species
  Y_[1] = 1.0 - YF_sol; // set second species

  if (flash) {
    this->solveTPY(this->T, this->P, Y_); // solve for T and P
  } else {
    this->rho = RYF / YF_sol; // set density
    this->SetMixture_TRY(this->T, this->rho, Y_);
  }

  if (std::fabs((this->rho * Y_[0] - RYF)/RYF) > tolerance) {
    printf("Error: SetMixture_TPRYF: Final RYF = %f, expected RYF = %f, using YF = %f\n", this->rho * Y_[0], RYF, YF_sol);
    printf("sol.first = %f, sol.second = %f, sol.middle=%f\n", sol.first, sol.second,0.5*(sol.first + sol.second));
    throw std::runtime_error("Final RYF does not match expected value in SetMixture_TPRYF");
  }
} // SetMixture_TPRYF

// ----------------------------------------------------------------------------
// @brief: Temperature, pressure and mass fraction; typically used in initial conditions
// @param T: temperature [K]
// @param P: pressure [Pa]
// @param Y: mass fractions of species [kg/kg]
void PengRobinson::SetMixture_TPY(const double T, const double P, const double* Y) {
  this->solveTPY(T, P, Y, true);
} // end of PengRobinson::SetMixture_TPY
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @brief: (mass-based) energy, density, and species; used after conservative variables are set
// @param e: energy per unit mass [J/kg]
// @param rho: density [kg/m^3]
// @param Y: mass fractions of species [kg/kg]
void PengRobinson::SetMixture_ERY(const double e, const double rho, const double* Y, const double T_guess, const double P_guess, const bool flash) {
  if (flash) {
    const double err_ERHO = 1e-4;
    const bool DeltaTmax = 10.0;
    const bool DeltaPmax = 1e5;
    const int Nmax = 2000;
    const double lambda_0 = 0.5;
    const bool accelerate = true; // use acceleration for the flash
    this->solveERHOY(e, rho, Y, this->T, this->P, this->betaV, this->x_TP, this->y_TP,
                      err_ERHO, DeltaTmax, DeltaPmax, Nmax, 
                      T_guess, P_guess, 
                      lambda_0, accelerate
                     );
  } else { // no flash
    throw std::runtime_error("PengRobinson::SetMixture_ERY: no-flash is not implemented for this method.");
  }
} // end of PengRobinson::SetMixture_ERY
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// @brief: pressure, density, and species; used after using quasi-conservative methods such as double-flux
// @param P: pressure [Pa]
// @param rho: density [kg/m^3]
// @param Y: mass fractions of species [kg/kg]
void PengRobinson::SetMixture_PRY(const double P, const double rho, const double* Y, const double T_guess, const bool flash) {
  if (flash) {
    this->solvePRHOY(P, rho, Y, this->T, this->betaV, this->x_TP, this->y_TP,
        1e-8, 1.0, 1000, T_guess);
  } else {
    throw std::runtime_error("PengRobinson::SetMixture_PRY: no-flash is not implemented for this method.");
  }
} //end of PengRobinson::SetMixture_PRY
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// \brief: Solve the Rachford-Rice equation for the Peng-Robinson EOS (2.2.37)
// ----------------------------------------------------------------------------
double PengRobinson::solveRachfordRice(const double* X_in, const double* x_in, const double* y_in) {
  // Ensure x and y are normalized
  double sum_x = 0.0;
  double sum_y = 0.0;

  double x [this->num_species];
  double y [this->num_species];
  double K [this->num_species];

  LOOP_k_N(this->num_species) {
    sum_x += x_in[k];
    sum_y += y_in[k];
  }

  if (sum_x == 0.0) return 1.0; // pure vapor
  if (sum_y == 0.0) return 0.0; // pure liquid

  LOOP_k_N(this->num_species) {
    x[k] = x_in[k] / sum_x; // normalize x
    y[k] = y_in[k] / sum_y; // normalize y
    K[k] = y[k] / x[k];
  }

  // Check if single phase
  double sum_XK = 0.;
  double sum_XoK = 0.;

  LOOP_k_N (this->num_species) {
    sum_XK += X_in[k] * K[k];
    sum_XoK += X_in[k] / K[k];
  }

  if ((sum_XK - 1. < 0.0) && (std::fabs(sum_XK - 1.) > 1e-15)) {
    // this mixture is pure liquid
    return 0.0; // exit
  } else if ((1. - sum_XoK > 1e-15) && (true)) {
    // hits mixture is pure vapor
    return 1.0;
  }

  // Returns tuple: F(beta), F'(beta), F''(beta)
  auto rr = [&](double beta) {
    double F   = 0.0;
    double dF  = 0.0;
    double d2F = 0.0;

    for (size_t i = 0; i < this->num_species; ++i) {
      double a = K[i] - 1.0;
      double denom = 1.0 + beta * a;
      // safeguard against division by zero
      if (denom <= 0.0) denom = std::numeric_limits<double>::min();

      double inv = 1.0 / denom;
      double zi  = X_in[i];

      F   += zi * a * inv;
      dF  += -zi * a * a * inv * inv;
      d2F +=  2.0 * zi * a * a * a * inv * inv * inv;
    }
    return std::make_tuple(F, dF, d2F);
  };
  // Max iterations
  const boost::uintmax_t maxit = 500;
  boost::uintmax_t it = maxit;

  // std::pair<double, double> sol = boost::math::tools::toms748_solve(fun, 0., 1., tol, it);
  // beta_out = 0.5 * (sol.first + sol.second);
  // Initial guess (clamped)
  double num = 0.0, den = 0.0;
  for (size_t i = 0; i < this->num_species; ++i) {
    double Ki = std::max(K[i], 1e-30); // guard against div by zero
    num += X_in[i] * (1.0 - 1.0 / Ki);
    den += X_in[i] * ((K[i] - 1.0) / Ki);
  }
  double beta0 = (den != 0.0) ? (num / den) : 0.5;
  beta0 = std::max(0.0, std::min(1.0, beta0)); // clamp to [0, 1]

  // Bracket for halley_iterate
  double min_beta = 0.0;
  double max_beta = 1.0;
  int digits = std::numeric_limits<double>::digits; // full double precision
  double sol = boost::math::tools::halley_iterate(rr, beta0, min_beta, max_beta, digits, it);

  // if (it == maxit) {
  //   // if we reach the maximum number of iterations, we have to return
  //   std::cerr << "Error: PengRobinson::SSI: "
  //             << "Halley's method did not converge after " << maxit << " iterations." << std::endl;
  //   throw (-1);
  // }

  double beta_out = sol;
  SAFEGUARD(beta_out, "beta_out");
  return beta_out;
} // end of PengRobinson::solveRachfordRice
// -------------------------------------------------------------------------------
} // namespace Physics
