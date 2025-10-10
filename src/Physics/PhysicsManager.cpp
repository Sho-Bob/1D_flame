#include "PhysicsManager.h"

#include <iostream>

Physics::Physics() {
  this->eos = new PengRobinson();
}

void Physics::initialize() {
  this->species_names = this->input->getStringArrayParam("species");
  int num_species = this->species_names.size();
  std::cout << "Number of species: " << num_species << std::endl;
  for (const auto& name : this->species_names) {
    std::cout << "Species: " << name << std::endl;
  }

  this->eos->init(this->species_names, 1e5);

  // set pointers
  this->rho = this->eos->GetRho_ptr();
  this->T = this->eos->GetT_ptr();
  this->p = this->eos->GetP_ptr();
  this->Y = this->eos->GetY_ptr();
  this->e = this->eos->GetE_ptr();

  this->sos = this->eos->getSos_ptr();
}

Physics::~Physics() {
  delete eos;
}
