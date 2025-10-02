#include "input.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

Input::Input(const std::string filename) : filename(filename) {
  this->config = toml::parse_file(filename); // create the toml table
}

std::string Input::getStringParam(const std::string& param_name) {
  std::string line;
  try {
    line = this->config[param_name].value_or("");
  } catch (const std::exception& e) {
    std::cerr << "Error reading parameter '" << param_name << "': " << e.what() << std::endl;
    line = "";
  }
  return line;
}

double Input::getDoubleParam(const std::string& param_name) {
  double value = 0.0;
  try {
    value = this->config[param_name].value_or(0.0);
  } catch (const std::exception& e) {
    std::cerr << "Error reading parameter '" << param_name << "': " << e.what() << std::endl;
    value = 0.0;
  }
  return value;
}


int Input::getIntParam(const std::string& param_name) {
  int value = 0;
  try {
    value = this->config[param_name].value_or(0);
  } catch (const std::exception& e) {
    std::cerr << "Error reading parameter '" << param_name << "': " << e.what() << std::endl;
    value = 0;
  }
  return value;
}
