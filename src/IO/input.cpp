#include "input.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

Input::Input(const std::string filename) : filename(filename) {
  this->config = toml::parse_file(filename); // create the toml table
}

std::vector<std::string> Input::getStringArrayParam(const std::string& param_name) {
  std::vector<std::string> result;
  try {
    auto arr = this->config[param_name].as_array();
    if (arr) {
      for (const auto& item : *arr) {
        if (item.is_string()) {
          result.push_back(item.value_or(""));
        } else {
          std::cerr << "Warning: Non-string item found in array for parameter '" << param_name << "'" << std::endl;
        }
      }
    } else {
      std::cerr << "Warning: Parameter '" << param_name << "' is not an array." << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error reading parameter '" << param_name << "': " << e.what() << std::endl;
  }
  return result;
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

std::string Input::getStringParam(const std::string& param_name, const std::string& default_value) {
  std::string line;
  try {
    line = this->config[param_name].value_or(default_value);
  } catch (const std::exception& e) {
    line = default_value;
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

double Input::getDoubleParam(const std::string& param_name, double default_value) {
  double value = default_value;
  try {
    value = this->config[param_name].value_or(default_value);
  } catch (const std::exception& e) {
    value = default_value;
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

int Input::getIntParam(const std::string& param_name, int default_value) {
  int value = default_value;
  try {
    value = this->config[param_name].value_or(default_value);
  } catch (const std::exception& e) {
    value = default_value;
  }
  return value;
}
