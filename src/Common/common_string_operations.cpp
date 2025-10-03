#include "common_string_operations.h"

#include <cstring> // for strtok
#include <cstdlib> // for atof
#include <sstream> // for istringstream

namespace common {

std::vector<double> Tokenize(const std::string stringIn,
    const char *delimiters) {
  // Member function ID string
  std::string memberID = "Tokenize";

  std::vector<double> outputVector;

  // flowControl while (for error handling)
  int flowControl = 1;
  while (flowControl == 1) {
    // Tokenize string read from file (spaces or commas are allowed(
    char *token = strtok((char *)stringIn.c_str(), delimiters);
    while (token != nullptr) {
      outputVector.push_back(std::atof(token));
      token = strtok(nullptr, delimiters);
    }

    // Continue flow normally (no unrecoverable errors up to this point)
    flowControl++;
  }

  return outputVector;
} // end Tokenize

std::vector<std::string> splitString(const std::string& str, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(str);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
} // end splitString

} // end namespace common

