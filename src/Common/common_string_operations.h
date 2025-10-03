#ifndef COMMON_STRING_OPERATIONS_H
#define COMMON_STRING_OPERATIONS_H

#include <string>
#include <vector>


namespace common {

std::vector<double> Tokenize(const std::string stringIn, const char *delimiters);
std::vector<std::string> splitString(const std::string& str, char delimiter = ',');

} // end namespace common

#endif // COMMON_STRING_OPERATIONS_H
