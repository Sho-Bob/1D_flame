#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <vector>
#include <toml++/toml.h>  // Ensure toml++ headers are in your include path

class Input {
  public:
    Input(const std::string filename);
    std::string getStringParam(const std::string& paramName);
    std::string getStringParam(const std::string& paramName, const std::string& defaultValue);

    std::vector<std::string> getStringArrayParam(const std::string& paramName);

    double getDoubleParam(const std::string& paramName);
    double getDoubleParam(const std::string& paramName, double defaultValue);

    int getIntParam(const std::string& paramName);
    int getIntParam(const std::string& paramName, int defaultValue);
  private:
    toml::table config;
    std::string filename;
};

#endif // INPUT_H

