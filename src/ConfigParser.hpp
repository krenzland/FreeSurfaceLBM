#ifndef CFD_LAB_CONFIGPARSER_HPP
#define CFD_LAB_CONFIGPARSER_HPP

#include "LBDefinitions.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>

class ConfigParser {
    using rawConfig_t = std::unordered_map<std::string, std::string>;

  public:
    ConfigParser(const std::string &filename, bool isVerbose) : isVerbose(isVerbose) {
        rawConfig = parseConfig(filename);
    }

    template <typename T> T parse(const std::string &paramName) {
        // We delegate this to the implementation struct because we need to
        // specialise the parse function
        // for the boundary flags and C++ doesn't allow this.
        return parseImpl<T>::parse(*this, paramName);
    }

  private:
    rawConfig_t parseConfig(const std::string &filename) {
        rawConfig_t map = rawConfig_t();
        auto ifs = std::ifstream(filename);
        if (!ifs) {
            throw std::logic_error("Failed to open config file " + filename + ".");
        }

        // First match is variable name, second match is value.
        auto format = std::regex(R"(^(\S*)\s*(\S*))", std::regex::ECMAScript);
        std::string line;
        while (std::getline(ifs, line)) {
            // First remove comments
            line = line.substr(0, line.find('#'));

            std::smatch match;
            std::regex_search(line, match, format);
            if (match.size() == 3) {
                if (match[0].str().size() == 0)
                    continue; // Skip empty lines
                map[match[1].str()] = match[2].str();
            }
        }
        return map;
    }

    template <typename R, typename = void> struct parseImpl {
        static R parse(ConfigParser &config, const std::string &paramName) {
            const auto it = config.rawConfig.find(paramName);
            if (it != config.rawConfig.end()) {
                std::istringstream ss(it->second);
                R value;
                ss >> value;
                if (config.isVerbose) {
                    std::cout << "config: " << std::setw(20) << paramName << ":" << std::setw(25)
                              << value << std::endl;
                }
                return value;
            } else {
                throw std::runtime_error("Parameter " + paramName + " not found or wrong format.");
            }
        }
    };

    template <typename S> struct parseImpl<flag_t, S> {
        static flag_t parse(ConfigParser &config, const std::string &paramName) {
            const auto configIt = config.rawConfig.find(paramName);
            if (configIt != config.rawConfig.end()) {
                const auto &value = configIt->second;
                const auto flagIt = stringToFlag.find(value);
                if (flagIt != stringToFlag.end()) {
                    const flag_t parsed = flagIt->second;
                    if (config.isVerbose) {
                        std::cout << "config: " << std::setw(20) << paramName << ":"
                                  << std::setw(25) << value << "( =  " << static_cast<int>(parsed)
                                  << " )" << std::endl;
                    }
                    return parsed;
                }
                {
                    throw std::runtime_error(value + " is not a valid flag field for boundary " +
                                             paramName + ".");
                }
            } else {
                throw std::runtime_error("Parameter " + paramName + " not found or wrong format.");
            }
        }
    };

    rawConfig_t rawConfig;
    const bool isVerbose;
};

#endif // CFD_LAB_CONFIGPARSER_HPP
