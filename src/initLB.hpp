#ifndef _INITLB_H_
#define _INITLB_H_

#include "LBDefinitions.hpp"
#include "LBMHelper.hpp"
#include "scenarios/Scenario.hpp"
#include <bits/unique_ptr.h>
#include <stdbool.h>
#include <vector>

/* reads the parameters for the simulation from a config file */
void readParameters(coord_t &length, double &tau, double &smagorinskyConstant,
                    boundary_t &boundaryConditions, int &timesteps, int &timestepsPerPlotting,
                    std::unique_ptr<Scenario> &scenario, char *parameterFile, bool verbose,
                    std::array<double, 3> &gravity);

/* initialises the particle distribution functions and the flagfield */
void initialiseCollideAndStreamFields(std::vector<double> &collideField,
                                      std::vector<double> &streamField);

void initialiseFlagField(std::vector<flag_t> &flagField, std::unique_ptr<Scenario> scenario,
                         boundary_t &boundaryConditions, bool verbose, const coord_t &length);

std::vector<double> initialiseMassField(std::vector<flag_t> &flags, const coord_t &length);

void initialiseInterface(std::vector<double> &distributions, std::vector<double> &mass,
                         std::vector<double> &density, const coord_t &length,
                         std::vector<flag_t> &flags);

#endif
