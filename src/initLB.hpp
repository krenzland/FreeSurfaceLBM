#ifndef _INITLB_H_
#define _INITLB_H_

#include "LBDefinitions.hpp"
#include "LBMHelper.hpp"
#include "scenarios/Scenario.hpp"
#include <stdbool.h>
#include <vector>
#include <bits/unique_ptr.h>

/* reads the parameters for the lid driven cavity scenario from a config file */
void readParameters(
        coord_t &length,                /* reads domain size. Parameter name: "length" */
        double &tau,                    /* relaxation parameter tau. Parameter name: "tau" */
        std::array<double, 3> &gravity, /* force acting on all cells */
        boundary_t &boundaryConditions, /* contains information about all bcs" */
        int &timesteps,                 /* number of timesteps. Parameter name: "timesteps" */
        int &timestepsPerPlotting,      /* timesteps between subsequent VTK plots.
                                   Parameter name: "vtkoutput" */
        std::unique_ptr<Scenario> &scenario,          /* name of the scenario. Parameter name: "scenario" */
        char *parameterFile, bool verbose);

/* initialises the particle distribution functions and the flagfield */
void initialiseCollideAndStreamFields(std::vector<double> &collideField,
                                      std::vector<double> &streamField);

void initialiseFlagField(std::vector<flag_t> &flagField, std::unique_ptr<Scenario> scenario, boundary_t &boundaryConditions,
                         bool verbose, const coord_t &length);

std::vector<double> initialiseMassField(std::vector<flag_t> &flags, const coord_t &length);

void initialiseInterface(std::vector<double> &distributions, std::vector<double> &mass,
                         std::vector<double> &density, const coord_t &length,
                         std::vector<flag_t> &flags);

#endif
