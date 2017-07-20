#ifndef _INITLB_H_
#define _INITLB_H_

#include "LBDefinitions.hpp"
#include "LBMHelper.hpp"
#include <stdbool.h>
#include <vector>

/* reads the parameters for the lid driven cavity scenario from a config file */
void readParameters(
    coord_t &length,                /* reads domain size. Parameter name: "length" */
    coord_t &procs,                 /* reads number of processes per dimension. Parameter name
                                       "xprocs, yprocs, zprocs" */
    double &tau,                    /* relaxation parameter tau. Parameter name: "tau" */
    boundary_t &boundaryConditions, /* contains information about all bcs" */
    int &timesteps,                 /* number of timesteps. Parameter name: "timesteps" */
    int &timestepsPerPlotting,      /* timesteps between subsequent VTK plots.
                                   Parameter name: "vtkoutput" */
    std::string &scenario,          /* name of the geometry file. Parameter name: "scenario" */
    char *parameterFile, bool verbose);

/* initialises the particle distribution functions and the flagfield */
void initialiseCollideAndStreamFields(std::vector<double> &collideField,
                                      std::vector<double> &streamField);

void initialiseFlagField(std::vector<flag_t> &flagField, const char *geometryFile,
                         const coord_t &length, const coord_t &offset, const coord_t coord,
                         const coord_t procs, boundary_t &boundaryConditions, bool verbose);

std::vector<double> initialiseMassField(std::vector<flag_t> &flags, const coord_t &length);

void initialiseInterface(std::vector<double> &distributions, std::vector<double> &mass,
                         std::vector<flag_t> &flags, const coord_t &length);

#endif
