# Compilation

Compile using `cmake . ` followed by `make`.
Cmake should find the MPI-Libraries automatically.

# Usage
Use the supplied params.dat (=Lid driven cavity).
To do so, run `mpirun -n 8 lbsim params.dat`.
The examples in the examples folder are not supported for the MPI-Version.

The results can be found in the folder results/.
