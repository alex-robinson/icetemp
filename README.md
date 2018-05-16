This repository is intended to serve as a sandbox for 
ice sheet thermodynamics solver development.

The idea is to perform 1D test cases that can be
compared to analytical or known solutions to ensure
that the thermodynamics solvers are:

    1. Validated
    2. Portable

Once they are confirmed to be working, since they
must be designed as stand-alone solvers, it will
be easy to plug them into any ice sheet model.

Getting started
#######

    1. Clone the repository 
    2. Modify the Makefile to match your systems settings.
       Dependencies: NetCDF library.
    3. Compile the code: `make test_icetemp`.
    4. Run the executable: `./test_icetemp.x`.
    5. Check the results.


