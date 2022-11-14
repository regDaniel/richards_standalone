# REQ-MODULAR

## Purpose

This package is written in order to investigate the numerical convergence or Richards Equation (Richards, 1931) and to test different parameterization options for soil hydraulic functions such as the Mualem (1976) - van Genuchten (1980) or the Rijtema (1969) model.

## Design Rules

This software was deliberately written in a modular fashion and with at least some ideas from functional programming. The purists may pardon me if at some places I am not very strict in this approach. The idea is to have a package that can easily be understood, maintained and extended by scientific programmers, PhD and MSc students alike.

The main idea is to have functions that are called from a main program, somewhat similar to the model used in the Numerical Modeling for Weather and Climate course. However, the approach here is more strict, as we try to avoid if-clauses concerning parameterization options in the individual functions. Instead, each option should get an individual function. This way users and maintainers should never have to read and decipher more than approximately 50 lines of code at the same time.

One more consideration concerns the programming paradigm. We chose a functional approach for two reasons: 1. It is natural to a program with no interaction by users during run time. 2. It doesn't require users and maintainers to be familiar with object-oriented programming and how OO programming is working in Python. This is particularly important as the target audience of the package is domain scientists and students such as myself.

To avoid setting global variables as options or to import large numbers of parameters in each file header, we make heavy use of dictionaries with groups of parameters concerning similar subjects. These parameters can then be unpacked wherever they are required, they should remain unaltered and it is easy to trace them throughout the whole code.

## File Hierarchy

### Main Folder

The main folder contains the python modules that run the model (the numerical implementation of Richards Equation) and some plotting routines. All files in the main folder are subject to a strict documentation policy which will be explained later in this document. If you violate the documentation policy, you will spend your afterlife trying to decipher Fortran code in the cosmo model.

### analysis

This folder contains tools for convergence analyses. The documentation policy is not as strict as for the main folder.

### run_scripts

This folder contains run scripts for experiments requiring multiple singe or multi-column runs with slightly different setups. This is for instance the case for parameter sensitivity or convergence experiments. The documentation policy is not as strict as for the main folder.

## Files in the main folder

### run_column.py

This is the main program for single column runs. It contains setup, initialization and time loop and makes the corresponding calls to functions of the other modules. It is important to understand the following design rule: All if statements concerning parameterization options are contained in the main programs (i.e. here and maybe in some lateral flow program in the future). Depending on the case, the main program will then call the appropriate functions. The author has a strong feeling that this way the code will be easier to understand as there will be no need to trace if-statements through all functions. The drawback of this approach is that some functions might look quite similar to each other if one does not have the time and patience to generalize the functions appropriately, which is almost certainly the case in this context. However the author feels that the readability argument is stronger in our context.

### solver.py

This package contains numerical solvers that are called by the main program on each time level. The implicit solvers depend on the scipy.linalg package. For reasonable setups, this is the most expensive part of the model, so if you are concerned with efficiency/ compute time/ gpuing things, this is probably the place to implement changes.

### fluxes.py

This package contains the soil hydraulic functions for different parameterization options and the estimation of ground runoff after Schlemmer et al. (2018).


### grid.py

Contains functions that set up the grid and that are required to calculate and apply metric terms in case the user chooses to run the model with a 'telescope-like' grid (see Regenass et al. (2021) for a description).

### init.py

Initialization of arrays.

### output.py

Functions that manage the storage of output. The output is written to netCDF files with (t, z) dimensions.

### parameter.py

This is the namelist. We tried to move away from setting global variables, because it can cause bugs that are hard to find and fix. Here, the parameter values are stored in dictionaries and may then be retrieved in the functions from these dictionaries. The performance penalty is marginal and the gain in safety and readability is probably worth it.
