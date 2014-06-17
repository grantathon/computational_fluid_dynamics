computational_fluid_dynamics
============================

Assignments for Scientific Computing Lab for CFD (IN 2186) at TU Munich SS14

Assignments
============================

1. Navier-Stokes Equations
2. Lattice Boltzmann Method
3. Arbitrary Geometries with Navier-Stokes Equations
4. Navier-Stokes Equations using MPI

Project
============================

The goal of this project is to use a Monte Carlo method and the Navier-Stokes equations to quantify the uncertainty of the effect of the Reynolds number on the separation point following a step within a 2D grid.  To approach maximum parallelization, the Monte Carlo method will launch multiple programs over multiple nodes using MPI while the individual programs being launched will use a hybrid of MPI and OpenMP.
