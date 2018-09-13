# Timothy - Tissue Modelling Framework

Timothy is a novel large scale parallel computational model allowing 3-D simulations of cell colonies growing and interacting with variable environment in previously unavailable tissue scale.

The cells are modelled as individuals located in the lattice-free 3-D space. The model incorporates cellular environment modelled in a continuous manner, mathematical description based on partial differential equations is formulated for selected important components of the environment. Discrete and continuous formulations are efficiently coupled in one model and allow considerations on different scales: sub-cellular, cellular and tissue scale.

High parallel scalability achieved allows simulation of up to 10^9 individual cells. This large scale computational approach allows for simulations to be carried out over realistic spatial scales up to 1cm in size i.e. the tissue scale.

# Implementation

Timothy is written in C and makes use of Message Passing Interface (MPI) library for internode communication and OpenMP for intranode communication. There are few prerequisites needed for the successful installation of Timothy:
* C compiler supporting OpenMP (tested compilers: GNU gcc, IBM xlc),
* MPI library (tested implementations: OpenMPI, IBM POE for Power7 AIX, IBM MPI library for Blue Gene/Q),
* Zoltan library (available at http://www.cs.sandia.gov/zoltan/),
* Hypre library (available at http://acts.nersc.gov/hypre/).
