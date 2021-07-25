# grav_cube_inversion
Grav-parallel is a C code written in parallel (with MPI) designed to model the gravity anomaly due to a body that can be represented by prisms. The code assumes that the prisms have a uniform top depth and uniform density contrast. The code models the depth to the bottom of each prism.

The gbox forward model is used. The inversion is done using the Ameoba algorthim, also called the Nedler-Meade simplex method.
