# spherocylinder_saxs_fit

A set of spherocylinders is used to fit SAXS data of Bik1 protein droplets.

Under certain conditions, Bik1 dimers do not stay in suspension but form concentrated protein droplets
that consist of a fractal network of Bik1 dimers.
A Bik1 dimer is modeled using three connected spherocylinders that can be freely oriented relative to
each other. The fractal network is formed by many Bik1 dimers. A large number of spherocylinders is 
used to represent a 'model Bik1 droplet'.

During the fit, the spherocylinders are translated and rotated to imporove the fit to the SAXS data. It 
is always checked that the spherocylinders do not overlap and that each triplet of spherocylinders 
representing a Bik1 dimer stays connected.

As there is no analytical expression for the scattering amplitude of a spherocylinder, the amplitudes of
the set of spherocylinders used in the model is calculated in the beginning to obtain interpolation 
functions that allow for a fast calculation of all scattering amplitudes.

File ??? contains the Julia code to fit example data given in file `cyl3_312_50k_0_000.jld`. This `jld`
file also contains the spherocylinder model and other data such as the q-vector orientations used in
the fit.
