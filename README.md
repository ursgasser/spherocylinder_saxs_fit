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

File `example_312.jl` contains the Julia code to fit example data given in file `cyl3_312_50k_0_000.jld`.
This `jld` file also contains the spherocylinder model and other data such as the q-vector orientations 
used in the fit.
To avoid the re-calculation of all scattering amplitudes after a change of the model, the scattering
amplitudes of the spherocylinder triplets are kept in memory and are re-calculated only after a change
of the model is accepted. A translation or re-orientation of a triplet is rejected, if it does not
improve the fit to the SAXS data.
After 10000 optimization steps, the improved model is saved in the `jld` file.

Run the example from the command line with `> julia --threads 8 example_312.jl`
