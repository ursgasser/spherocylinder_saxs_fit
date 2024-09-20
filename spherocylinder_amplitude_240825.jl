# spherocylinder scattering amplitude using interpolation

using Interpolations, SpecialFunctions

"""
this returns an interpolation function for the scattering amplitude of a spherocylinder.
"""
function amplitude_spherocylinder_interpol(qs::Array{Float64,1}, rad::Float64, len::Float64; Nz=500, Nalpha=10000)
    #   alpha is the argument of interpolation function and qs, rad, and len are fixed.
    #   usage: > spherocyl = amplitude_spherocylinder_interpol(q, R, L)
    #          > alpha = 0.734
    #          > fq = spherocol(alpha)  # the scattering amplitude of the spherocylinder
    # alpha: angle between cylinder axis and q vector (0 ... π/2)
    # qs: q values to calculate the scattering amplitude
    # rad: cylinder radius
    # len: length of cylinder (without spherical caps)
    
    # use rad as unit length
    q = qs * rad
    L = len/rad/2.0  # half the cylinder length
    z = [x for x in range(0.0, 1.0, Nz)];
    emz2 = 1.0 .- z.^2;
    s1mz2 = sqrt.(emz2);

    als = range(0.0, π/2, Nalpha)
    fqal = Vector{Vector{Float64}}(undef, length(als));
    alock = Threads.SpinLock()
    Threads.@threads for ial = 1:length(als)
        local al = als[ial]  # alpha: the angle between cylinder axis and q vector
        # jinc is the function 2*J_1(π*x) / (π*x)
        # scattering amplitude for this alpha...
        local h = (2π/(Nz-1)) * sum( ones(Float64,length(q)) .* emz2' .* 
            jinc.((sin(al)/π) * q .* s1mz2') .* cos.(cos(al) * q * (L.+z)'), dims=2 ) .+
            (2π*L) * sinc.(cos(al) * q * (L/π)) .* jinc.(sin(al)/π * q);  # add cylinder
        lock(alock)
            fqal[ial] = dropdims(h,dims=2)
        unlock(alock)
    end
    fqal *= rad^3  # switch back to length unit used in input
    return( linear_interpolation(als, fqal) )  # interpolation of vectors
end

#==============
# code to test the function given above

using PyPlot
include("/Users/ursgasser/varia/julia/spherocylinder_amplitude_240825.jl")
include("/Users/ursgasser/varia/julia/formfactors.jl")

R = 1.8
L = 0. #100.

qmin, qmax = 0.01, 10.0
q = exp.( range(log(qmin), log(qmax), 200) );

amplitude_spherocylinder = amplitude_spherocylinder_interpol(q, R, L);

al = π/2 #0.0 #π/2
loglog(q, amplitude_spherocylinder(al).^2)
plot(qmin*[1,10], (π*R^2*L + 4π/3*R^3)^2 * [1,1])
qpa = cos(al) * q;
qpe = sin(al) * q;
loglog(q, amplitude_cylinder(qpa, qpe, L, R).^2)
plot(qmin*[1,10], (π*R^2*L)^2 * [1,1])

plot(q, (4π/3*R^3)^2 * formfactor_sphere(q, R))
===============#


"""
calculate scattering amplitude for a group of spherocylinders with given positions (pos), lengths (ilen)
and angles relative to the q-vector (alpha).
"""
function amplitude_spherocyl_group_interpol(q::Array{Float64,1}, qv::Array{Float64,1}, pos::Array{Float64,2}, alpha::Array{Float64,1},
    ilen::Array{Int,1}, ehw::Float64)
    # q: q values
    # qv: unit vector giving q direction
    # pos: positions of cylinders in the group
    # alpha: angles between cylinders and q vector
    # ilen: index of interpolation function (for the length of the spherocylinder)
    # ehw: exponent for hanning window
    ampl = zeros(ComplexF64, length(q));  # scattering amplitude
    for i = 1:length(alpha)
        ampl += exp.(ehw .- Complex(0.0,1.0)*q*sum(qv.*pos[i,:])) .*
            amp_spherocyl[i,ilen[i]](alpha[i])
    end
    return(ampl)
end
