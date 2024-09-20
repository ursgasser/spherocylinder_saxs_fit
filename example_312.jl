# fit SAXS data with many spherocylinders
# interpolation functions are used for the spherocylinder scattering amplitude to be fast

using HDF5, JLD
#using PyPlot, HDF5, JLD
include("spherocylinder_overlap.jl")
include("spherocylinder_amplitude_240825.jl")

sample_name = "cyl3_312_50k_"

# load data from file
d = load("cyl3_312_50k_0_000.jld")
q = d["q"];
Nq = length(q)
iq = d["iq"];
iqerr = d["iqerr"];
bg = d["bg"];
qh1 = d["qh1"];
qh2 = d["qh2"];
i15 = d["i15"];
pos = d["pos"];  # center-of-mass positions of cylinders
Npos = length(pos[:,1])
th = d["th"];
cth = cos.(th);
sth = sin.(th);
ph = d["ph"];
cph = cos.(ph);
sph = sin.(ph);
len = d["len"];
rad = d["rad"];
cylg = d["pcyl"];
thavg = d["thavg"];
cthavg = cos.(thavg);
sthavg = sin.(thavg);
phavg = d["phavg"];
cphavg = cos.(phavg);
sphavg = sin.(phavg);
Navg = length(thavg)
sig_hw = d["sig_hw"]  # sigma for hanning window
sig4_hw = 2.0 * sig_hw^2

# plot3D(pos[:,1], pos[:,2], pos[:,3], ",")

Nparts = length(cylg[1,:])
Nobj = length(cylg[:,1])
Nmov = 2000  # number of proteins in optimize group
Nopt = 4000  # change optimization group after this number of optimization steps
imax = 2000000  # number of optimization steps

ov = [sth.*sph sth.*cph cth];  # orientation of cylinders
mass_cyl = π*rad.^2 .* len;  # mass of cylinders
mass_cylg = sum(mass_cyl[cylg], dims=2);  # mass of cylinder groups
cm_cg = zeros(Float64, Nobj,3);  # center of mass of cylinder groups
ehw = zeros(Float64, Nobj);  # exponent for hanning window
for i = 1:Nobj
    cm_cg[i,:] = sum(mass_cyl[cylg[i,:]].*ones(Float64,1,3) .* pos[cylg[i,:],:], dims=1) ./
        mass_cylg[i]  # center of mass
    ehw[i] = -sum(cm_cg[i,:].^2) / sig4_hw  # exponent for hanning window
end
lengths = [155.0, 76.0, 50.0]
ilen = Int64.(round.(len[cylg] .- ones(Int64,Nobj) .* lengths')) .+ 1;  # index of interpol. function

qvavg = zeros(Float64, Navg, 3);  # q vectors for averaging
for i = 1:Navg
    qvavg[i,:] = [sthavg[i]*sphavg[i], sthavg[i]*cphavg[i], cthavg[i]]
end

# interpolation functions for spherocylinders with different length
amp_spherocyl = Array{Any,2}(undef, 3,6)
for j = 1:3
    L = lengths[j]
    for i = 1:6
        amp_spherocyl[j,i] = amplitude_spherocylinder_interpol(q, 12.0, L+i-1);
    end
end

# initial chi^2
alock = Threads.SpinLock()
F2q = zeros(Float64, Nq);  # F^2(q) result from model
# calculate F^2(q)
@time Threads.@threads for iavg = 1:Navg
    # angle between cylinder and q vector
    alpha = acos.(abs.( (sthavg[iavg] * sth) .* (cphavg[iavg] * cph .+ sphavg[iavg] * sph) .+
        cthavg[iavg] * cth ))
    qv = [sthavg[iavg]*sphavg[iavg], sthavg[iavg]*cphavg[iavg], cthavg[iavg]]
    aa = zeros(ComplexF64, Nq);  # amplitude of a configuration
    for icg = 1:Nobj
        #println(icg)
        aa += amplitude_spherocyl_group_interpol(q, qv, pos[cylg[icg,:],:], alpha[cylg[icg,:]],
            ilen[icg,:], ehw[icg])
    end
    lock(alock)
        F2q .+= abs.(aa).^2
    unlock(alock)
end
F2q ./= Navg;
# find scale factor
w = findall(x -> qh1 <= x <= qh2, q);
fact = i15 / sum(F2q[w]) * length(w)
chi20 = sum((iq .- fact*F2q .- bg).^2 ./ iqerr.^2) / Nq


r_max = [300.0, 300.0, 300.0];  # max movement
gr_start, gr_end = 1, Nmov
# array for the amplitudes of cylinders
ampl = zeros(ComplexF64, Nq, Navg, Nmov);  # the group to be optimized
ampl_r = zeros(ComplexF64, Nq, Navg);  # amplitude of the other, fixed objects
flag_reject = false

@time for iopt = 1:imax
    if mod(iopt, 1000) == 0
        print("iopt = ", iopt, "\n")
    end
    
    if mod(iopt, Nopt) == 1
        # change the set of objects to be optimized
        global gr_start = mod(gr_start-1+Nmov, Nobj)+1
        global gr_end = mod(gr_end-1+Nmov, Nobj)+1
        # calculate amplitudes for those to be optimized
        ampl_r .= 0.0  # reset amplitudes of cylinders that stay fixed
        Threads.@threads for iavg = 1:Navg
            # factors for q_parallel and q_perpendicular
            fpa = qvavg[iavg,1] .* ov[:,1] .+ qvavg[iavg,2] .* ov[:,2] .+ qvavg[iavg,3] .* ov[:,3]
            alpha = acos.(abs.(fpa))
            ##fpe = sqrt.(1.0 .- fpa.^2)
            qv = qvavg[iavg,:]
            # amplitude of all in the group to be optimized...
            for ig = gr_start:gr_end
                h = amplitude_spherocyl_group_interpol(q, qv, pos[cylg[ig,:],:], 
                    alpha[cylg[ig,:]], ilen[ig,:], ehw[ig])
                lock(alock)
                    ampl[:,iavg,ig-gr_start+1] = h
                unlock(alock)
            end
            # get amplitude of all not in the group...
            for ig = 1:Nobj
                if gr_start <= ig <= gr_end
                    continue
                end
                h = amplitude_spherocyl_group_interpol(q, qv, pos[cylg[ig,:],:], 
                    alpha[cylg[ig,:]], ilen[ig,:], ehw[ig])
                lock(alock)
                    ampl_r[:,iavg] += h
                unlock(alock)
            end
        end
        println("   changed group to ", gr_start, " ... ", gr_end)
    end  # change of group
    
    # choose particle in the optimize-group to change
    im = rand(range(1,Nmov))  # index of object and in ampl array
    ig = gr_start + im -1  # index of object in cylg array
    # choose kind of variation
    mode = rand(range(1,Nparts+1))  # 1: position, 2,3,4: orientation

    if mode == 1
        # change position
        ∆p = 2 * r_max .* (rand(3) .- 0.5)
        p0 = pos[cylg[ig,:],:]  # store old positions
        ampl0 = ampl[:,:,im]  # store old amplitude
        cm_cg0 = cm_cg[ig,:]  # store old center of mass
        ehw0 = ehw[ig]  # store old exponent for hanning window
        pos[cylg[ig,:],:] .+= ones(Float64,Nparts).*∆p'
    elseif 2 <= mode <= 4
        # change orientation of one cylinder
        ic = mode-1  # index of cylinder in object
        k = cylg[ig,ic]  # index in pos, ov arrays
        cthnew = rand(Float64)  # new cos(theta)
        sthnew = sqrt(1.0-cthnew^2)  # new sin(theta)
        sphnew, cphnew = sincos(2π*rand(Float64))  # new sin(phi), cos(phi)
        ov0 = ov[k,:]  # store old orientation
        ov[k,:] = [sthnew*sphnew, sthnew*cphnew, cthnew]
        p0 = pos[cylg[ig,:],:]  # store old positions

        # get new part positions
        # in the order 1,2,3, the cylinders are connected at p+l/2*ov of the previous cylinder
        if ic == 1
            # new pos for part 1
            pos[cylg[ig,ic],:] = spherocylinder_connect(pos[cylg[ig,2],:], -ov[cylg[ig,2],:], len[cylg[ig,2]],
                -ov[cylg[ig,1],:], len[cylg[ig,1]], 1.01*(rad[cylg[ig,2]]+rad[cylg[ig,2]]) )
        elseif ic == 2
            # new pos for parts 2 and 3
            for i = ic:Nparts
                pos[cylg[ig,i],:] = spherocylinder_connect(pos[cylg[ig,i-1],:], ov[cylg[ig,i-1],:], len[cylg[ig,i-1]],
                    ov[cylg[ig,i],:], len[cylg[ig,i]], 1.01*(rad[cylg[ig,i-1]]+rad[cylg[ig,i]]) )
            end
        elseif ic == 3
            # new pos for part 3
            pos[cylg[ig,ic],:] = spherocylinder_connect(pos[cylg[ig,2],:], ov[cylg[ig,2],:], len[cylg[ig,2]],
                ov[cylg[ig,3],:], len[cylg[ig,3]], 1.01*(rad[cylg[ig,2]]+rad[cylg[ig,3]]) )
        else
            println("ERROR: unknown case with ic=",ic)
        end        
        ampl0 = ampl[:,:,im]  # store old amplitude
        cm_cg0 = cm_cg[ig,:]  # store old center of mass
        ehw0 = ehw[ig]  # store old exponent for hanning window
    elseif mode > 4
        # ???
    end

    # check for collision with changed cylinders
    global flag_reject = false
    if mode == 1
        ich1, ich2 = 1, Nparts  # translation: check all
    elseif mode == 2
        ich1, ich2 = 1, 1
    elseif mode == 3
        ich1, ich2 = 2, Nparts
    elseif mode == 4
        ich1, ich2 = 3, Nparts
    end
    for ich = ich1:ich2
        if flag_reject
            break
        end
        Threads.@threads for i = 1:Npos
            if flag_reject
                break
            end
            if i == cylg[ig,ich]
                continue
            end
            h = spherocylinders_collision(pos[cylg[ig,ich],:], ov[cylg[ig,ich],:], len[cylg[ig,ich]],
                pos[i,:], ov[i,:], len[i], rad[cylg[ig,ich]]+rad[i])
            if h > 0.0
                lock(alock)
                    global flag_reject = true
                unlock(alock)
            end
        end
    end

    if !flag_reject
        # there are no collisions
        # recalculate amplitudes for changed cylinders...
        if mode == 1
            ##mass = π*rad[cylg[ig,:]].^2 .* len[cylg[ig,:]]  # weights of parts for center-of-mass calculation
            ##cm = sum(mass_cyl[cylg[ig,:]].*ones(Float64,1,3) .* p0, dims=1) ./ mass_cylg[ig]  # center of mass
            ##ehw0 = -sum(cm.^2) / sig4_hw  # exponent for hanning window
            cm_cg[ig,:] = sum(mass_cyl[cylg[ig,:]].*ones(Float64,1,3) .* pos[cylg[ig,:],:],
                dims=1) ./ mass_cylg[ig]  # center of mass
            ehw[ig] = -sum(cm_cg[ig,:].^2) / sig4_hw  # exponent for hanning window
            Threads.@threads for iavg = 1:Navg
                qv = qvavg[iavg,:]
                f = exp.((ehw[ig]-ehw0) .- Complex(0.0,1.0)*sum(qv.*∆p)*q)
                h = f .* ampl0[:,iavg]  # translation
                lock(alock)
                    ampl[:,iavg,im] = h
                unlock(alock)
            end
        elseif 2 <= mode <= 4
            # recalculate amplitude for the group
            cm_cg[ig,:] = sum(mass_cyl[cylg[ig,:]].*ones(Float64,1,3) .*
                pos[cylg[ig,:],:], dims=1) ./ mass_cylg[ig]  # center of mass
            ehw[ig] = -sum(cm_cg[ig,:].^2) / sig4_hw  # exponent for hanning window
            Threads.@threads for iavg = 1:Navg
                qv = qvavg[iavg,:]
                fpa = qv[1] * ov[cylg[ig,:],1] .+ qv[2] * ov[cylg[ig,:],2] .+
                    qv[3] * ov[cylg[ig,:],3]
                alpha = acos.(abs.(fpa))
                ##fpe = sqrt.(1.0 .- fpa.^2)
                h = amplitude_spherocyl_group_interpol(q, qv, pos[cylg[ig,:],:], 
                    alpha, ilen[ig,:], ehw[ig])
                lock(alock)
                    ampl[:,iavg,im] = h
                unlock(alock)
            end
        else
            println("ERROR: unknown mode:", mode)
        end

        # calculate F^2(q)
        F2q .= 0.0
        Threads.@threads for iavg = 1:Navg
            h = abs.(sum(ampl[:,iavg,:], dims=2) .+ ampl_r[:,iavg]).^2
            lock(alock)
                F2q .+= h
            unlock(alock)
        end
        F2q ./= Navg
        # find scale factor
        global w = findall(x -> qh1 <= x <= qh2, q);
        global fact = i15 / sum(F2q[w]) * length(w)
        global chi2 = sum((iq .- fact*F2q .- bg).^2 ./ iqerr.^2) / Nq
        # accept or reject the change
        if chi2 < chi20
            print("iopt=", iopt, "  changing obj ", ig, "  chi2: ", chi2,
                "  mode: ", mode, "\n")
            global chi20 = chi2
        else
            global flag_reject = true
        end
    end  # flag_reject is now set

    if flag_reject
        if mode == 1
            pos[cylg[ig,:],:] = p0
            cm_cg[ig,:] = cm_cg0
            ehw[ig] = ehw0
            ampl[:,:,im] = ampl0
        elseif 2 <= mode <= 4
            pos[cylg[ig,:],:] = p0
            cm_cg[ig,:] = cm_cg0
            ehw[ig] = ehw0
            ov[k,:] = ov0
            ampl[:,:,im] = ampl0
        else
            println("ERROR: unknown mode:", mode)
        end
    end

    if mod(iopt, 100) == 0
        print("iopt=", iopt, "  chi20=", chi20, "  chi2=", chi2, "\n")
        GC.gc()  # prevent excessive memory use
    end

    if mod(iopt, 10000) == 0
        # save fit result
        global cth = ov[:,3];  # cos(theta)
        global th = acos.(cth);
        #sth = sqrt.(1.0 .- cth.^2);  # sin(theta)
        global ph = atan.(ov[:,1],ov[:,2]);  # phi
        #cph = cos.(ph);  # cos(phi)
        #sph = sin.(ph);  # sin(phi)
        s1,s2 = split(string(chi20),".")
        out_name = sample_name * s1*"_"*s2[1:3] * ".jld"
        save(out_name, Dict("q"=> q, "iq"=> iq, "iqerr"=> iqerr, "bg"=> bg,
            "qh1"=> qh1, "qh2"=> qh2, "i15"=> i15,
            "pos"=> pos, "th"=> th, "ph"=> ph, "len"=> len, "rad"=> rad,
            "pcyl"=> cylg, "thavg"=> thavg, "phavg"=> phavg, "sig_hw"=> sig_hw))
    end
end



### plot the fit ###
#=
using PyPlot
loglog(q,iq,".")
w = findall(x -> x>=qh1 && x<=qh2, q);
f = i15 / sum(F2q[w]) * length(w)
loglog(q, f*F2q .+ bg)
=#
