using LinearAlgebra

"""
get minimum distance between two lines in 3d:
lines_dist(u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, b::Array{Float64,1})
    u: unit vector defining direction of line 1
    v: unit vector defining direction of line 2
    a: vector giving a point on line 1
    b: vector giving a point on line 2
"""
function lines_dist(u::Array{Float64,1}, v::Array{Float64,1},
    a::Array{Float64,1}, b::Array{Float64,1})
    d = b .- a;  # vector connecting both lines
    cv = cross(u,v);  # vector perpendicular to both lines

    ##if sum(cv.^2) == 0.0
    if abs(sum(u.*v)) == 1.0  # better to detect antiparallel line
        # parallel orientation
        d .-= sum(d.*u) * u
        return(sqrt(sum(d.^2)), true)
    else
        # orientation not parallel
        cv ./= sqrt(sum(cv.^2))
        return(abs(sum(d.*cv)), false)
    end
end


"""
get parllel and perpendicular distance between a line and a point
para_perp_dist(a::Array{Float64,1}, u::Array{Float64,1}, b::Array{Float64,1})
    # a: point on line
    # u: unit vector giving direction of line
    # b: point
"""
function para_perp_dist(a::Array{Float64,1}, u::Array{Float64,1}, b::Array{Float64,1})
    # a: point on line
    # u: direction unit vector of line
    # b: point
    d = b .- a;  # vector connecting point and line
    c1 = sum(d.*u)
    cv = cross(u,d);  # to get perpendicular component of d
    return( c1, sqrt(sum(cv.^2)) )
end


"""
check for overlap of two spherocyliners
it is assumed that the spherocylinders have the same diameter.
I THINK IT IS ASSUMED THAT BOTH SPEROCYLINDERS HAVE A LENGTH > DIAMETER

spherocylinders_collision(p1::Array{Float64,1}, u1::Array{Float64,1}, l1::Float64, p2::Array{Float64,1}, u2::Array{Float64,1}, l2::Float64, diam::Float64)
p1, p2: center-of-mass positions of the cylinders
u1, u2: orientation unit vectors of the cylinders
l1, l2: lengths of the cylinders
diam: diameter of cylinders
return value: 1.0: overlap, 0.0: no overlap
"""
function spherocylinders_collision(p1::Array{Float64,1}, u1::Array{Float64,1}, l1::Float64,
    p2::Array{Float64,1}, u2::Array{Float64,1}, l2::Float64, diam::Float64)
    # p1, p2: center-of-mass positions of the cylinders
    # u1, u2: orientation unit vectors of the cylinders
    # l1, l2: lengths of the cylinders
    # diam: diameter of cylinders

    hl1 = l1/2.0  # half lengths
    hl2 = l2/2.0
    # get min distance of the central lines of the cylinders
    dmin, flag_parallel = lines_dist(u1,u2, p1,p2)
    if dmin >= diam
        # no collision possible
        #println("no overlap! dmin=", dmin)
        return(0.0)
    end

    # central lines show that collision is possible

    if flag_parallel
        # parallel and min. distance smaller than diameter -> collision possible
        ∆ = sum((p2-p1).*u1)
        if abs(∆) < hl1+hl2
            # centers are close enough for collision
            return(1.0)
        elseif abs(∆) < hl1+hl2+diam
            # check collision of spherical caps
            e1 = p1 + sign(∆)*hl1*u1  # end of cyl1
            e2 = p2 + hl2*u2  # end of cyl2
            if sum((e2-e1).^2) < diam^2
                # collision
                return(1.0)
            end
            e2 = p2 - hl2*u2  # end of cyl2
            if sum((e2-e1).^2) < diam^2
                # collision
                return(1.0)
            end
            # no collision
            return(0.0)
        else
            return(0.0)
        end
    end

    d = p1 - p2
    if sum(d.^2) > (hl1+hl2+diam)^2
        # too far apart, no collision possible
        return(0.0)
    end

    # cylinders are close enough and central lines suggest a collision
    # cylinders are not parallel

    # calculate points at min distance
    # x1 = pos1 + t1*uv1; x2 = pos2+t2*uv2
    # d = pos2 .- pos1
    cv = cross(u1,u2)  # perpendicular to u1 and u2
    acv = sqrt(sum(cv.^2))
    ##cv ./= acv  # unit vector
    d = p2-p1
    t1 = sum(cross(d-dmin/acv*cv, u2) .* cv) / acv^2
    t2 = sum(cross(d-dmin/acv*cv, u1) .* cv) / acv^2
    if -hl1 < t1 < hl1 && -hl2 < t2 < hl2
        # both points at min distance are in cylinders
        #println("overlap! x1 and x2 are in: dmin=", dmin, "  ",t1, "  ",t2)
        return(1.0)
    elseif -hl1 < t1 < hl1
        # not both points at min distance are inside the cylinders
        # closest approach is in cyl1
        # get distance from end of cyl2 to cyl1
        for fsign in [1.0, -1.0]  # check both ends of cyl2
            ep = p2 + fsign*hl2*u2  # end point
            dpara, dperp = para_perp_dist(p1, u1, ep)
            if abs(dpara) < hl1
                # end point is next to other cylinder
                if dperp < diam
                    return(1.0)
                end
            else
                # end point ep is beyond end of other cylinder
                # calculate distance between end points
                d = ep - (p1 + sign(dpara)*hl1*u1)
                if sqrt(sum(d.^2)) < diam
                    return(1.0)
                end
            end
        end
        # no collision detected
        return(0.0)
    elseif -hl2 < t2 < hl2
        # not both points at min distance are inside the cylinders
        # closest approach is in cyl2
        # get distance from end of cyl1 to cyl2
        for fsign in [1.0, -1.0]  # check both ends of cyl1
            ep = p1 + fsign*hl1*u1  # end point
            dpara, dperp = para_perp_dist(p2, u2, ep)
            if abs(dpara) < hl2
                # end point is next to other cylinder
                if dperp < diam
                    return(1.0)
                end
            else
                # end point ep is beyond end of other cylinder
                # calculate distance between end points
                d = ep - (p2 + sign(dpara)*hl2*u2)
                if sqrt(sum(d.^2)) < diam
                    return(1.0)
                end
            end
        end
        # no collision detected
        return(0.0)
    end

    # min distance is outside both cylinders

    # calculate points where infinite cylinder hulls touch
    cos12 = sum(u1.*u2)  # angle between cylinders
    alpha = acosd(cos12)
    if alpha > 90.0
        ff = -1.0
    else
        ff = 1.0
    end
    # find points where surfaces of infinite cylinders touch
    # y1 = x1 + xi*uv1; y2 = x2 + xi*uv2
    # z1 = x1 - xi*uv1; z2 = x2 - xi*uv2
    xi = sqrt( (diam^2 - dmin^2) / (2.0*(1.0 - ff*cos12)) )
    if -hl1 < t1+xi < hl1 && -hl2 < t2+ff*xi < hl2
        # surfaces touch within cylinders at y1,y2: collision
        return(1.0)
    elseif -hl1 < t1-xi < hl1 && -hl2 < t2-ff*xi < hl2
        # infinite surfaces touch within cylinders at z1,z2: collision
        return(1.0)
    elseif -hl1 < t1+xi < hl1 || -hl1 < t1-xi < hl1
        # infinite surfaces touch within cyl1
        # get distance from end of cyl2 to cyl1
        for fsign in [1.0, -1.0]  # check both ends of cyl2
            ep = p2 + fsign*hl2*u2  # end point
            dpara, dperp = para_perp_dist(p1, u1, ep)
            if abs(dpara) < hl1
                # end point is next to other cylinder
                if dperp < diam
                    return(1.0)
                end
            else
                # end point ep is beyond end of other cylinder
                # calculate distance between end points
                d = ep - (p1 + sign(dpara)*hl1*u1)
                if sqrt(sum(d.^2)) < diam
                    return(1.0)
                end
            end
        end
        # no collision detected
        return(0.0)
    elseif -hl2 < t2+ff*xi < hl2 || -hl2 < t2-ff*xi < hl2
        # infinite surfaces touch within cyl2
        # get distance from end of cyl1 to cyl2
        for fsign in [1.0, -1.0]  # check both ends of cyl1
            ep = p1 + fsign*hl1*u1  # end point
            dpara, dperp = para_perp_dist(p2, u2, ep)
            if abs(dpara) < hl2
                # end point is next to other cylinder
                if dperp < diam
                    return(1.0)
                end
            else
                # end point ep is beyond end of other cylinder
                # calculate distance between end points
                d = ep - (p2 + sign(dpara)*hl2*u2)
                if sqrt(sum(d.^2)) < diam
                    return(1.0)
                end
            end
        end
        # no collision detected
        return(0.0)
    else
        # hulls of infinite cylinders cross at points outside the cylinders
        # do hulls cross at opposite ends of the cylinders?
        yz1 = (t1+xi)*(t1-xi)
        yz2 = (t2+ff*xi)*(t2-ff*xi)
        if yz1 < 0.0 && yz2 < 0.0
            # hulls cross at opposite ends -> collision possible
            # check ends of cylinders
            e = p1 - hl1*u1
            cp = cross(e-p2, u2)
            dp = sum((e-p2).*u2)
            if sum(cp.^2) < diam^2
                # collision possible
                if -hl2 < dp < hl2
                    # collision!
                    return(1.0)
                end
                d = e - (p2+sign(dp)*hl2*u2) # end-to-end distance
                if sum(d.^2) < diam^2
                    # collision!
                    return(1.0)
                end 
            end

            e = p1 + hl1*u1
            cp = cross(e-p2, u2)
            dp = sum((e-p2).*u2)
            if sum(cp.^2) < diam^2
                # collision possible
                if -hl2 < dp < hl2
                    # collision!
                    return(1.0)
                end
                d = e - (p2+sign(dp)*hl2*u2) # end-to-end distance
                if sum(d.^2) < diam^2
                    # collision!
                    return(1.0)
                end 
            end

            e = p2 - hl2*u2
            cp = cross(e-p1, u1)
            dp = sum((e-p1).*u1)
            if sum(cp.^2) < diam^2
                # collision possible
                if -hl1 < dp < hl1
                    # collision!
                    return(1.0)
                end
                d = e - (p1+sign(dp)*hl1*u1) # end-to-end distance
                if sum(d.^2) < diam^2
                    # collision!
                    return(1.0)
                end 
            end

            e = p2 + hl2*u2
            cp = cross(e-p1, u1)
            dp = sum((e-p1).*u1)
            if sum(cp.^2) < diam^2
                # collision possible
                if -hl1 < dp < hl1
                    # collision!
                    return(1.0)
                end
                d = e - (p1+sign(dp)*hl1*u1) # end-to-end distance
                if sum(d.^2) < diam^2
                    # collision!
                    return(1.0)
                end 
            end
            # no collision detected
            return(0.0)
        end
        # no collision
        return(0.0)
    end
    #  what case is this?
    println("WARNING! what is this?")
    return(-1.0)
end



"""
Get the center of mass of a spherocylinder connecting to a given spherocylinder.
The spherocylinder with orientation u2 connects at end p1+len/2*u1.

spherocylinder_connect(p1::Array{Float64,1}, u1::Array{Float64,1}, l1::Float64, u2::Array{Float64,1}, l2::Float64, dist::Float64)
p1: center-of-mass of the 1st cylinder
u1, u2: orientation unit vectors of the cylinders
l1, l2: lengths of the cylinders
dist: distance of end points: choose radius1+radius2+dr
      dr>0 guarantees that the spherocylinders do not overlap.
"""
function spherocylinder_connect(p1::Array{Float64,1}, u1::Array{Float64,1}, l1::Float64,
    u2::Array{Float64,1}, l2::Float64, dist::Float64)
    # p1: center-of-mass of the 1st cylinder
    # u1, u2: orientation unit vectors of the cylinders
    # l1, l2: lengths of the cylinders
    # dist: distance of end points: choose radius1+radius2+dr
    #       dr>0 guarantees that the spherocylinders do not overlap.
    al = acos(sum(u1.*u2))  # angle between cylinders
    vs = cross(u2,u1)  # perp. vector
    us = vs ./ sqrt(sum(vs.^2))
##    dv = cos((π+al)/2) * u1 + sin((π+al)/2) * cross(u1,us)  # vector towards end-point of cyl2
    dv = cos(al/2.0) * u1 + sin(al/2.0) * cross(u1,us)  # vector towards end-point of cyl2
##    p2 = p1 - l1/2.0*u1 + dist*dv + l2/2.0*u2
    p2 = p1 + l1/2.0*u1 + dist*dv + l2/2.0*u2
    return(p2)
end
