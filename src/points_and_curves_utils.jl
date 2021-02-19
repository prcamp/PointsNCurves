using LinearAlgebra

function boxit(crv::Curve,wbuff::Number,dbuff::Number)
    if x(crv[end]) < x(crv[1])
        wbuff = -wbuff
    end
    p1 = crv[end]+Point(wbuff,0)
    p2 = crv[end]+Point(wbuff,dbuff)
    p3 = crv[1]+Point(-wbuff,dbuff)
    p4 = crv[1]+Point(-wbuff,0)
    push!(crv,p1)
    push!(crv,p2)
    push!(crv,p3)
    push!(crv,p4)
    ClosedCurve(crv)
end

function occlude_sequence(crvs::ClosedCurves)
    # mask from front to back
    n = length(crvs)
    ncrvs = Curves()
    ocrvs = crvs
    push!(ncrvs,Curve(crvs[1]))
    for k in 2:n
        ncrv = Curve(crvs[k])
        o = reverse(ocrvs[1:k-1])
        ncrv = occlude_sequence(ncrv,o)
        ncrvs = vcat(ncrvs,ncrv)
    end
    return ncrvs
end

function occlude_sequence(crvs::Curves)
    # mask from front to back
    n = length(crvs)
    t = linspace(0,1,n)
    dbuff = miny(crvs)-1000
    wbuff0 = 900
    wbuff1 = 300
    ncrvs = Curves()
    ocrvs = ClosedCurves()
    for k in 1:n
        wbuff = (1-t[k])*wbuff0 + t[k]*wbuff1
        push!(ocrvs,boxit(crvs[k],wbuff,dbuff))
    end
    push!(ncrvs,crvs[1])
    for k in 2:n
        ncrv = crvs[k]
        o = flipdim(ocrvs[1:k-1],1)
        ncrv = occlude_sequence(ncrv,o)
        ncrvs = vcat(ncrvs,ncrv)
    end
    return ncrvs
end

function mask_sequence(crvs::Curves,crv::ClosedCurve)
    # mask from front to back
    n = length(crvs)
    ncrvs = Curves()
    push!(ncrvs,crvs[1])
    ncrvs = Curves()
    for k = 1:n
        ncrv,~ = masked_occluded(crvs[k],crv)
        ncrvs = vcat(ncrvs,ncrv)
    end
    return ncrvs
end

function occlude_sequence(crvs::Curves,ocrvs::ClosedCurves)
    # mask from front to back
    ncrvs = Curves()
    push!(ncrvs,crvs[1])
    for c in crvs
        ncrv = occlude_sequence(c,ocrvs)
        ncrvs = vcat(ncrvs,ncrv)
    end
    return ncrvs
end

function occlude_sequence(crvs::Curves,ocrv::ClosedCurve)
    oseq = ClosedCurves()
    push!(oseq,ocrv)
    occlude_sequence(crvs,oseq)
end

function occlude_sequence(crv::Curve,crvs::ClosedCurves)
    # mask from front to back
    ncrv = crv
    for o in crvs
        if length(ncrv) > 0
            ~,ncrv = masked_occluded(ncrv,o)
        end
    end
    return ncrv
end

# function apply_comparison_to_coordinates(c::Curve,z::Char, f::Function)
#    #=
#     f is a function which takes two inputs and returns one of them.
#     e.g. f = min or f = max
#    =#
#    a = c[1].x
#    for i in 2:length(c)
#       x = min(x,c[i].x)
#    end
#    return x
# end


# function inshape_cn(p::Point, c::Curve) #= crossing number test for a point in a polygon
#      The ray is taken to be the x-axis passing through the point p
#           Input:   p = a point,
#                      c[] = vertex points of a polygon V[n+1] with V[n]=V[0]
#           Return:  false = outside, true = inside
#          This code is patterned after [Franklin, 2000]
#     =#
#
#     # Initialize the crossing number:
#     crossing_counter = 0
#     n = length(c)
#     for i in 1:n
#       i0 = i
#       (i == length(c)) ? i1 = 1 : i1 = i+1
#        if ( ((c[i0].y <= p.y) & (c[i1].y >  p.y))  ||   # an upward crossing
#             ((c[i0].y >  p.y) & (c[i1].y <= p.y)) )      #  a downward crossing
#             # compute  the actual edge-ray intersect x-axis
#             edgerayvert = (p.y  - c[i0].y) / (c[i1].y - c[i0].y)
#             if (p.x <  c[i0].x + edgerayvert * (c[i1].x - c[i0].x)) # p.x < intersect
#                  crossing_counter +=1 # a valid crossing of y=P.y right of P.x
#             end
#          end
#     end
#     isin = (mod(crossing_counter,2) == 1)
# end
#
# function intersectionpoints(c1,c2)
#    #=
#    returns index pairs of
#    =#
#
#
#
# end

# Even if no points of c1 lie inside of c2 there may still be points of intersection.
# Also consecutive points can't be guarenteed to

function masked_occluded(c::Curves,c2::ClosedCurve)
    masked = Curves()
    occluded = Curves()
    for crv in c
        m,o = masked_occluded(crv,c2)
        masked = vcat(masked,m)
        occluded = vcat(occluded,o)
    end
    return masked, occluded
end

function masked_occluded(c1::ClosedCurve,c2::ClosedCurve)
    # Partitions c1 into a set of curves.
    c1 = Curve(c1)
    m,o = masked_occluded(c1,c2)
end

function masked_occluded(c1::Curve,c2::ClosedCurve)
    #=
    c1 will occlude c2, i.e. only the parts of c2 lying outside of c1 will remain
    The result is a set of curves of type Curves.

    =#
    #println("masking & occluding...")
    if length(c1)>1
        tt = 0.5
        masked = Curves()
        occluded = Curves()

        m1x = minx(c1)
        M1x = maxx(c1)
        m1y = miny(c1)
        M1y = maxy(c1)
        m2x = minx(c2)
        M2x = maxx(c2)
        m2y = miny(c2)
        M2y = maxy(c2)
        mayintersect = (m1x <= M2x) & (m2x <= M1x) & (m1y <= M2y) & (m2y <= M1y)

        if mayintersect
            #= For each edge of the partitioned curve c2,
            I must find the segments of c1 intersecting them. There may be
            multiples. This operation is O(n*m), but hopefully it will be mostly short circuited
            =#
            isin = inshape(c1[1],c2)
            if isin
                push!(masked,Curve())
            else
                push!(occluded,Curve())
            end

            #println("they might intersect...")
            # numintersects = 0
            for i in 1:length(c1)-1
                s1 = Segment(c1[i],c1[i+1])
                # now find all the segments of c2 which might intersect s1:
                intersectingsegments = Int[]
                intersectionpoints = Point[]
                scale = Float64[]
                for j in 1:c2.length
                    s2 = Segment(c2[j],c2[j+1])
                    doesintersect,s,t,p = segmentintersect(s1,s2)
                    if doesintersect
                        # println("does intersect: \n",
                        # "s is $s \n",
                        # "t is $t \n",
                        # "p is ($(p.x),$(p.y)) \n",
                        # )
                        push!(intersectingsegments,j)
                        # ignoring double intersections for the moment
                        push!(intersectionpoints,p)
                        push!(scale,s)
                    end
                end
                # Now deal with the intersections and push results onto the curves
                # first sort according to scale:
                # println("now to deal with those intersections...")
                idx = sortperm(scale)
                scale = scale[idx]
                intersectingsegments = intersectingsegments[idx]
                intersectionpoints = intersectionpoints[idx]
                # assuming non overlapping points:
                # first segment:
                if isin
                    push!(masked[end],s1.s)
                else
                    push!(occluded[end],s1.s)
                end
                #println("scale is $(length(scale)) elements long")
                for j in 1:length(scale)
                    if j < length(scale)
                        mp = (1-tt)*intersectionpoints[j] + tt*intersectionpoints[j+1]
                    else
                        mp = (1-tt)*intersectionpoints[j] + tt*s1.t
                    end
                    if isin
                        push!(masked[end],intersectionpoints[j])
                        if !strictlyinshape(mp,c2)
                            push!(occluded,Curve())
                            push!(occluded[end],intersectionpoints[j])
                            isin = false
                        end
                    else
                        push!(occluded[end],intersectionpoints[j])
                        if strictlyinshape(mp,c2)
                            push!(masked,Curve())
                            push!(masked[end],intersectionpoints[j])
                            isin = true
                        end
                    end
                end
                if isin
                    push!(masked[end],s1.t)
                else
                    push!(occluded[end],s1.t)
                end
            end
        else
            push!(occluded, c1)
        end
        masked = cleancurve(masked)
        occluded = cleancurve(occluded)
    else
        masked = Curves()
        occluded = Curves()
    end
    return masked, occluded
end

function simplify(c::ClosedCurve)
    # return a curve with no self intersections with outer boundary same as c
    t = .99
    splfd = Curve()
    interiorinds = strictlyinshape(c)
    outinds = setdiff(1:c.length,interiorinds)
    idx = outinds[1]
    p = c[idx]
    push!(splfd,p)
    isclosed = false
    state = idx
    while ! isclosed
        (pa, state) = next(c, state)
        qa = c[state]
        s1 = Segment(pa,qa)
        intersectingsegments = Int[]
        intersectionpoints = Point[]
        scale = Float64[]
        for j in 1:c.length
            if state != j && state+1 != j && state-1 != j
                s2 = Segment(c[j],c[j+1])
                doesintersect,s,t,r = segmentintersect(s1,s2)
                if doesintersect
                    # println("does intersect: \n",
                    # "s is $s \n",
                    # "t is $t \n",
                    # "p is ($(r.x),$(r.y)) \n",
                    # )
                    push!(intersectingsegments,j)
                    # ignoring double intersections for the moment
                    push!(intersectionpoints,r)
                    push!(scale,s)
                end
            end
        end
        # println("now to deal with those intersections...")
        inds = sortperm(scale)
        scale = scale[inds]
        intersectingsegments = intersectingsegments[inds]
        intersectionpoints = intersectionpoints[inds]
        if length(scale)>0
            r = intersectionpoints[1]
            j = intersectingsegments[1]
            push!(splfd,r)

            pb = c[j]
            qb = c[j+1]
            s1a = Segment(pa,r)
            s1b = Segment(r,qa)
            s2a = Segment(pb,r)
            s2b = Segment(r,qb)
            if ! strictlyinshape(t*s1b,c)

            else
                qn = qa
                push!(splfd,qn)
            end
        end
        if qn == p
            isclosed = true
        end
    end
    return ClosedCurve(splfd[1:end-1])
end

function checkcurvesbounds(crvs::Curves)
    initc = crvs[1]
    for i in 2:length(crvs)
        nxtc = crvs[i]
        (initc[end] != nxtc[1]) && begin
        # @show length(crvs)
        # @show i
        # @show initc[end]
        # @show nxtc[1]
        return false
    end
    initc = nxtc
end
true
end

function unify(c1::Curves,c2::Curves)::Curves

    # sanity check:

    !checkcurvesbounds(c1) && begin
    error("c1: bad bounds")
end
!checkcurvesbounds(c2) && begin
error("c2: bad bounds")
end
# if !isempty(c1) && !isempty(c2)
# @show c1[1][1]
# @show c2[1][1]
# @show c1[end][end]
# @show c2[end][end]
# end
# @show length(c1)
# @show length(c2)
# println("-----------ok-----------")
if (c1 == c2)
    # println("kkkkkkkkkkkkkkkkkkkk")
    return c1
else

    if length(c1)==1 && length(c2) > 0
        return c2
    elseif length(c2)==1 && length(c1) > 0
        return c1
    end

    # remove any empty curves:
    # filter!(crv -> !isempty(crv), c1)
    # filter!(crv -> !isempty(crv), c2)

    # Assuming that c1 and c2 are partitions of the same curve.
    if length(c1) == 0
        # println("lllllllllllllllaaaaaaaaa")
        # @show length(c1)
        # @show length(c2)
        return c2
    end
    if length(c2) == 0
        # println("uuuuuuuuuuuhhhhhhhhhhhhhh")
        return c1
    end

    # Now there are curves in there

    (c1[1][1] != c2[1][1]) && begin
    # @show c1[1][1]
    # @show c2[1][1]
    # @show length(c1)
    # @show length(c2)
    map(c -> length(c),c1)
    map(c -> length(c),c2)
    error("unify: curves must agree on where to start")
end
i1 = c1[1]
i2 = c2[1]

if length(i1) == length(i2)
    if i1[end] == i2[end]
        # then they're the same. unify the rest of it
        # println("hhhhhhhhhiiiiiiiiii")
        return pushfirst!(unify(c1[2:end],c2[2:end]), i1)
    else
        # Now they should only differ in the last entry:
        if (length(c1) < 2 || length(c2) < 2)
            # println("mmmmmmmmmmmmmmmmmmmmm")
            return c1
        end
        # error("unify: curves must agree on they're last entry")


        # assuming curves join nicely at endpoints within each curve set.
        # so i1[end] = n1[1]

        # basically need to check which one comes first:

        sp11 = i1[end-1]
        sp12 = i2[end-1]

        sp21 = c1[2][2]
        sp22 = c2[2][2]

        p0,p1,p2,p3 = sort_colinear_points(sp11,sp12,sp21,sp22)

        # sp11 - sp21 and sp12 - sp22 should be colinear


        ip1 = i1[end]
        ip2 = i2[end]
        if ip1 == ip2
            # println("eeeeeeeeeehhhhhhhhhh")
            return pushfirst!(unify(c1[2:end],c2[2:end]),i1)
        end


        one_is_less = insegment(p0,ip1,ip2)
        if one_is_less == 0
            # then ip1 < ip2
            # i1 is initial.
            intercurve = Curve([ip1,ip2])
            # and now must modify n1 to start at ip2
            n1 = c1[2]
            n12 = pushfirst!(n1[2:end],ip2)
            c11 = if length(c1) >= 3
                pushfirst!(c1[3:end],n12)
            else
                Vector([n12])
            end
            # println("oooooohohhhhhohohoh")
            return pushfirst!(pushfirst!(unify(c11,c2[2:end]),intercurve),i1)
        elseif one_is_less == 1
            # then ip2 < ip1
            # i2 is initial.
            intercurve = Curve([ip2,ip1])
            # and now must modify n2 to start at ip1
            n2 = c2[2]
            n22 = pushfirst!(n2[2:end],ip1)
            c22 = if length(c2) >= 3
                pushfirst!(c2[3:end],n22)
            else
                Vector([n22])
            end
            # println("aaaaaaahhhhhhhhhhh")
            return pushfirst!(pushfirst!(unify(c1[2:end],c22),intercurve),i2)
        else
            one_is_less = insegment(ip1,ip2,p3)
            if one_is_less == 0
                # then ip1 < ip2
                # i1 is initial.
                intercurve = Curve([ip1,ip2])
                # and now must modify n1 to start at ip2
                n1 = c1[2]
                n12 = pushfirst!(n1[2:end],ip2)
                c11 = if length(c1) >= 3
                    pushfirst!(c1[3:end],n12)
                else
                    Vector([n12])
                end
                # println("oooooohohhhhhohohoh")
                return pushfirst!(pushfirst!(unify(c11,c2[2:end]),intercurve),i1)
            elseif one_is_less == 1
                # then ip2 < ip1
                # i2 is initial.
                intercurve = Curve([ip2,ip1])
                # and now must modify n2 to start at ip1
                n2 = c2[2]
                n22 = pushfirst!(n2[2:end],ip1)
                c22 = if length(c2) >= 3
                    pushfirst!(c2[3:end],n22)
                else
                    Vector([n22])
                end
                # println("aaaaaaahhhhhhhhhhh")
                return pushfirst!(pushfirst!(unify(c1[2:end],c22),intercurve),i2)

            else
                one_is_less = insegment(p3,ip1,ip2)
                if one_is_less == 0
                    # then ip1 < ip2
                    # i1 is initial.
                    intercurve = Curve([ip1,ip2])
                    # and now must modify n1 to start at ip2
                    n1 = c1[2]
                    n12 = pushfirst!(n1[2:end],ip2)
                    c11 = if length(c1) >= 3
                        pushfirst!(c1[3:end],n12)
                    else
                        Vector([n12])
                    end
                    # println("oooooohohhhhhohohoh")
                    return pushfirst!(pushfirst!(unify(c11,c2[2:end]),intercurve),i1)
                elseif one_is_less == 1
                    # then ip2 < ip1
                    # i2 is initial.
                    intercurve = Curve([ip2,ip1])
                    # and now must modify n2 to start at ip1
                    n2 = c2[2]
                    n22 = pushfirst!(n2[2:end],ip1)
                    c22 = if length(c2) >= 3
                        pushfirst!(c2[3:end],n22)
                    else
                        Vector([n22])
                    end
                    # println("aaaaaaahhhhhhhhhhh")
                    return pushfirst!(pushfirst!(unify(c1[2:end],c22),intercurve),i2)

                else
                    error("should be unreachable")
                end
            end
        end
    end
elseif length(i1) < length(i2)
    # then i1 is initial. chop i2:
    n2 = pushfirst!(i2[length(i1):end],i1[end])
    # @show i1[end] == c1[2][1]
    # @show i1[end]
    # @show c1[2][1]
    # @show i1  == c1[1]
    # @show i2[end] == c2[2][1]
    # @show i2[end]
    # @show c2[2][1]
    # @show i2  == c2[1]

    c22 = if length(c2) >= 2
        pushfirst!(c2[2:end], n2)
    else
        Vector([n2])
    end
    # println("bbbbbbbbbbaaaaaaaaaaaa")
    return pushfirst!(unify(c1[2:end],c22),i1)
else
    # then i2 is initial. chop i1:
    # @show i2[length(i1)] == c1[2][2]
    # @show i1[end] == c1[2][1]
    # @show i1[end]
    # @show c1[2][1]
    # @show i1  == c1[1]
    # @show i2[end] == c2[2][1]
    # @show i2[end]
    # @show c2[2][1]

    n1 = pushfirst!(i1[length(i2):end],i2[end])
    c11 = if length(c1) >= 2
        pushfirst!(c1[2:end], n1)
    else
        # println("short c1")
        # @show length(c1)
        Vector([n1])
    end
    # println("aaaaaaaaakkkkkkkkkkkkk")

    return pushfirst!(unify(c11,c2[2:end]),i2)

end
end

end

function collect_self_intersection_points(crvs::Curves)
    intersection_indices = Vector{Tuple{Int,Int}}()
    crv_to_inds = Dict{Int,Set{Int}}()
    function add_to_dict(i,j)
        if !haskey(crv_to_inds,i)
            crv_to_inds[i] = Set{Int}()
        end
        push!(crv_to_inds[i],j)
    end
    intind = 1
    for (i,c) in enumerate(crvs)
        for j in i+1:length(crvs)
            if c[end] == crvs[j][end]
                push!(intersection_indices, (i,j))
                add_to_dict(i , intind)
                add_to_dict(i+1 , intind)
                add_to_dict(j , intind)
                add_to_dict(j+1 , intind)
                intind += 1
            end
        end
    end
    intersection_indices, crv_to_inds
end

function self_intersect_curve!(crv::Curve, point2ind::Dict{Point,Set{Int}}=Dict{Point,Set{Int}}())
    # returns the curve modified along with a dictionary from points of intersections
    # to where they occur as indices in the curve

    hadintersect, c1,c2,c3,i,j = self_intersect_curve_first_with_indices(crv)
    if hadintersect
        intpnt = c2[1]
        # assumption is that point2ind doesn't have intpnt yet
        push!(point2ind,Set{Int}())
        #### TODO!!!!!!!!!!!!!
    end
end

function find_self_intersection_points(crv::Curve)
    p0 = if crv[1]==crv[end]
        crv[1]
    else
        Point(-10000*rand()-10000,-10000*rand()-10000)
    end
    points = Set{Point}()
    indices = Vector{Tuple{Int,Int}}()
    inds = Dict{Int,Vector{Tuple{Number,Point}}}()
    scalers = Vector{Tuple{Number,Number}}()
    if length(crv) < 2
        points, indices, scalers, inds
    else
        for i = 1:length(crv)-3
            seg0 = Segment(crv[i],crv[i+1])
            for j = i+2:length(crv)-1
                seg1 = Segment(crv[j],crv[j+1])
                doesintersect, s, t, p = if i != 1 && j+1 != length(crv)
                    segmentintersect(seg0,seg1, include_ends=true)
                else
                    segmentintersect(seg0,seg1, include_ends=false)
                end
                if doesintersect && p != p0
                    if !haskey(inds,i)
                        inds[i] = Vector{Tuple{Number,Point}}()
                    end
                    if !haskey(inds,j)
                        inds[j] = Vector{Tuple{Number,Point}}()
                    end
                    push!(points,p)
                    push!(indices,(i,j))
                    push!(scalers,(s,t))
                    push!(inds[i],(s,p))
                    push!(inds[j],(t,p))
                end
            end
        end
        points, indices, scalers, inds
    end
end

function compute_self_intersections(crv::Curve)
    points, indices, scalers, inds2points = find_self_intersection_points(crv)
    # Now I want to insert the points into a new curve.
    newCurve = Curve()
    pointIndices = Dict{Point,Set{Int}}()
    for p in points
        pointIndices[p] = Set{Int}()
    end
    ind = 1
    for i = 1:length(crv)-1
        p = crv[i]
        push!(newCurve,p)
        if in(p, points)
            push!(pointIndices[p],ind)
        end
        ind += 1
        if haskey(inds2points,i)
            srtinds = sortperm(map(spnt -> spnt[1],inds2points[i]))
            toconcat = map(indpp -> indpp[2], filter(indp -> indp[1] != 0 && indp[1] != 1, inds2points[i][srtinds]))
            for j = 1:length(toconcat)
                pnt = toconcat[j]
                if !in(pnt,points)
                    error("intersection point not in known set of intersection points")
                end
                push!(newCurve,pnt)
                push!(pointIndices[pnt],ind)
                ind += 1
            end
        end
    end
    push!(newCurve,crv[end])
    newCurve, pointIndices
end

showitNbreak = false

function closed_curve_to_faces(c::ClosedCurve, testsclr = .000001)
    # test scaler

    testpoints = Curve()
    faces = Vector{ClosedCurve}()
    newcrv, pointIndices = compute_self_intersections(Curve(c))

    intpoints = keys(pointIndices)
    len = length(newcrv)-1
    face_verts = Vector{Set{Tuple{Point,Int}}}()
    seenAtVerts = Dict{Point,Set{Int}}()
    if isempty(pointIndices)
        push!(faces,newcrv)
        push!(face_verts,Set{Tuple{Point,Int}}())
    else
        vertOrds = compute_face_vertices(newcrv, pointIndices)
        for p in keys(pointIndices)
            seenAtVerts[p] = Set{Int}()
        end
        # Need to go into each intersection point, and find faces.
        for (p,v) in vertOrds
            for i in 1:length(v)
                showitNbreak && @show p

                showitNbreak && @show i
                showitNbreak && printdict(seenAtVerts)
                showitNbreak && @show v


                initedge = v[i]

                showitNbreak && @show initedge
                nbredge = v[nextindx(i,length(v),1)]#initedge[3])]
                showitNbreak && @show nbredge
                nbrpnt = newcrv[nbredge[2]] + (-1) * p
                nbrpnt = 1/(norm(nbrpnt)) * nbrpnt

                thisidx = initedge[2]
                thisdir = initedge[3]
                thispnt = newcrv[thisidx] + (-1) * p
                thispnt = 1/norm(thispnt) * thispnt


                # test if the face is in the polygon:
                testp = 1/2 * (thispnt + nbrpnt)
                testp = 1/norm(testp) * testp
                testp = p + testsclr * testp
                push!(testpoints,testp)
                isin = strictlyinshape(testp,ClosedCurve(newcrv))
                showitNbreak && @show isin
                # print("break or continue: ")
                # msg = readline()
                # if msg=="b"
                #     return faces, face_verts, testpoints
                # end
                if !isin
                    showitNbreak && println("that one isn't in")
                    push!(seenAtVerts[p],i)
                else
                    if !in(i,seenAtVerts[p])
                        # make a face!
                        showitNbreak && println("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

                        push!(seenAtVerts[p],i)


                        face = Curve()
                        push!(face,p)
                        # point to baseidx and next idx
                        thisfaceverts = Set{Tuple{Point,Int}}()
                        push!(thisfaceverts,(p,i))
                        # this is the index of the intersection point:
                        initidx = nextindx(thisidx,len,-thisdir)

                        vidx = i #findfirst(edge -> edge[2]==thisidx,vertOrds[p])


                        curridx = initidx
                        currdir = thisdir
                        currpoint = p

                        showitNbreak && println("makin a face")
                        showitNbreak && @show currpoint
                        showitNbreak && @show vidx


                        isclosed = false
                        # println("here we go again")
                        while !isclosed
                            # println("------------------------------------------------")
                            # march along the curve until you reach the next intersection point:
                            isintpoint = false
                            # @show newcrv

                            lastidx = curridx

                            while !isintpoint
                                # @show
                                lastidx = curridx
                                curridx = nextindx(curridx,len,currdir)
                                currpoint = newcrv[curridx]
                                isintpoint = in(currpoint,intpoints)
                                push!(face,currpoint)
                            end
                            showitNbreak && @show lastidx
                            showitNbreak && @show curridx
                            showitNbreak && @show currpoint
                            # @show vertOrds[currpoint]
                            # get the edge info of the direction we're coming in on

                            intpointedges = vertOrds[currpoint]
                            vidx = findfirst(edge -> edge[2]==lastidx,intpointedges)

                            push!(thisfaceverts,(currpoint,vidx))
                            # now we know we're at an intersection point, and the last index indicates
                            # where we were coming from
                            # (lastidx != previntpointidx) && error("I thought we'd have lastidx == previntpointidx")
                            # @show prevpointidx
                            # @show face

                            # @show intpointedges
                            prevpointidx_in_vertOrd = vidx
                            # for j in 1:length(intpointedges)
                            #     edgeinfo = intpointedges[j]
                            #     if edgeinfo[2] == lastidx
                            #         prevpointidx_in_vertOrd = j
                            #         break
                            #     end
                            # end
                            # if prevpointidx_in_vertOrd != vidx
                            #     error("those should really be the same I would think")
                            # end
                            # @show intpointedges
                            # @show map(idx -> newcrv[idx[2]],intpointedges)
                            # @show prevpointidx_in_vertOrd

                            # prevpointidx_in_vertOrd == 0 && error("that index should be nonzero")

                            # push!(thisfaceverts,(currpoint,prevpointidx_in_vertOrd))
                            # Now find the next edge:

                            curredgeinfo = intpointedges[vidx]



                            # @show curredgeinfo
                            # @show prevpointidx

                            nidx = nextindx(vidx,length(intpointedges),-1)#initedge[3])
                            # @show nidx
                            push!(seenAtVerts[currpoint], nidx)
                            nextedgeout = intpointedges[nidx]
                            nextpoint = newcrv[nextedgeout[2]]
                            # @show nextedgeout
                            # @show initedge
                            # @show nbredge
                            # @show curredgeinfo
                            # @show newcrv[nextedgeout[2]]
                            # @show newcrv[initedge[2]]
                            # @show newcrv[curredgeinfo[2]]
                            # @show nextedgeout == initedge

                            # check to see if we're done.
                            # @show faces
                            if curredgeinfo == nbredge  #|| initedge == nextedgeout
                                showitNbreak && println("face made")
                                showitNbreak && println("yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy")
                                push!(faces,ClosedCurve(face))
                                push!(face_verts,thisfaceverts)
                                isclosed = true
                                face = Curve()
                                thisfaceverts = Set{Tuple{Point,Int}}()
                            else
                                # push!(face,nextpoint)
                                currdir = nextedgeout[3]
                                lastidx = curridx
                                curridx = nextindx(nextedgeout[2],len,-currdir)
                                vidx = nidx
                                # (nidx != vidx) && error("they look like they should be equal")
                                push!(thisfaceverts,(currpoint,vidx))
                            end
                        end
                    end
                end
            end
        end
    end

    ufaces = Vector{ClosedCurve}()
    fverts = typeof(face_verts)()
    for (f,v) in zip(faces,face_verts)
        if !in(v,fverts)
            push!(ufaces,f)
            push!(fverts,v)
        end
    end
    ufaces, fverts, testpoints
end

function compute_face_vertices(c::Curve,pointIndices::Dict{Point,Set{Int}})
    # determine the outgoing edges and their angles
    vertOrds = Dict{Point,Vector{Tuple{Number,Int,Int}}}()
    # @show newcrv
    # @show pointIndices
    len = length(c)
    for (p,inds) in pointIndices
        # do some Initializations:
        # @show p
        thisvec = Vector{Tuple{Number,Int,Int}}()
        # @show inds
        for id in inds
            backidx = nextindx(id,len,-1)
            backpnt = c[backidx]
            backedge = backpnt + (-1) * p
            btheta = pointθ(backedge)
            push!(thisvec,(btheta,backidx,-1))
            nextidx = nextindx(id,len,1)
            nextpnt = c[nextidx]
            nextedge = nextpnt + (-1) * p
            ntheta = pointθ(nextedge)
            push!(thisvec,(ntheta,nextidx,1))
        end
        angles = map(idthdir -> idthdir,thisvec)
        if length(angles) != length(unique(angles))
            error("assuming that each outgoing edge from an intersection point has a unique angle")
        end
        srtinds = sortperm(angles)
        thisvec = thisvec[srtinds]
        vertOrds[p] = thisvec
    end
    vertOrds
end

function compute_neighbors(face_verts::Vector{Set{Tuple{Point,Int}}})::Vector{Set{Int}}
    outVec = Vector{Set{Int}}()
    for i in 1:length(face_verts)
        push!(outVec,Set{Int}())
    end
    for i in 1:length(face_verts)-1
        for j in i+1:length(face_verts)
            if !isempty(intersect(face_verts[i],face_verts[j]))
                push!(outVec[i],j)
                push!(outVec[j],i)
            end
        end
    end
    outVec
end

set2vec(s::Set) = filter(i -> in(i,s),s.dict.keys)

function curve_union(c1::ClosedCurve,c2::ClosedCurve,poly_fill_type=PolyFillTypeNonZero)::Vector{ClosedCurve}
    p1 = topath(c1)
    p2 = topath(c2)
    clp = Clip()
    add_path!(clp,p1,PolyTypeSubject,true)
    add_path!(clp,p2,PolyTypeClip,true)
    res,os = execute(clp,ClipTypeUnion, poly_fill_type, poly_fill_type)
    !res && error("failed to take the union of the curves")
    map(p -> ClosedCurve(frompath(p)),os)
end

function curve_union(cs::Vector{ClosedCurve})::Vector{ClosedCurve}
    if length(cs)==1 || length(cs) == 0
        return cs
    else
        c1 = cs[1]
        c2 = cs[2]
        uc = curve_union(c1,c2)
        for c3 in cs[3:end]
            nc = Vector{ClosedCurve}()
            for c4 in uc
                # @show typeof(c4)
                nc = vcat(nc,curve_union(c3,c4))
            end
            uc = nc
        end
        uc
    end
end

function nextindx(i,len,dir)
    if i > 1 && i < len
        i + dir
    elseif i == 1
        if dir == 1
            2
        else
            len
        end
    else
        if dir == -1
            len - 1
        else
            1
        end
    end
end

function fuse_faces(faces::Vector{ClosedCurve},neighbors::Vector{Set{Int}},i::Int,j::Int)
    if i == j
        return faces, neighbors, i
    else
        newface = curve_union(faces[i],faces[j])
        if (length(newface) > 1)
            # println("faces are assumed to be neighbors. got $i and $j")
            return faces, neighbors, 0
        end

        newface = newface[1]
        i0 = min(i,j)
        i1 = max(i,j)
        faces[i0] = newface
        neighbors[i0] = union(neighbors[i0],neighbors[i1])
        deleteat!(faces,i1)
        deleteat!(neighbors,i1)
        # now need to remap the indices above i1 within neighbors:
        remapfcn = k -> (k > i1) ? k - 1 : (k == i1) ? i0 : k
        for nbr in neighbors
            replace!(remapfcn, nbr)
        end
        neighbors[i0] = setdiff(neighbors[i0], Set(i0))
        return faces, neighbors, i0
    end
end

function fuse_faces(faces::Vector{ClosedCurve},neighbors::Vector{Set{Int}},indices::Vector{Int})
    if isempty(indices)
        return faces, neighbors, 0
    else
        filter!(i -> 1 <= i <= length(faces),indices)
        faces_to_merge = map(i -> faces[i],indices)
        newface = curve_union(faces_to_merge)
        if (length(newface) > 1)
            # println("faces are assumed to be neighbors.")
            # @show indices
            # plt(newface)
            return faces, neighbors, 0
        end
        newface = newface[1]

        i0 = minimum(indices)

        new_neighbors = Set{Int}()

        for i in indices
            new_neighbors = union(new_neighbors, neighbors[i])
        end
        faces[i0] = newface
        neighbors[i0] = new_neighbors
        inds_to_drop = sort(setdiff(indices,[i0]))
        if length(unique(inds_to_drop)) < length(inds_to_drop)
            error("that shouldn't happen")
        end
        deleteat!(faces,inds_to_drop)
        deleteat!(neighbors,inds_to_drop)

        # now need to remap the indices above i1 within neighbors:
        function remapfcn(k)
            if k < i0
                k
            elseif in(k,inds_to_drop)
                i0
            else
                takeoff = length(filter(i -> i <= k, inds_to_drop))
                k - takeoff
            end
        end
        # !isempty(indices) && map!(remapfcn,indices)
        for nbr in neighbors
            replace!(remapfcn, nbr)
        end
        neighbors[i0] = setdiff(neighbors[i0], Set(i0))
        faces, neighbors, i0
    end
end

function drop_face(faces::Vector{ClosedCurve}, neighbors::Vector{Set{Int}}, i0::Int)
    if i0 > 0 && i0 <= length(faces)
        deleteat!(faces,i0)
        deleteat!(neighbors,i0)
        for n in neighbors
            filter!(j -> j != i0, n)
        end
        inds_to_drop = [i0]
        function remapfcn(k)
            if k < i0
                k
            elseif k == i0
                i0
            else
                takeoff = 1
                k - takeoff
            end
        end
        # !isempty(indices) && map!(remapfcn,indices)
        for nbr in neighbors
            replace!(remapfcn, nbr)
        end
    else
        error("nope this should never happen that the input i doesn't correspond to a face")
    end
    if length(faces) != length(neighbors)
        # @show length(faces)
        # @show length(neighbors)
        error("we should never have this:")

    end
    faces, neighbors
end

function fuse_faces_depthwise_with_branching(faces::Vector{ClosedCurve},neighbors::Vector{Set{Int}},depth::Int, branches::Int)
    faces,neighbor,i0 = fuse_faces_depthwise(faces,neighbors,depth)
    for i in 2:branches
        j0 = i0
        faces,neighbor,i0 = fuse_faces_depthwise(faces,neighbors,depth,i0)
        if i0 == 0
            i0 = j0
        end
    end
    faces, neighbor, i0
end

function fuse_faces_depthwise(faces::Vector{ClosedCurve},neighbors::Vector{Set{Int}},depth::Int, initial_face::Int = 0)
    # choose an initial index:
    j = (initial_face==0) ? ceil(Int,length(faces)*rand()) : initial_face
    indices = Vector{Int}([j])
    # now from there fuse neighbor faces without backtracking until a depth is achieved
    seen_nbrs = Set{Int}()
    if j > length(faces)
        return faces, neighbors, 0
    end
    # println(j)
    for i in 1:depth
        if j <= length(neighbors)
            optionaldirs = setdiff(neighbors[j],seen_nbrs)
            if isempty(optionaldirs)
                # println("ran out of directions")
                break
            end
            k = rand(optionaldirs)
            push!(indices, k)
            seen_nbrs = union(seen_nbrs,neighbors[j])
            push!(seen_nbrs,j)
            j = k
        else
            # @show j
            # @show length(neighbors)
            break
        end
        # println(j)
    end
    i0 = minimum(indices)
    faces, neighbors, i0 = fuse_faces(faces,neighbors,indices)
    faces, neighbors, i0
end

function randomly_fuse_faces(faces::Vector{ClosedCurve},neighbors::Vector{Set{Int}},itrs::Int)
    for  i = 1:itrs
        j = ceil(Int,length(faces)*rand())
        k = rand(neighbors[j])
        faces, neighbors = fuse_faces(faces, neighbors, j, k)
    end
    faces, neighbors
end


function self_intersect_curve(crv::Curve)
    hadintersect, c1,c2,c3 = self_intersect_curve_first(crv)
    if !hadintersect
        Vector([c2])
    else
        # println("in it now")

        # now we know that there isn't any intersection c1 with c1.
        # possible other intersections are
        # between:
        #  c2 and c2
        #  c3 and c2
        #  c3 and c3

        c22 = self_intersect_curve(c2)

        c23, c32 = intersect_curves(c2,c3)

        c33 = self_intersect_curve(c3)
        # @show c22[1][1]
        # @show c23[1][1]
        # @show c22[end][end]
        # @show c23[end][end]
        midsec = unify(c22,c23)

        endsec = unify(c32,c33)

        # println("post first intersect_curves")

        # println("post inner self_intersect_curve")
        # println("finished computing end section")
        # endendsec = Curves()
        # for cA in enumerate(endsec)
        #     endendsec = vcat(endendsec, self_intersect_curve(cA))
        # end
        return pushfirst!(vcat(midsec,endsec),c1)
    end
end


function self_intersect_curve_first_with_indices(crv::Curve)
    # assuming the curve doesn't have multiple segments intersecting at the same point
    l = length(crv)
    (l <= 2) && begin
    return false, Curve(), crv, Curve(), 0, 0
end
pinit1 = crv[1]
for i in 2:l-1
    seg0 = Segment(pinit1,crv[i])
    pinit1 = crv[i]
    pinit2 = crv[i]
    for j in (i+1):l
        seg1 = Segment(pinit2, crv[j])
        pinit2 = crv[j]
        doesintersect, s, t, p = segmentintersect(seg0,seg1, include_ends=false)
        if doesintersect
            crv1 = push!(crv[1:i-1],p)
            crv2 = push!(vcat([p],crv[i:j-1]),p)
            crv3 = vcat([p],crv[j:end])
            return true, crv1, crv2, crv3, i, j+1
        end
    end
end
return false, Curve(), crv, Curve(), 0, 0
end

function self_intersect_curve_first(crv::Curve)
    hadintersect, crv1, crv2, crv3, ~, ~ = self_intersect_curve_first_with_indices(crv)
end


# function intersect_curves(crv1::Vector{Curve},crv2::Vector{Curve})
#     for (i,c1) in enumerate(crv1)
#         for (j,c2) in enumerate(crv2)
#             dointersect, c11, c22 = intersect_curves(c1,c2)
#             crv1 = vcat(crv)


function collect_intersection_indices(crvs1::Curves,crvs2::Curves)
    intersection_indices = Vector{Tuple{Int,Int}}()
    if length(crvs1) < 2 || length(crvs2) < 2
        intersection_indices
    else
        for (i,c1) in enumerate(crvs1)
            for (j,c2) in enumerate(crvs2)
                if (c1[end]==c2[end])
                    push!(intersection_indices,(i,j))
                end
            end
        end
        intersection_indices
    end
end

function intersect_curves(crv1::Curve,crv2::Curve)
    # ignores self intersections
    hadintersect, c1, c2 = intersect_curves_first(crv1,crv2)
    if !hadintersect
        return c1, c2
    else
        # can assume that the lengths of c1 and c2 are both > 1
        # infact assume that length(c1)==length(c2)==2
        (length(c1) != 2 || length(c2) != 2) && begin
        error("curve sets should both have length 2")
    end

    c11 = c1[1]
    c12 = c1[2]

    c21 = c2[1]
    c22 = c2[2]

    # we know there's no intersection between c11 and c21
    initc1, endc2 = intersect_curves(c11,c22)
    endc1, initc2 = intersect_curves(c12,c21)
    midc1,midc2 = intersect_curves(c12,c22)
    # and now we need to check for intersections in the ends:
    ends1 = unify(endc1,midc1)
    ends2 = unify(endc2,midc2)
    # ends1 = Curves()
    # ends2 = Curves()
    # for c111 in endc1
    #     for c222 in endc2
    #         next1, next2 = intersect_curves(c111,c222)
    #         ends1 = vcat(ends1,next1)
    #         ends2 = vcat(ends2,next2)
    #     end
    # end
    return vcat(initc1,ends1), vcat(initc2,ends2)
end
end

function intersect_curves_first(crv1::Curve, crv2::Curve)
    l1 = length(crv1)
    l2 = length(crv2)
    if (l1 <= 1 || l2 <= 1)
        false, Curves([crv1]), Curves([crv2])
    else
        mx1 = minx(crv1)
        my1 = miny(crv1)
        Mx1 = maxx(crv1)
        My1 = maxy(crv1)
        mx2 = minx(crv2)
        my2 = miny(crv2)
        Mx2 = maxx(crv2)
        My2 = maxy(crv2)
        # test for potential intersection:
        # if (Mx2 <= mx1 || My2 <= my1 || Mx1 <= mx2 || My1 <= my2)
        #     false, Curves([crv1]), Curves([crv2])
        # else
        if true
            pinit1 = crv1[1]
            for i in 2:l1
                seg1 = Segment(pinit1,crv1[i])
                pinit1 = crv1[i]
                pinit2 = crv2[1]
                for j in 2:l2
                    seg2 = Segment(pinit2, crv2[j])
                    pinit2 = crv2[j]
                    doesintersect, ~, ~, p = segmentintersect(seg1,seg2, include_ends=false)
                    if doesintersect
                        c11 = push!(crv1[1:i-1],p)
                        c12 = pushfirst!(crv1[i:end],p)

                        c21 = push!(crv2[1:j-1],p)
                        c22 = pushfirst!(crv2[j:end],p)
                        return true, Curves([c11,c12]),Curves([c21,c22])
                    end
                end
            end
            false, Curves([crv1]), Curves([crv2])
        end
    end
end

#
# function join(c1::ClosedCurve,c2::ClosedCurve)
#   inds = setdiff(1:c1.length,inshape(c1,c2))
#   if isempty(inds)
#     joined = join(c2,c1)
#   else
#     joined = Curve()
#     startidx = inds[1]
#     startpnt = c1[inds[1]]
#     p = startpnt
#     push!(joined,p)
#     closed = false
#     while ! closed
#
#
#     end
#     return joined
#   end
#
#   push!(c,c1)
#   push!(c,c2)
#   joined = Curve()
#   # first find a point of c1 not lying in c2: or none exists, a point in c2 not in c1.
#   haveapoint = false
#   startingindex = 0
#   crv = 1
#   for i in 1:c1.length
#     if ! inshape(c1[i],c2)
#       haveapoint = true
#       startingindex = i
#       break
#     end
#   end
#   if haveapoint == false
#     crv = 2
#     for i in 1:c1.length
#       if ! inshape(c1[i],c2)
#         haveapoint = true
#         startingindex = i
#         break
#       end
#     end
#   end
#   if haveapoint == false
#     error("No vertices of either curve lie outside of the other curve! \n Congradulation, because I'm pretty sure that's impossible!")
#   else
#     startingpoint = c[crv][startingindex]
#     p = startingpoint
#     push!(joined,p)
#
#
#
#
#
#
#
#   end
# end


# # find the points of c1 which may lie in c2, and vice versa
# subinds1 = find(p -> ((m2x <= x(p) <= M2x) & (m2y <= y(p) <= M2y)),c1)
# subinds2 = find(p -> ((m1x <= x(p) <= M1x) & (m1y <= y(p) <= M1y)),c2)
#
# # Now I need to find which segments actually intersect
# # Just going to be lazy for now:
# intinds1 = Array{Int,2}()
# intinds2 = Array{Int,2}()
# intpnts = Point[]
# for i in subinds1
#    for j in subinds2
#       flag,pnt,s,t = segmentintersect(c1[[i-1,i]],c2[[j-1,j]])
#       if flag
#          intinds1 = vcat(intinds1,[i-1 i])
#          intinds2 = vcat(intinds2,[j-1 j])
#          push!(intpnts,pnt)
#       end
#       flag,pnt,s,t = segmentintersect(c1[[i,i+1]],c2[[j-1,j]])
#       if flag
#          intinds1 = vcat(intinds1,[i-1 i])
#          intinds2 = vcat(intinds2,[j-1 j])
#          push!(intpnts,pnt)
#       end
#       flag,pnt,s,t = segmentintersect(c1[[i-1,i]],c2[[j,j+1]])
#       if flag
#          intinds1 = vcat(intinds1,[i-1 i])
#          intinds2 = vcat(intinds2,[j-1 j])
#          push!(intpnts,pnt)
#       end
#       flag,pnt,s,t = segmentintersect(c1[[i,i+1]],c2[[j,j+1]])
#       if flag
#          intinds1 = vcat(intinds1,[i-1 i])
#          intinds2 = vcat(intinds2,[j-1 j])
#          push!(intpnts,pnt)
#       end
#    end
# end

# sequences of points between points of actual intersection either remiain,
# or are excised. Given an intersection find out what part is in and which is out.
# Need to deal with the case of multiple intersections, i.e. one line which passes through
# several on the other curve. Only intersections between the curves matter

#
# else
#    push!(c3,c2)
# end
# return c3
#end
#
# function mask(c1::Curve,c2::Curve)
#    #=
#    c1 will be mask of c2, i.e. only the parts of c2 lying in c1 will remain
#    The result will be a set of curves of type Curves
#    =#
#
#    m1x = minx(c1)
#    M1x = maxx(c1)
#    m1y = miny(c1)
#    M1y = maxy(c1)
#    m2x = minx(c2)
#    M2x = maxx(c2)
#    m2y = miny(c2)
#    M2y = maxy(c2)
#
# end

function isleft(p0::Point,p1::Point,p2::Point)
    #=    Input:  three points P0, P1, and P2
    Return: 1 for P2 left of the line through P0 and P1
    0 for P2  on the line
    -1 for P2  right of the line
    =#
    sign( (p1.x - p0.x) * (p2.y - p0.y) - (p2.x -  p0.x) * (p1.y - p0.y) )
end

function strictlyinshape( p::Point, c::ClosedCurve) #= winding number test for a point in a polygon
    Input:   P = a point,
    V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
    Return:  wn = the winding number (=0 only when P is outside)
    =#

    # Initialize the winding number:
    winding_number = 0
    n = c.length
    for i in 1:n
        i0 = i
        i1 = i+1
        if (c[i0].y < p.y)
            if (c[i1].y  > p.y)      # an upward crossing
                if (isleft( c[i0], c[i1], p) > 0)  # p left of  edge
                    winding_number += 1            # have  a valid up intersect
                end
            end
        else                         # start y > P.y (no test needed)
            if (c[i1].y  < p.y) #     // a downward crossing
                if (c[i0].y > p.y)
                    if (isleft( c[i0], c[i1], p) < 0)  #// P right of  edge
                        winding_number -= 1    # have  a valid down intersect
                    end
                end
            end
        end
    end
    isin = (winding_number != 0)
end

function isinteriorsection(c::Curve, cc::ClosedCurve)
    # assuming that c is a subsegment of cc bounded by points of self intersection
    if length(c) <= 1
        false
    else
        isit = strictlyinshape(.5*(c[1] + c[2]),cc)
        if length(c) > 2
            isit = isit || strictlyinshape(.5*(c[2] + c[3]),cc)
        end
        isit
    end
end


function strictlyinshape(c::ClosedCurve)
    inds = find(p -> strictlyinshape(p,c),c)
end

function inshape( p::Point, c::ClosedCurve) #= winding number test for a point in a polygon
    Input:   P = a point,
    V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
    Return:  wn = the winding number (=0 only when P is outside)
    =#

    # Initialize the winding number:
    winding_number = 0
    n = c.length
    for i in 1:n
        i0 = i
        i1 = i+1
        # (i == length(c)) ? i1 = 1 : i1 = i+1
        if p == c[i]
            return true
        end
        if (c[i0].y <= p.y)
            if (c[i1].y  > p.y)      # an upward crossing
                if (isleft( c[i0], c[i1], p) > 0)  # p left of  edge
                    winding_number += 1            # have  a valid up intersect
                end
            end
        else                         # start y > P.y (no test needed)
            if (c[i1].y  <= p.y) #     // a downward crossing
                if (isleft( c[i0], c[i1], p) < 0)  #// P right of  edge
                    winding_number -= 1    # have  a valid down intersect
                end
            end
        end
    end
    isin = (winding_number != 0)
end
function inshape(c1::ClosedCurve,c2::ClosedCurve)
    #=
    Find the points of c1 which lie in the interior of c2:
    =#

    # We can use the find function from interfaces:
    inds = find(c -> inshape(c,c2),c1)
end

function mightintersect(c1::Curve,c2::Curve)

    m1x = minx(c1)
    M1x = maxx(c1)
    m1y = miny(c1)
    M1y = maxy(c1)
    m2x = minx(c2)
    M2x = maxx(c2)
    m2y = miny(c2)
    M2y = maxy(c2)

    (m1x <= M2x) && (m2x <= M1x) && (m1y <= M2y) && (m2y <= M2x)
    #intervalintersect([m1x,M1x],[m2x,M2x]) & intervalintersect([m1y,M1y],[m2y,M2y])
end

function mightintersect(c1::Segment,c2::Segment)

    m1x = minx(c1)
    M1x = maxx(c1)
    m1y = miny(c1)
    M1y = maxy(c1)
    m2x = minx(c2)
    M2x = maxx(c2)
    m2y = miny(c2)
    M2y = maxy(c2)

    return (m1x <= M2x) && (m2x <= M1x) && (m1y <= M2y) && (m2y <= M1y)
    #intervalintersect([m1x,M1x],[m2x,M2x]) & intervalintersect([m1y,M1y],[m2y,M2y])
end

function orientation(p::Point,q::Point,r::Point)
    #=
    Got this from geeksforgeeks. slides here: http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    0 => p, q and r are colinear
    1 => Clockwise
    -1 => Counterclockwise
    =#
    dtrm = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    if abs(dtrm) < eps()
        dtrm = 0
    end
    round(Int, sign(dtrm))
end

function insegment(p::Point,q::Point,r::Point)
    # Assuming that p,q,r are colinear
    # determine if q is before, in or after p -> r
    (orientation(p,q,r) != 0) && begin
    @show ((q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y))
    error("insegment: p, q, and r are assumed to be colinear")
end
basevec = r + (-1)*p
compvec = q + (-1)*p
dotprod=x(basevec)*x(compvec) + y(basevec) * y(compvec)
if (dotprod < 0)
    -1
elseif dotprod == 0
    0
elseif norm(basevec)^2 >= dotprod
    0
else
    1
end
end

function sort_colinear_points(p::Point,q::Point,r::Point,s::Point)
    # assuming all four points are colinear
    # will sort such that p < r
    pqr = insegment(p,q,r)

    if pqr == 0
        # p < q < r
        qrs = insegment(q,r,s)
        if qrs == 0
            # q < r < s
            p,q,r,s
        elseif qrs == 1
            # q < s < r
            p,q,s,r
        else
            # qrs == -1
            # s < q < r
            # figure out order of psq:
            psq = insegment(p,s,q)
            if psq == -1
                # s < p < q
                s,p,q,r
            elseif psq == 0
                p,s,q,r
            else
                error("sortpoints: contradiction. s > q and s < q")
            end
        end
    elseif pqr == -1
        # then q < p < r
        prs = insegment(p,r,s)
        if prs == 0
            q,p,r,s
        elseif prs == 1
            # s is between p and r
            q,p,s,r
        else
            # prs == -1
            # so  s < p < r
            # and q < p < r
            # so determine relationship between q and s
            # sqp
            sqp = insegment(s,q,p)
            if sqp == 0
                # then s < q < p
                s,q,p,r
            elseif sqp == 1
                # then s < p < q
                s,p,q,r
            else
                # sqp == -1
                # q < s < p
                q,s,p,r
            end
        end
    else
        # pqr == 1
        # so p < r < q
        rqs = insegment(r,q,s)
        if rqs == 0
            p,r,q,s
        elseif rqs == 1
            p,r,s,q
        else
            # rqs == -1
            # s < r < q
            psr = insegment(p,s,r)
            if psr == 0
                p,s,r,q
            elseif psr == -1
                s,p,r,q
            else
                error("sortpoints: contradiction")
            end
        end
    end
end





function onsegment(p::Point,q::Point,r::Point)
    #=
    check if q lies on the line segment from p to r
    Really this just checks if q lives in the box

    =#
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true
    else
        return false
    end
end

function segmentintersect(s1::Segment,s2::Segment; include_ends = true)
    #=
    Assuming s1 and s2 have length 2
    Also assuming genericity for the moment
    returns:
    (bool, s), where bool is true or false and if true s
    is

    =#
    if mightintersect(s1,s2)
        if s1.s == s2.s
            doesintersect = true && include_ends
            p=s1.s
            s = 0
            t = 0
        elseif s1.s == s2.t
            doesintersect = true && include_ends
            p=s1.s
            s = 0
            t = 1
        elseif s1.t == s2.s
            doesintersect = true && include_ends
            p=s1.t
            s = 1
            t = 0
        elseif s1.t == s2.t
            doesintersect = true && include_ends
            p=s1.s
            s = 0
            t = 0
        else
            o1 = orientation(s1.s, s1.t, s2.s)
            o2 = orientation(s1.s, s1.t, s2.t)
            o3 = orientation(s2.s, s2.t, s1.s)
            o4 = orientation(s2.s, s2.t, s1.t)

            # Find the four orientations needed for general and special cases
            # General case
            o = (o1 != o2 && o3 != o4) ? 1 : 0
            if o == 1

                doesintersect = true
                A = [(s1.t.x-s1.s.x) (s2.s.x-s2.t.x)
                (s1.t.y-s1.s.y) (s2.s.y-s2.t.y)]
                w = [s2.s.x-s1.s.x
                s2.s.y-s1.s.y]
                dA = det(A)
                iA = (1/dA)*[A[2,2] -A[1,2]
                -A[2,1]  A[1,1]]
                v = iA*w
                #if (0 <= v[1] <= 1) & (0 <= v[2] <= 1)
                s = v[1]
                t = v[2]
                p = s*s1.t+(1-s)*s1.s

                # else
                #    if o1 == 0 & onsegment(s1.s,s2.s,s1.t)
                #       if onsegment(s1.s, s2.t, s1.s)
                #          doesintersect = 2
                #          # find the point closest to s1.s:
                #          ds = d(s1.s,s2.s)
                #          dt = d(s1.s,s2.t)
                #          if ds < dt
                #          end
                #       end
                #
                #    elseif o1 == 0 & onsegment(s1.s, s2.t, s1.s)
                #
                #
                #       s =
                #         doesintersect = false
                #         p = Point(0.,0.)
                #         s = 0
                #         t = 0
                #       end
            else
                doesintersect = false
                p = Point(0.,0.)
                s = 0
                t = 0
            end
        end

        #elseif

    else

        doesintersect = false
        p = Point(0.,0.)
        s = 0
        t = 0

    end
    return doesintersect, s, t, p
end

# elseif o1 == 1 && onsegment()# special case 1
#
# // p1, q1 and p2 are colinear and p2 lies on segment p1q1
# if (o1 == 0 && onSegment(p1, p2, q1)) return true;
#
# // p1, q1 and p2 are colinear and q2 lies on segment p1q1
# if (o2 == 0 && onSegment(p1, q2, q1)) return true;
#
# // p2, q2 and p1 are colinear and p1 lies on segment p2q2
# if (o3 == 0 && onSegment(p2, p1, q2)) return true;
#
#  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
# if (o4 == 0 && onSegment(p2, q1, q2)) return true;
#
# return false; // Doesn't fall in any of the above cases
#
#

# m1x = minx(s1)
# M1x = maxx(s1)
# m1y = miny(s1)
# M1y = maxy(s1)
# m2x = minx(s2)
# M2x = maxx(s2)
# m2y = miny(s2)
# M2y = maxy(s2)
#
# mayintersect = (m1x <= M2x) & (m2x <= M1x) & (m1y <= M2y) & (m2y <= M2x)
# if mayintersect
#    A = [(s1.t.x-s1.s.x) (s2.s.x-s2.t.x)
#         (s1.t.y-s1.s.y) (s2.s.y-s2.t.y)]
#    dA = det(A)
#    if abs(dA) > .00001
#       iA = (1/dA)*[A[2,2] -A[1,2]
#                   -A[2,1]  A[1,1]]
#       v = iA*[s2.s.x-s1.s.x
#               s2.s.y-s1.s.y]
#       if (0 <= v[1] <= 1) & (0 <= v[2] <= 1)
#          doesintersect = true
#          s = v[1]
#          t = v[2]
#          p = s*s1.t+(1-s)*s1.s
#       else
#          doesintersect = false
#          p = Point(0.,0.)
#          s = 0
#          t = 0
#       end
#    else
#       doesintersect = false
#       p = Point(0.,0.)
#       s = 0
#       t = 0
#    end
# else
#    doesintersect = false
#    p = Point(0.,0.)
#    s = 0
#    t = 0
# end
#    return doesintersect,p,s,t
# end


struct CurveGraph
    curves::Curves
    intersection_points::Vector{Point}
    neighbors::Dict{Int,Vector{Int}}
    edges::Dict{Tuple{Int,Int},Tuple{Int,Tuple{Int,Int}}}
    CurveGraph(crvs::Curves) = curve2curvegraph(crvs)
end

function normalize_curve_segment(c::Curve, ddom::Number=1, mag::Number=1)
    # c should be such that c[1] != c[end]
    # map c[1] -> (0,0)
    # rotate so that c[end] -> (x,0)
    # scale so that c[end] -> (1,0)
    # scale y so that max(maxy(c),-miny(c)) == mag

    # take care of degenerate cases
    (length(c) == 0) && (return c)
    (length(c) == 1) && (return Curve([Point(0,0)]))
    (c[1] == c[end]) && (error("cannot normalize closed curve. should have c[1] != c[end]"))

    c = c + (-1) * c[1]
    dm = c[end] + (-1) * c[1]
    ndm = norm(dm)
    c = (1/ndm) * c
    udm = c[end] + (-1) * c[1]
    ψm = pointθ(udm)
    rψ = [cos(ψm) -sin(ψm); sin(ψm) cos(ψm)]
    # scly = [1 0; 0 mag/max(maxy(c),miny(c))]
    crv = Curve()
    for p in c
        push!(crv,rψ * p)
    end
    c = crv
    # c =  c
    # now scale x and y by ddom and mag
    sclr = [1 0; 0 mag]
    sclr * c
end

function align_curve(motif::Curve,p1::Point,p2::Point)
    # motif assumed to be normalized as in normalize_curve_segment
    dp = p2 + (-1)*p1
    ndp = norm(dp)
    udp = (1/ndp)*dp
    ψp = pointθ(udp)
    rp = [cos(ψp) -sin(ψp); sin(ψp) cos(ψp)]
    mout = Curve()
    for p in motif
        push!(mout, ndp * rp * p + p1)
    end
    # mout = map(p ->  , motif)::Curve
    mout
end



function substitute(motif::Curve, c::Curve, ddom::Number=1, mag::Number=1, offset::Number=0.0)
    (length(c) < 2) && error("substitute input curve must have length >= 2 ")
    if !(0 <= offset <= 1)
        error("substitute: must have 0 <= offset <= 1")
    end
    outcrv = Curve()

    offsetlen = offset * ddom

    prm = parameterize_curve(c)
    # p1 = prm.fcn.prmz1.contents
    prmz2 = prm.fcn.prmz2

    # compute the previous node
    # get the two interval points:
    p2 = prmz2([offsetlen])[1]
    # dc = (c[1] + (-1)*p2)
    # nrmdc = norm(dc)
    # xoffsetlen = nrmdc

    nummotifs = floor(Int,(prm.len - offsetlen + 0.000001) / ddom)
    lastmotift = ddom * nummotifs + offsetlen
    # worry about end points next
    prmcrv = prmz2(collect(linspace(offsetlen,lastmotift,nummotifs+1)))

    # normalize the motif:
    motif = normalize_curve_segment(motif,ddom,mag)

    # if offsetlen > 0
    #
    #     # we know that norm(x(dc)) <= offsetlen
    #     nrmdc > offsetlen && begin
    #         @show x(dc)
    #         @show offsetlen
    #         error("the distance between initial point and start of prmcrv should be offsetlen")
    #     end
    #     ndc = (ddom / nrmdc) * dc
    #     p0 = p2 + ndc
    #
    #     # use normdc rather than offsetlen
    #
    #
    #     # properly aligned curve:
    #     fcrv = align_curve(motif,p0,p2)
    #
    #     # but we might as well do the computations first and then shift:
    #     acrv = motif #align_curve(motif,Point(-ddom,0),Point(0,0))
    #     # we want to go backwards
    #     reverse!(acrv)
    #     reverse!(fcrv)
    #     # Now need to march forward until we go past offsetlen
    #     idx = 1
    # initp = acrv[1]
    #
    # for i in 2:length(acrv)
    #     p = acrv[i]
    #     idx = i
    #     if abs(x(initp) - x(p)) > xoffsetlen
    #         break
    #     end
    # end
    #
    # # so now we know that the point is between idx and idx - 1.
    # fctr = abs(xoffsetlen - x(acrv[idx-1])) / abs(x(acrv[idx]) - x(acrv[idx-1]))
    # !(0 <= fctr <= 1) && begin
    #     @show acrv[idx-1]
    #     @show acrv[idx]
    #     @show fctr
    #     error("bad factor")
    # end
    # # marched forward in reverse so
    # # now fctr is amt to go from idx-1 to idx
    #
    # pnt = (1-fctr)*fcrv[end-idx-1] + fctr * fcrv[end-idx]
    #
    # icrv = pushfirst!(fcrv[end-idx],pnt)
    #
    # # take off the last point
    # outcrv = vcat(outcrv,icrv[1:end-1])
    # else
    #     # I'll need to replace this with better end point handling:
    if !isempty(prmcrv)
        push!(outcrv,prmcrv[1])
        # end
        initp = prmcrv[1]
        for i = 1:nummotifs
            nextp = prmcrv[i+1]
            outcrv = vcat(outcrv,align_curve(motif,initp,nextp)[2:end])
            initp = nextp
        end

        # Now I've got to deal with the last curve:

        # if (outcrv[end] != c[end])
        # the next curve needs to go from outcrv[end] to c[end]

        # end
    end
    outcrv

end





# First let's assume that no more than two curves intersect at a point, and that
# they always intersect in the middle of segments:
# function curve2curvegraph(crvs::Curves)
#
# end
