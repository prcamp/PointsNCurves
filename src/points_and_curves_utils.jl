
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
    const tt = 0.5
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
  const t = .99
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
          println("does intersect: \n",
          "s is $s \n",
          "t is $t \n",
          "p is ($(r.x),$(r.y)) \n",
          )
          push!(intersectingsegments,j)
          # ignoring double intersections for the moment
          push!(intersectionpoints,r)
          push!(scale,s)
        end
      end
    end
    println("now to deal with those intersections...")
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

  (m1x <= M2x) & (m2x <= M1x) & (m1y <= M2y) & (m2y <= M2x)
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
  o = round(Int, sign((q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)))
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

function segmentintersect(s1::Segment,s2::Segment)
  #=
  Assuming s1 and s2 have length 2
  Also assuming genericity for the moment
  returns:
  (bool, s), where bool is true or false and if true s
  is

  =#
  if mightintersect(s1,s2)
    if s1.s == s2.s
      doesintersect = true
      p=s1.s
      s = 0
      t = 0
    elseif s1.s == s2.t
      doesintersect = true
      p=s1.s
      s = 0
      t = 1
    elseif s1.t == s2.s
      doesintersect = true
      p=s1.t
      s = 1
      t = 0
    elseif s1.t == s2.t
      doesintersect = true
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
