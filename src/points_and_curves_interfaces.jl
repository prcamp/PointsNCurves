import Base.+, Base.*,Base.getindex,Base.convert,
    Base.iterate,Base.setindex!, Base.firstindex,
    Base.lastindex

const P = Float64

struct Point # {T}
   x::P #T
   y::P #T
   Point(x::Number,y::Number)=new(convert(P,x),convert(P,y))
end

Base.isequal(p1::Point,p2::Point) = (p1.x==p2.x)&&(p1.y==p2.y)
norm(p::Point) = sqrt(p.x^2+p.y^2)
struct Polar
  r::P
  θ::P
end
function Polar(p::Point)
  rad = norm(p)
  if rad == 0
    return Polar(0,0)
  else
    theta = pointθ(p)
    return Polar(rad,theta)
  end
end
norm(rt::Polar)=rt.r
polar(p::Point=Point(1,0);center::Point=Point(0.0)) = Polar(p-center)

pointθ(p::Point) = π - (π/2)*(1+sign(p.x))*(1-sign(p.y^2)) - (π/4)*(2+sign(p.x))*sign(p.y) - sign(p.x * p.y) * atan((abs(p.x) - abs(p.y))/(abs(p.x) + abs(p.y)))
dpointθ(p::Point) = - sign(p.x * p.y)/(norm(p)^2) * Point(-sign(p.x)*abs(p.y), sign(p.y) * abs(p.x))

dpointr(p::Point) = (1/norm(p)) * p
Point(rt::Polar) = rt.r*Point(cos(rt.theta),sin(rt.theta))

# Given a function f(x,y) = g(r(x,y),θ(x,y)), we have ∇f = dpolar(p) * ∇g
function dpolar(p::Point)
    dr = dpointr(p)
    dθ = dpointθ(p)
    [dr.x dθ.x; dr.y dθ.y]
end

struct Segment
   s::Point
   t::Point
end

# type ArcInstruct
# struct ArcInstruct
#   params::Vector{AbstractString}
#   code
# end
# type XYArc
# struct XYArc
#   center::Point
#   θ0::Number
#   θ1::Number
# end

const Curve = Vector{Point}

#= the vertices are indexed from 0 to length-1,
   so curve[length]=curve[0]
   this is the only logical thing to do.
=#
struct ClosedCurve
   vertices::Vector{Point}
   length::Int
   ClosedCurve() = new(Vector{Point}(),0)
   ClosedCurve(v::Vector{Point},n::Int) = new(v,n)
   ClosedCurve(v::Vector{Point}) = begin
       if !isempty(v)
           if v[1] == v[end]
               return ClosedCurve(v[1:end-1])
           else
               return new(v,length(v))
           end
       else
           return ClosedCurve()
       end
   end
   ClosedCurve(n::Int) = ClosedCurve(Vector{Point}(undef,n))
end

struct ParameterizedCurve
    curve::Curve
    len::Float64
    lens::Vector{Float64}
    fcn::Function
end

function curve_length(c::Union{Curve,ClosedCurve})#::(Float64,Vector{Float64})
    c = typeof(c) == Curve ? c : Curve(c)
    l = 0.0
    lens = Vector{Float64}()
    if length(c) < 2
        length(c) == 0 ? (l, lens) : (l,push!(lens,l))
    else
        push!(lens, l)
        p0 = c[1]
        for i = 2:length(c)
            p1 = c[i]
            l += norm(p1 + (-1)*p0)
            push!(lens,l)
            p0 = p1
        end
        l, lens
    end
end

function interpolator(lens::Vector{Float64})
    # lens must be increasing
    # output is a function that takes a t, and returns the index of lens s.t.
    # lens[idx] <= t < lens[idx+1]
    # or just the last index if t > lens[end]
    # as well as a float representing where t lies on the segment from
    # lens[idx] to lens[idx+1]
    len = lens[end]
    n = length(lens)
    if n < 2
        error("no point in parameterizing a point")
    end
    interp(t::Float64) = begin
        if t < 0
            1, 0
        elseif t > len
            n-1, 1
        else
            idx1 = findfirst(l -> l > t, lens)
            if idx1 != nothing
                idx0 = idx1-1
                l1 = lens[idx1]
                l0 = lens[idx0]
                idx0, (t - l0)/(l1 - l0)
            else
                n-1, 1
            end
        end
    end
    interp(t::Vector{Float64}) = begin
        indices = Vector{Int}()
        intrps = Vector{Float64}()
        for tt in t
            idx, prm = interp(tt)
            push!(indices,idx)
            push!(intrps, prm)
        end
        indices, intrps
    end
    interp(0.0)
    interp([0.0])
    interp
end

function parameterize_curve(c::Union{Curve,ClosedCurve})
    c = typeof(c) == Curve ? c : Curve(c)
    if length(c) < 2
        if length(c) == 0
            crv(t::Float64) = nothing
            crv(t::Vector{Float64})::Curve = Curve()
            ParameterizedCurve(c, 0, [0] ,crv)
        else
            crv(t::Float64) = c[1]
            crv(t::Vector{Float64})::Curve = map(tt -> c[1], t)
            ParameterizedCurve(c, 0, [0], t -> c)
        end
    end
    len, lens = curve_length(c)
    N = length(c)
    intrp = interpolator(lens)
    intrptr = intrp.interp.contents
    function prmz1(t::Float64)
        idx, prm = intrptr(t)
        if idx == N
            c[end]
        else
                p0 = c[idx]
            p1 = c[idx+1]
            (1-prm)*p0 + prm*p1
        end
    end

    function prmz2(t::Vector{Float64})
        prmc = Curve()
        for tt in t
            idx, prm = intrptr(tt)
            if idx == N
                push!(prmc,c[end])
            else
                p0 = c[idx]
                p1 = c[idx+1]
                push!(prmc,(1-prm)*p0 + prm*p1)
            end
        end
        prmc
    end
    p = t -> typeof(t)==Float64 ? prmz1(t) : prmz2(t)
    ParameterizedCurve(c,len,lens,p)
end

const ClosedCurves = Vector{ClosedCurve}

function Curve(c::Array{Float64,2})
        crv = Curve()
        for i in 1:size(c,1)
            push!(crv, Point(c[i,1],c[i,2]) )
        end
        return crv::Curve
    end

function Curve(c::Array{Int,2})
        crv = Curve()
        for i in 1:size(c)[1]
            push!(crv, Point(c[i,1],c[i,2]) )
        end
        return new(crv)
    end

function Curve(c::ClosedCurve)
    verts = c.vertices::Curve
    if verts[1] == verts[end]
        return verts
    else
        return [verts..., verts[1]]
    end
end

const Curves = Vector{Curve}

# Methods:
Base.iterate(c::ClosedCurve,state=0) =
  if state > c.length
      return nothing
  elseif state == 0
      return (c.vertices[1], 1)
  elseif state == c.length
      return (c.vertices[1], c.length+1)
  else
      return (c.vertices[state+1],state+1)
  end
Base.eltype(::Type{ClosedCurve}) = Point
Base.length(c::ClosedCurve) = c.length
Base.firstindex(c::ClosedCurve) = 0
Base.lastindex(c::ClosedCurve) = c.length - 1
Base.insert!(c::ClosedCurve,idx::Int,p::Point) = begin
  c = Curve(c)
  insert!(c,idx+1,p)
  ClosedCurve(c[1:end-1])
end

using Plots
plt(c::Curve) = plot(x(c),y(c))
plt(c::ClosedCurve) = plt(Curve(c))
plt!(c::Curve) = plot!(x(c),y(c))
plt!(c::ClosedCurve) = plt!(Curve(c))
function plt(C::Vector{Curve})
    plt(C[1])
    for c in C[2:end-1]
        plt!(c)
    end
    plt!(C[end])
end
plt(C::Vector{ClosedCurve}) = plt(map(c -> Curve(c),C))
function plt!(C::Vector{Curve})
    plt!(C[1])
    for c in C[2:end-1]
        plt!(c)
    end
    plt!(C[end])
end
plt!(C::Vector{ClosedCurve}) = plt!(map(c -> Curve(c),C))
pltface(c::ClosedCurve) = plot(Shape(x(Curve(c)),y(Curve(c))),opacity=.8)
pltface!(c::ClosedCurve) = plot!(Shape(x(Curve(c)),y(Curve(c))),opacity=.8)
function pltface(C::Vector{ClosedCurve})
    if isempty(C); println("nothing to do"); return; end
    pltface(C[1])
    pltface!(C[2:end])
end

function pltface!(C::Vector{ClosedCurve})
    for c in C[1:end-1]
        pltface!(c)
    end
    pltface!(C[end])
end

scttr(p::Point,sym::Symbol=:+) = scatter([x(p)],[y(p)],marker=sym)
scttr!(p::Point,sym::Symbol=:+) = scatter!([x(p)],[y(p)],marker=sym)
function scttr!(c::Curve,sym::Symbol=:+)
    for p in c[1:end-1]
        scttr!(p,sym)
    end
    scttr!(c[end],sym)
end
function scttr(c::Curve,sym::Symbol=:+)
    scttr(c[1],sym)
    if length(c) > 1
        scttr!(c[2:end],sym)
    end
end


function cleancurve(crvs::Curves)
  ncrvs = Curves()
  for c in crvs
    push!(ncrvs,cleancurve(c))
  end
  return ncrvs
end
function cleancurve(c::Curve)
  crv = Curve()
  v0 = c[1]
  push!(crv,v0)
  for i in 2:length(c)
    if c[i] != v0
      push!(crv,c[i])
      v0 = c[i]
    end
  end
  return crv
end

# Indexing:
Base.getindex(c::ClosedCurve,i::Int) = begin
   if 0 <= i <= c.length-1
      return c.vertices[i+1]
   else
      i = mod(i,c.length)+1
      return c.vertices[i]
   end
end
Base.getindex(c::ClosedCurve, i::Number) = c[convert(Int, i)]
function Base.getindex(c::ClosedCurve, I)
   if length(I) > 2
      return ClosedCurve([c[i] for i in I])
   elseif length(I)==2
      return Segment(c[I[1]],c[I[2]])
   end
end
mean(x::Vector{Float64}) = sum(x)/length(x)
mean(c::ClosedCurve) = Point(mean(x(c)),mean(y(c)))
import Base.size
Base.size(c::ClosedCurve) = (c.length,)
# Base.setindex
# Base.endof

using Clipper
include("clipper_interface.jl")

import Clipper.area
area(c::ClosedCurve) = abs(area(topath(Curve(c))))/10^8

function closedcurve(c::Array{P,2})
   # convert a n x 2 array into a curve
   crv = Point[]
   for i in 1:size(c)[1]
      push!(crv, Point(c[i,1],c[i,2]) )
   end
   return ClosedCurve(crv)
end

function Base.convert(::Type{Curve}, crv::ClosedCurve)
  c = crv.vertices
  (!isequal(c[end],c[1])) && (push!(c,c[1]))
  return c
end

function curves2array(crvs::Curves)::Array{Array{Number,2},1}
  arr = Array{Number,2}[]
  for c in crvs
    push!(arr,hcat(x(c),y(c)))
  end
  return arr
end

function d(p1::Point, p2::Point)
   norm(p1+(-1*p2))
end

# Should also do a curve plus a point, and a scalar times a curve,
# and a curve + a curve
x(p::Point) = p.x
y(p::Point) = p.y
x(c::Curve) = map(p -> x(p), c)
y(c::Curve) = map(p -> y(p), c)
x(c::ClosedCurve) = map(p -> x(p), c.vertices)
y(c::ClosedCurve) = map(p -> y(p), c.vertices)
x(s::Segment) = [s.s.x,s.t.x]
y(s::Segment) = [s.s.y,s.t.y]

function point2vec(p::Point)
  Float64[x(p),y(p)]
end

minx(c::Curve) = minimum(x(c))#mapreduce(x,min,c)
maxx(c::Curve) = maximum(x(c))#mapreduce(x,max,c)
miny(c::Curve) = minimum(y(c))#mapreduce(y,min,c)
maxy(c::Curve) = maximum(y(c))#mapreduce(y,max,c)
minx(crvs::Curves) = mapreduce(c->minx(c),min,crvs)
maxx(crvs::Curves) = mapreduce(c->maxx(c),max,crvs)
miny(crvs::Curves) = mapreduce(c->miny(c),min,crvs)
maxy(crvs::Curves) = mapreduce(c->maxy(c),max,crvs)
minx(c::ClosedCurve) = minimum(x(c))#mapreduce(x,min,c)
maxx(c::ClosedCurve) = maximum(x(c))#mapreduce(x,max,c)
miny(c::ClosedCurve) = minimum(y(c))#mapreduce(y,min,c)
maxy(c::ClosedCurve) = maximum(y(c))#mapreduce(y,max,c)
minx(crvs::ClosedCurves) = mapreduce(c->minx(c),min,crvs)
maxx(crvs::ClosedCurves) = mapreduce(c->maxx(c),max,crvs)
miny(crvs::ClosedCurves) = mapreduce(c->miny(c),min,crvs)
maxy(crvs::ClosedCurves) = mapreduce(c->maxy(c),max,crvs)
minx(s::Segment) = minimum(x(s))
miny(s::Segment) = minimum(y(s))
maxx(s::Segment) = maximum(x(s))
maxy(s::Segment) = maximum(y(s))

+(p1::Point,p2::Point) = Point(p1.x+p2.x,p1.y+p2.y)
function +(crv::ClosedCurve, pnt::Point)
   ncrv = Point[]
   for p in crv
      push!(ncrv,p+pnt)
   end
   ClosedCurve(ncrv,crv.length)
end
function +(crv::Curve, pnt::Point)
   ncrv = Point[]
   for p in crv
      push!(ncrv,p+pnt)
   end
   return ncrv
end
+(pnt::Point, crv::ClosedCurve) = +(crv::ClosedCurve, pnt::Point)
function +(c1::ClosedCurve,c2::ClosedCurve)
   if c1.length != c2.length
      error("To add two curves they must be the same length")
   else
      crv = Point[]
      for i in 0:c1.length-1
         push!(crv,c1[i]+c2[i])
      end
   end
   ClosedCurve(crv)
end
function +(crvs::Union{Curves,ClosedCurves},p::Point)
  crvs2 = typeof(crvs)()
  for c in crvs
    c = c + p
    push!(crvs2,c)
  end
  return crvs2
end

*(a::Number,p::Point) = Point(a*p.x,a*p.y)
function *(a::Number, crv::ClosedCurve)
   ncrv = Point[]
   for p in crv
      push!(ncrv,a*p)
   end
   ClosedCurve(ncrv)
end
function *(a::Number,s::Segment)
  p = (1-a)*s.s + a*s.t
end
*(A::Array{T,2},p::Point) where T <: Number = Point((A*point2vec(p))...)
function *(A::Array{T,2},c::Union{Curve,ClosedCurve}) where T <: Number
  Tp = typeof(c)
  crv = Curve()
  for p in c
    push!(crv,A*p)
  end
  return Tp(crv)
end
function *(A::Array{T,2},c::Union{Curves,ClosedCurves}) where T <: Number
  crvs = similar(c)
  for (i,crv) in enumerate(c)
    crvs[i]=A*crv
  end
  return crvs
end



function genrandpoint(p0::Point=Point(-1,-1),p1::Point=Point(1,1),n::Int=1)
    xx = x([p0,p1])
    yy = y([p0,p1])
    mx = minimum(xx)
    Mx = maximum(xx)
    my = minimum(yy)
    My = maximum(yy)
    if n == 1
        px = (Mx-mx)*rand()+mx
        py = (My-my)*rand()+my
        Point(px,py)
    else
        px = (Mx-mx).*rand(n).+mx
        py = (My-my).*rand(n).+my
        Curve([px py])
    end
end

function genrandcrv(n::Int = 25)
   x = sort(rand(n))
   y = rand(n)
   crv = hcat(x,y)
   crv = curve(crv)
end

function genrandcirc(n::Int = 10)
   theta = collect(range(0,2*pi;length=n+1)[1:n]).+2*pi*rand()
   scl = 2*rand()
   x = cos.(theta)
   y = sin.(theta)
   crv = hcat(x,y)
   crv = closedcurve(crv)
   crv = scl*crv + genrandpoint()
end

function gencirc(n::Int=10)
   theta = collect(range(0,2*pi;length=n+1)[1:n]).+2*pi*rand()
   x = cos.(theta)
   y = sin.(theta)
   crv = closedcurve(hcat(x,y))
end

function genstar(n=5,theta0 = 2*pi*rand())
  theta = collect(range(0,2*pi;length=n+1)[1:n]).+theta0
  halfind = round(Int,ceil(n/2))
  inds = vcat(2*(1:halfind).-1,2*(1:halfind-1))
  theta = theta[inds]
  x = cos.(theta)
  y = sin.(theta)
  crv = closedcurve(hcat(x,y))
end

function gensquare(width::Number=1,center::Point=Point(0,0))
  square = .5*ClosedCurve([Point(-1,-1)
                           Point(-1,1)
                           Point(1,1)
                           Point(1,-1)])
  square = width*square + center
end

function gentestsquare()
  crv = Curve()
  push!(crv,Point(0,0))
  push!(crv,Point(0,1))
  push!(crv,Point(1,1))
  push!(crv,Point(1,0))
  push!(crv,Point(0.25,0))
  push!(crv,Point(.5,.5))
  push!(crv,Point(0.75,0))
  crv = ClosedCurve(crv)
end

linspace(a,b,n) = range(a,b; length = n)
