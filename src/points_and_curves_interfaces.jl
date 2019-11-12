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
    if p.x != 0
      theta = 2*atan(p.y/p.x)
      return Polar(rad,theta)
    elseif sign(p.y)>0
      return Polar(rad,pi/2)
    else
      return Polar(rad,-pi/2)
    end
  end
end
norm(rt::Polar)=rt.r
polar(p::Point=Point(1,0);center::Point=Point(0.0)) = Polar(p-center)

Point(rt::Polar) = rt.r*Point(cos(rt.theta),sin(rt.theta))

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

function Curve(c::Array{Float64,2})
        crv = Curve()
        for i in 1:length(c)
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
    c.vertices::Curve
end

const Curves = Vector{Curve}

#= the vertices are indexed from 0 to length-1,
   so curve[length]=curve[0]
   this is the only logical thing to do.
=#
struct ClosedCurve
   vertices::Vector{Point}
   length::Int
   ClosedCurve(v::Vector{Point},n::Int) = new(v,n)
   ClosedCurve(v::Vector{Point}) = new(v,length(v))
end
const ClosedCurves = Vector{ClosedCurve}

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

function curves2array(crvs::Curves)
  arr = Any[]
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
      for i in 1:c1.length
         push!(crv,c1[i]+c2[i])
      end
   end
   ClosedCurve(crv)
end
function +(crvs::Curves,p::Point)
  crvs2 = Curves()
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



function genrandpoint()
   x = 2*rand()-1
   y = 2*rand()-1
   Point(x,y)
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

function genstar(n,theta0 = 2*pi*rand())
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
