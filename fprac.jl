# fprac.jl
include("src/points_and_curves_interfaces.jl")
include("src/points_and_curves_utils.jl")
include("src/p&c_recursors.jl")
include("src/MakeSVGDrawingCmds.jl")
w = 100
sw = 300

function differential(c::Curve)
  if length(c) < 2
    Curve([])
  else
    d = Curve([])
    for i in 2:length(c)
      push!(d,c[i]+(-1*c[i-1]))
    end
    d
  end
end

function circle(θ1,θ2,rad=1.5*w,dir=1)
  function circle_points(n::Int)
    θ = range(θ1,θ2;length=n+2)[2:end-1]
    pnts = map(t -> rad*Point(cos(t),sin(t)),θ)
    if dir > 0
      return differential(pnts)
    else
      return differential(reverse(pnts))
    end
  end
end

function outer_circle(θ1,θ2,dir=1)
  circle(θ1,θ2,1.5*w,dir)
end
function inner_circle(θ1,θ2,dir=1)
  circle(θ1,θ2,.5*w,dir)
end

add_by!(c::Curve,incr::Point) = push!(c,c[end]+incr)

#letter_f =
begin
  c = Curve([Point(0,0)])
  add_by!(c,Point(0,sw-w))
  add_by!(c,Point(-w,0))
  add_by!(c,Point(0,w))
  add_by!(c,Point(w,0))
  add_by!(c,Point(0,.5*sw))
  oc = outer_circle(π,0,1)(20)
  for p in oc
    add_by!(c,p)
  end
  add_by!(c,Point(-w,0))
  ic = inner_circle(0,π)(18)
  for p in ic
    add_by!(c,p)
  end
  add_by!(c,Point(0,-.5*sw))
  add_by!(c,Point(w,0))
  add_by!(c,Point(0,-w))
  add_by!(c,Point(-w,0))
  add_by!(c,Point(0,-(sw-w)))
  letter_f = ClosedCurve(c)

  c = Curve([Point(2.1*w,sw)])
  add_by!(c,Point(w,0))
  oc = outer_circle(π/2,-π/2)(20)
  for p in oc
    add_by!(c,p)
  end
  add_by!(c,Point(0,-sw))
  add_by!(c,Point(-w,0))
  letter_p = ClosedCurve(c)

  c = Curve([Point(4.5*w,0)])
  add_by!(c,Point(0,.5*sw))
  oc = outer_circle(π,0,1)(20)
  for p in oc
    add_by!(c,p)
  end
  add_by!(c,Point(-w,0))
  ic = inner_circle(0,π)(18)
  for p in ic
    add_by!(c,p)
  end
  add_by!(c,Point(0,-.5*sw))
  letter_r = ClosedCurve(c)

  c = Curve([Point(8.8*w,0)])
  oc = outer_circle(3*π/2,π/2,1)(20)
  for p in oc
    add_by!(c,p)
  end
  add_by!(c,Point(w,0))
  add_by!(c,Point(0,-3*w))
  letter_a = ClosedCurve(c)

  c = Curve([Point(12.8*w,2*sqrt(2)*w)])
  oc = outer_circle(π/4,7*π/4,1)(80)
  for p in oc
    add_by!(c,p)
  end
  add_by!(c,(1/sqrt(2))*Point(-w,w))
  ic = inner_circle(7*π/4,π/4)(80)
  for p in ic
    add_by!(c,p)
  end
  add_by!(c,(1/sqrt(2))*Point(w,w))
  letter_c = ClosedCurve(c)


  DrawSVGCurves("test_fprac.svg",[letter_f,letter_p,letter_r,letter_a,letter_c])
end
