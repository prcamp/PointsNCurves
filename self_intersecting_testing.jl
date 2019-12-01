include("src/points_and_curves_interfaces.jl")
include("src/points_and_curves_utils.jl")
include("src/p&c_recursors.jl")
include("src/MakeSVGDrawingCmds.jl")

fn = "star2.svg"
scale = 1000
# center =

function gensimplestar(numpoints = 5, phase=0)
    angles = range(0+phase,2*pi+phase;length=2*numpoints+1)
    c = Curve()
    innerrad = .25
    outerrad = 1
    for i in 1:length(angles)-1
        θ = angles[i]
        if iseven(i)
            push!(c,innerrad*Point(cos(θ),sin(θ)))
        else
            push!(c,outerrad*Point(cos(θ),sin(θ)))
        end
    end
    ClosedCurve(c)
end

function genstars(numstars=10,numpnts=5,pow=1,center_variance=0)
    C = ClosedCurves()
    scales = range(1,scale^(1/pow); length=numstars).^pow
    cents = range(0, 0; length=numstars)
    for (c,scl) in zip(cents,scales)
        c = sqrt(center_variance).*randn(2)
        p = Point(c[1],c[2])
        push!(C, (scl*gensimplestar(numpnts,pi*rand())+ p))
    end
    return C
end

function genstarseq(numstars=10,numpnts=5,pow=1,center_variance=0)
    C = genstars(numstars,numpnts,pow,center_variance)

    @time crvs = occlude_sequence(C)
    # m,o = masked_occluded(C[1],C[2])
    return crvs
end
using Clipper
include("clipper_interface.jl")

function occlude_sequence_clpr(crvs::Vector{Vector{IntPoint}})
    out = Vector{Vector{IntPoint}}()
    # paths = map(c -> topath(c),crvs)
    # clp = Clip()
    push!(out,crvs[1])
    # add_path!(clp,crvs[1],PolyTypeClip,true)

    for i = 2:length(crvs)
        clp = Clip()
        add_paths!(clp,reverse(crvs[1:i-1]),PolyTypeClip,true)
        add_path!(clp, crvs[i],PolyTypeSubject,true)
        #
        # clp2 = clp
        # add_path!(clp2,crvs[i],PolyTypeSubject,true)
        res, os = execute(clp, ClipTypeDifference, PolyFillTypeEvenOdd, PolyFillTypeNonZero)
        if res == false
            @warn "oh no res was false!"
        end
        out = vcat(out, os)
        # add_path!(clp,crvs[i],PolyTypeClip,true)
    end
    return out
end

function genstarseq_clpr(numstars=10,numpnts=5,pow=1,center_variance=0)
    C = genstars(numstars,numpnts,pow,center_variance)
    Strs = map(c -> topath(c),C)
    @time stseq = occlude_sequence_clpr(Strs)
    crvs = map(c -> frompath(c),stseq)
    # m,o = masked_occluded(C[1],C[2])
    return crvs
end
global i = 1
function genstardrwng(numstars=10, numpnts=3,pow=1, center_var=0)

    crvs = genstarseq(numstars,numpnts,pow,center_var)
    DrawSVGCurves(fn,crvs)
    crvs
end


function smallest_step(c::ClosedCurve)
    smallest_step(Curve(c))
end
function smallest_step(c::Curve)
    if isempty(c)
        return 10000000
    elseif length(c) == 1
        return 0
    else
        mindist = d(c[1],c[2])
        for i in 2:length(c)
            mindist = min(mindist,d(c[i],c[i-1]))
        end
        return mindist
    end
end

function smallest_step(c::Curves)
    if isempty(c)
        return 10000000
    else
        mindist = smallest_step(c[1])
        for i in 2:length(c)
            mindist = min(smallest_step(c[i]),mindist)
        end
        return mindist
    end
end

function smallest_step(c::ClosedCurves)
    if isempty(c)
        return 10000000
    else
        mindist = smallest_step(c[1])
        for i in 2:length(c)
            mindist = min(smallest_step(c[i]),mindist)
        end
        return mindist
    end
end



@time genstardrwng(30,3,.6);
#@time genstardrwng(300,20,.6)
# @time genstardrwng(400,20,.6)
#@time genstardrwng(100,40,.6)
# @time genstardrwng(60,80,.6)
# @time genstardrwng(600,20,.6)
# @time genstardrwng(100,3,.6,1000)
# @time genstardrwng(1000,3,.6,1000)
# best:
# @time genstardrwng(10000,3,.6,1000)

# @time genstardrwng(1000,6,.6,1000)
# @time genstardrwng(1000,5,.6,1000)
# @time genstardrwng(2000,6,.6,1000)

# @time genstardrwng(2000,6,.6,10000)

# @time genstardrwng(4000,6,.6,10000)
