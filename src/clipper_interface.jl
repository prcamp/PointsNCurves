using Clipper

global sigdigs = 8
global resolution_scale = 2

point2intpoint(p::Point) = IntPoint(p.x,p.y,resolution_scale,sigdigs)
intpoint2point(ip::IntPoint) = begin
    crds = tofloat(ip,resolution_scale,sigdigs)
    Point(crds[1],crds[2])
end
function topath(c::Curve)
        path = Vector{IntPoint}()
        for p in c
            push!(path,point2intpoint(p))
        end
        return path
end
function topath(c::ClosedCurve)
    topath(Curve(c)[1:end-1])
end

function frompath(path::Vector{IntPoint})
    c = Curve()
    for ip in path
        push!(c,intpoint2point(ip))
    end
    c
end
