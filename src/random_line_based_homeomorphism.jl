
# generate a random homeomorphism taking [-1,1]^2 -> R^2
# using lines
include("points_and_curves_interfaces.jl")
function generate_random_lines()
    width = 1/4
    θ = 2*π*rand()
    v = Point(cos(θ),sin(θ))
    r = () -> 2*rand()-1
    p1 = Point(r(),r())
    p2 = width*(1/norm(p1))*Point(-y(p1),x(p1)) + p1
    p1, p2, v
end

# gives distance from the point p to the line p1 + t v
function dist_to_line(p::Point,p1::Point,v::Point)
    p1top = p1 + (-1)*p
    sclr = -(x(p1top)*x(v) + y(p1top) * y(v))
    norm(p1top + sclr * v)
end

function generate_random_betweeness_function()
    p1,p2,v = generate_random_lines()
    generate_random_betweeness_function(p1,p2,v)
end

function generate_random_betweeness_function(p1::Point,p2::Point,v::Point)
    dl = dist_to_line(p1,p2,v)
    s = (p::Point) -> begin
        dp1 = dist_to_line(p,p1,v)
        dp2 = dist_to_line(p,p2,v)
        if (dp1 < dl) && (dp2 < dl)
            2*min(dp1,dp2)
        else
            0
        end
    end
end

function generate_random_purturbation()
    p1,p2,v = generate_random_lines()
    scl = 5
    s = generate_random_betweeness_function(p1,p2,v)
    (p::Point) -> p + scl * s(p) * v
end


function generate_random_homeomorphism(n::Int)
    homs = []
    for i in 1:n
        purt = generate_random_purturbation()
        push!(homs,purt)
    end

    (p::Point) -> begin
        p0 = Point(0,0)
        for prt in homs
            p0 += (1/n)*prt(p)
        end
        p + p0
    end
end
