mutable struct RandomSurface
    sfreq::Union{Int,Nothing}
    tfreq::Union{Int,Nothing}
    coeffs::Union{Vector{Vector{Float64}},Nothing}
    f::Union{Function,Nothing}
    gf::Union{Function,Nothing}
    RandomSurface() = new(nothing,nothing,nothing,nothing,nothing)
end

triangle(t::Number) =
  (abs(t) < 0.00001 || abs(t-pi) < 0.00001 || abs(t-2*pi) < 0.00001) ? 0 : (0 <= t <= pi ? (t <= pi/2 ? t/(pi/2) : -t/(pi/2) + 2) : (abs(t - pi) < 0.00001 ? 0 : (pi <= t <= 2*pi ? (t <= 3*pi/2 ? -t/(pi/2) + 2 : t/(pi/2) - 4 ) : triangle(mod2pi(t)))))

square(t::Number) =
  (abs(t) < 0.00001 || abs(t-pi) < 0.00001 || abs(t-2*pi) < 0.00001) ? 0 : (0 <= t <= pi ? 1 : (abs(t - pi) < 0.00001 ? 0 : (pi <= t <= 2*pi ? -1 : square(mod2pi(t)))))

#
# mutable struct RandomLoop
#     freq::Union{Int,Nothing}
#     xcoeffs::Union{Vector{Float64},Nothing}
#     ycoeffs::Union{Vector{Float64},Nothing}
#     f::Union{Function,Nothing}
#     RandomLoop() = new(nothing,nothing,nothing,nothing)
# end
#

function RandomSurface(sfreq::Int,tfreq::Int)::RandomSurface
    rs = RandomSurface()
    rs.sfreq = sfreq
    rs.tfreq = tfreq
    rs.coeffs = Vector{Vector{Float64}}(undef,rs.tfreq)
    for i in 1:rs.tfreq
        rs.coeffs[i] = randn(rs.sfreq)
    end
    rs.f = generate_random_surface(rs.coeffs)
    rs.gf = generate_random_surface_gradient(rs.fctl.coeffs)
    rs
end

function generate_random_function(coeffs::Vector{Float64}; fcn=sin)
    t -> mapreduce(i -> (coeffs[i]/i)*fcn.(i*π.*t),+,1:length(coeffs))
end

function generate_random_function(n::Int; fcn = sin)
    α = randn(n)
    t -> mapreduce(i -> (α[i]/i)*fcn.(i*π.*t),+,1:n)
end

function generate_random_curve(coeffsx::Vector{Float64},coeffsy::Vector{Float64}; fcn = sin)
  cx = generate_random_function(coeffsx; fcn = fcn)
  cy = generate_random_function(coeffsy; fcn = fcn)
  c(t) = map(s -> Point(cx(s),cy(s)),t)
  c
end

function generate_random_surface(f::Vector{Function}; fcn = sin)
    (s,t) -> mapreduce(jf -> (jf[2](s)/jf[1])*fcn.(jf[1]*π.*t),+,enumerate(f))
end

function generate_random_surface(s_coeff_count::Int,t_coeff_count::Int; sfcn = sin, tfcn = sin)
    coeffs = Vector{Function}(undef,t_coeff_count)
    for i in 1:t_coeff_count
        coeffs[i] = generate_random_function(s_coeff_count; fcn = sfcn)
    end
    generate_random_surface(coeffs; fcn = tfcn)
end

function generate_random_surface(coeffs::Vector{Vector{Float64}}; sfcn = sin, tfcn = sin)
    t_coeff_count = length(coeffs)
    coeffuncs = Vector{Function}(undef,t_coeff_count)
    for i in 1:t_coeff_count
        coeffuncs[i] = generate_random_function(coeffs[i]; fcn = sfcn)
    end
    generate_random_surface(coeffuncs; fcn = tfcn)
end
