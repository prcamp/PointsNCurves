
struct spatial_hash
    horizontal_bounds::Tuple{Number,Number}
    vertial_bounds::Tuple{Number,Number}
    horizontal_bin_count::Int
    vertical_bin_count::Int
    curve_to_bin::Dict{Union{ClosedCurve,Curve},Tuple{Int,Int}}
    bin_to_curve::Dict{Tuple{Int,Int},Union{ClosedCurve,Curve}}
end

function spatial_hash(
        horizontal_bounds::Tuple{Number,Number},
        vertical_bounds::Tuple{Number,Number},
        horizontal_bin_count::Int,
        vertical_bin_count::Int,
        curves::Union{ClosedCurves,Curves}
    )::spatial_hash

    curve_to_bin = Dict{Union{ClosedCurve,Curve},Set{Tuple{Int,Int}}}()
    bin_to_curve = Dict{Tuple{Int,Int},Set{Union{ClosedCurve,Curve}}}()
    spatial_hash(
        horizontal_bounds,
        vertical_bounds,
        horizontal_bin_count,
        vertical_bin_count,
        curve_to_bin,
        bin_to_curve
    )
end

function collect_bins
