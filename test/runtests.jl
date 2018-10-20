using PointsNCurves
using Base.Test

function test(run::Function, name::AbstractString; verbose=true)
    test_name = "Test $(name)"
    if verbose
        println(test_name)
    end
    try
        run()
    catch e
        println("Error running: ", test_name)
        rethrow(e)
    end
end

include("unit_tests.jl")