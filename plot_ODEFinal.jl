include("ODEFinal.jl")
include("helper_functions.jl")
include("visualConventions.jl")
using DifferentialEquations
using Plots
using DataFrames

# each ODE has 5 variables (none dependent directly on time) and 16 parameters

curr_init = 5

# plots a single run of a fully specified simulation
function plotRun(model, params, init, timespan)
    problem = ODEProblem(model, init, timespan, params)
    solution = solve(problem)
    plot(solution)
end

# plots 5 heat maps, 3 for the raw links, and 1 for the combined 
function LiaoTypeHeatMap(model, params, initReduced)
end

# plots 5 boundary maps, 3 for the raw links, and 1 for the combined 
function LiaoTypeBoundaryMap(model, params, initReduced, persistenceThreshold)
end
