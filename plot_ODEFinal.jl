include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))
include(joinpath(@__DIR__, "CONSTANTS", "visualConventions.jl"))
using DifferentialEquations
using Plots
using DataFrames

# each ODE has 5 variables (none dependent directly on time) and 16 parameters

# plots a single run of a fully specified simulation
function plotRun(model, params, init, timespan)
    problem = ODEProblem(model, init, timespan, params)
    solution = solve(problem)
    plot(solution)
end

# simulate model accross grid and return data 
# NOTE: <baseParams> consists of the first 11 nuisance parameters
# NOTE: <resourceParams> are the last 3 parameters that specify the resource life-history
function LiaoTypeGrid(model, baseParams, resourceParams, grain, init, timespan)
    # create grid and empty dataframe
    grid = create_unit_grid(grain)
    data = DataFrame(
        Availability=Float64[], Connectivity=Float64[],
        DensityR=Float64[], Density2=Float64[], Density3=Float64[],
        Density23=Float64[], Density13=Float64[]
    )
    # populate dataframe
    for point in grid
        if !is_valid(point)
            push!(data, (point[1], point[2], NaN, NaN, NaN, NaN, NaN))
        else
            landscapeParams = convert_landscape_params([point[1], point[2]])
            curr_params = [baseParams; landscapeParams; resourceParams]
            curr_problem = ODEProblem(model, init, timespan, curr_params)
            curr_solution = solve(curr_problem)
            sol_end = curr_solution[end]
            push!(data, (point[1], point[1], sol_end[1], sol_end[2], sol_end[3] + sol_end[4], sol_end[3], sol_end[4]))
        end
    end
    return (data)
end

# plots 5 heat maps, 4 for the raw links + 1 for the combined top predator density
function LiaoTypeHeatMap(model, baseParams, resourceParams, grain, init, timespan)
    # get simGridData
    simGridData = LiaoTypeGrid(model, baseParams, resourceParams, grain, init, timespan)
    # setup
    grid_length = 0:grain:1
    n = length(grid_length)
    Z_R = reshape(simGridData.DensityR, n, n)
    Z_2 = reshape(simGridData.Density2, n, n)
    Z_3 = reshape(simGridData.Density3, n, n)
    Z_23 = reshape(simGridData.Density23, n, n)
    Z_13 = reshape(simGridData.Density13, n, n)
    # create heatmaps
    pR = heatmap(grid_length, grid_length, Z_R, title="Resource Density (P₁)",
        aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
    p2 = heatmap(grid_length, grid_length, Z_2, title="Consumer Density (P₂)",
        aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
    p3 = heatmap(grid_length, grid_length, Z_3, title="Tertiary Density (P₃)",
        aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
    p23 = heatmap(grid_length, grid_length, Z_23, title="2-3 Link Density (P₍₂₃₎)",
        aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
    p13 = heatmap(grid_length, grid_length, Z_13, title="1-3 Link Density (P₍₁₃₎)",
        aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
    # arrange plots
    plot(pR, p2, p3, p23, p13, size=(1600, 1000))
end



# plots 5 boundary maps, 3 for the raw links, and 1 for the combined 
#function LiaoTypeBoundaryMap(model, params, initReduced, persistenceThreshold)
#end


#########################################################################################################################
# TESTING
#########################################################################################################################
