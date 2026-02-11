include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))
include(joinpath(@__DIR__, "CONSTANTS", "visualConventions.jl"))
using Revise
using DifferentialEquations
using Plots
using DataFrames

# Note: each ODE has 6 variables (none dependent directly on time) and 16 parameters
# Note: of a coexistance pattern, proportional persistance represents the area of the Liao Type plot it occupies

################################################################################################################
# DATA GENERATING FUNCTIONS
################################################################################################################

# simulate model accross grid and return data 
# NOTE: <baseParams> consists of the first 11 nuisance parameters
# NOTE: <resourceParams> are the last 3 parameters that specify the resource life-history
function LiaoTypeGrid(model, baseParams, resourceParams, grain, init, timespan)
    # input check
    if !LiaoTypeGridValidInput(baseParams, resourceParams, grain, init, timespan)
        error("{LiaoTypeGrid}: Invalid input.")
    end
    # create grid and empty dataframe
    grid = create_unit_grid(grain)
    data = DataFrame(
        Availability=Float64[], Connectivity=Float64[],
        DensityR=Float64[], Density2=Float64[], Density3=Float64[],
        Density23=Float64[], Density13=Float64[]
    )
    # populate dataframe
    for point in grid
        if !is_gridPoint_valid(point)
            push!(data, (point[1], point[2], NaN, NaN, NaN, NaN, NaN))
        else
            UF = transform_goodLandscape_params([point[1], point[2]])
            curr_params = [baseParams; UF; resourceParams]
            curr_problem = ODEProblem(model, init, timespan, curr_params)
            curr_solution = solve(curr_problem)
            sol_end = curr_solution[end]
            # if the end solution is mathematically invalid, throw error
            if !is_valid_ODE_end(point, sol_end)
                error("{LiaoTypeGrid}: ODE solution violates model variable constraints.")
            end
            push!(data, (point[1], point[1], sol_end[1], sol_end[2], sol_end[3] + sol_end[4], sol_end[3], sol_end[4]))
        end
    end
    return (data)
end

# runs LiaoTypeGrid, then checks if the links are persisting at a given point
function LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan)
    # get raw grid data simGridData + threshold
    rawGridData = LiaoTypeGrid(model, baseParams, resourceParams, grain, init, timespan)
    n = nrow(rawGridData)
    t = PERSISTENCE_THRESHOLD
    # return a dataframe augmented with binary persistence data
    newData = DataFrame(
        Availability=rawGridData.Availability, Connectivity=rawGridData.Connectivity,
        DensityR=rawGridData.DensityR, Density2=rawGridData.Density2, Density3=rawGridData.Density3,
        Density23=rawGridData.Density23, Density13=rawGridData.Density13,
        PersistenceR=discretizeDensity(rawGridData.DensityR, t),
        Persistence2=discretizeDensity(rawGridData.Density2, t),
        Persistence3=discretizeDensity(rawGridData.Density3, t),
        Persistence23=discretizeDensity(rawGridData.Density23, t),
        Persistence13=discretizeDensity(rawGridData.Density13, t))
    return newData
end

# generates a single point vector, to be plotted in <ProportionalPersistencePlot>
function ProportionalPersistencePoint(model, baseParams, resourceParams, grain, init, timespan)
    # get augmented sim data, add distribution IDs, and save it's number of rows
    simGridData = LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan)
    simGridData.SpeciesDistributionID = assignSpeciesDistributionID(simGridData)
    n_valid_landscapes = count(!isnan, simGridData.SpeciesDistributionID)
    # for each unique distribution ID, caclulate its proportional persistence
    pointVector = zeros(length(SPECIES_DISTRIBUTION_IDS))
    for i in eachindex(SPECIES_DISTRIBUTION_IDS)
        curr_ID = SPECIES_DISTRIBUTION_IDS[i]
        pointVector[i] = count(simGridData.SpeciesDistributionID .== curr_ID) / n_valid_landscapes
    end
    # output check
    if abs(sum(pointVector) - 1) > grain^2 # threshold chosen arbitrarily
        error("{ProportionalPersistencePoint}: Point vector sum is very different from 1")
    end
    return pointVector
end

# generates a single point vector, to be plotted in <WeightedProportionalPersistencePlot>
# Note: this function weighs the proportion of each coexistance pattern by the mean equilibrium resource density accross it's range
function WeightedProportionalPersistencePoint(model, baseParams, resourceParams, grain, init, timespan)
    # get augmented sim data, add distribution IDs, and save it's number of rows
    simGridData = LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan)
    simGridData.SpeciesDistributionID = assignSpeciesDistributionID(simGridData)
    n_valid_landscapes = count(!isnan, simGridData.SpeciesDistributionID)
    # for each unique distribution ID, caclulate its proportional persistence
    weightedPointVector = zeros(length(SPECIES_DISTRIBUTION_IDS))
    for i in eachindex(SPECIES_DISTRIBUTION_IDS)
        curr_ID = SPECIES_DISTRIBUTION_IDS[i]
        filtered_data = simGridData[simGridData.SpeciesDistributionID.==curr_ID, :]
        # handle if this species coexistance pattern does not exist in the given data
        if nrow(filtered_data) == 0
            weightedPointVector[i] = 0.0
        else
            weightedPointVector[i] = sum(filtered_data.DensityR) / n_valid_landscapes
        end
    end
    return weightedPointVector
end

################################################################################################################
# GENUINE PLOTTING FUNCTIONS
################################################################################################################

# plots a single run of a fully specified simulation
function plotRun(model, params, init, timespan)
    # input check
    if plotRunValidInput(params, init, timespan)
        problem = ODEProblem(model, init, timespan, params)
        solution = solve(problem)
        sol_end = solution[end]
        # output checks
        println("Final Pair Sums Less than P1?: ", sol_end[5] + sol_end[6] < sol_end[1])
        println("Minimum Dynamic Variable is: ", round(minimum(sol_end); digits=4))
        # plot solution
        plt = plot(
            solution, size=(1600, 1000),
            label=LABELS_plotRun, lw=LINEWIDTHS_plotRun, alpha=LINEOPACITY_plotRun
        )
        display(plt)
    else
        error("{plotRun}: Invalid input.")
    end
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

# plots 5 binary persistence maps, 4 for the raw links + 1 for the combined top predator density
function LiaoTypePersistenceMap(model, baseParams, resourceParams, grain, init, timespan)
    # get simGridData
    simGridData = LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan)
    # setup
    grid_length = 0:grain:1
    n = length(grid_length)
    Z_R = reshape(simGridData.PersistenceR, n, n)
    Z_2 = reshape(simGridData.Persistence2, n, n)
    Z_3 = reshape(simGridData.Persistence3, n, n)
    Z_23 = reshape(simGridData.Persistence23, n, n)
    Z_13 = reshape(simGridData.Persistence13, n, n)
    # create heatmaps
    pR = heatmap(grid_length, grid_length, Z_R, title="Resource Persistence (P₁)",
        aspect_ratio=1, c=cgrad(:grays, rev=true, scale=:linear), clims=(0, 1), colorbar=false)
    p2 = heatmap(grid_length, grid_length, Z_2, title="Consumer Persistence (P₂)",
        aspect_ratio=1, c=cgrad(:grays, rev=true, scale=:linear), clims=(0, 1), colorbar=false)
    p3 = heatmap(grid_length, grid_length, Z_3, title="Tertiary Persistence (P₃)",
        aspect_ratio=1, c=cgrad(:grays, rev=true, scale=:linear), clims=(0, 1), colorbar=false)
    p23 = heatmap(grid_length, grid_length, Z_23, title="2-3 Link Persistence (P₍₂₃₎)",
        aspect_ratio=1, c=cgrad(:grays, rev=true, scale=:linear), clims=(0, 1), colorbar=false)
    p13 = heatmap(grid_length, grid_length, Z_13, title="1-3 Link Persistence (P₍₁₃₎)",
        aspect_ratio=1, c=cgrad(:grays, rev=true, scale=:linear), clims=(0, 1), colorbar=false)
    # arrange plots
    plot(pR, p2, p3, p23, p13, size=(1600, 1000))
end

# assigns each point on the LiaoType plot a color. Each colors is a unique species combination
function LiaoTypeSpeciesRichnessMap(model, baseParams, resourceParams, grain, init, timespan)
    # get augmented sim data
    simGridData = LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan)
    # create diversity dataset
    speciesDiversityData = DataFrame(
        Availability=simGridData.Availability,
        Connectivity=simGridData.Connectivity,
        SpeciesDistributionID=assignSpeciesDistributionID(simGridData)
    )
    # plot on a heatmap
    grid_length = 0:grain:1
    n = length(grid_length)
    Z_S = reshape(speciesDiversityData.SpeciesDistributionID, n, n)
    p = heatmap(grid_length, grid_length, Z_S, title="Species Richness Across Landscape Types",
        aspect_ratio=1, c=cgrad(:viridis, scale=:linear), clims=(0, 1), colorbar=false)
    plot(p, size=(1600, 1000))
end

# Plot the "proportional persistence" of all 5 possible coexistance patterns as a function of <param_i>
# y-axis: interpreted as the fraction of all possible habitat types in which this coexistance pattern is observed
# Note: it is assumed that <param_i> is the index of a parameter that ONLY varies between 0 and 1
function ProportionalPersistencePlot(param_i, model, baseParams, resourceParams, grain, init, timespan)
    # check input
    if !(param_i in CANDIADATE_PP_PARAMETER_INDICES)
        error("{ProportionalPersistencePlot}: Parameter index chosen is invalid or not yet accounted for.")
    end
    # generate plotting data
    plottingData = DataFrame(ParameterValue=Float64[], PP_1=Float64[],
        PP_2=Float64[], PP_3=Float64[], PP_4=Float64[], PP_5=Float64[])
    for parameter_value in 0.025:PP_GRAIN:0.5
        # alter model parmeter inputs as required
        baseParams1 = copy(baseParams)
        resourceParams1 = copy(resourceParams)
        if param_i < 12
            baseParams1[param_i] = parameter_value
        elseif param_i >= 14 && param_i <= 16
            resourceParams1[param_i-14+1] = parameter_value
        else
            error("{ProportionalPersistencePlot}: Parameter index chosen is invalid.")
        end
        # run model and push PP results to plotting data
        pp = ProportionalPersistencePoint(model, baseParams1, resourceParams1, grain, init, timespan)
        push!(plottingData, (ParameterValue=parameter_value, PP_1=pp[1], PP_2=pp[2],
            PP_3=pp[3], PP_4=pp[4], PP_5=pp[5]))
    end
    # plot
    colors = ["#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"]
    plt = plot(plottingData.ParameterValue, plottingData.PP_1, label="PP₁",
        xlabel=get_label_from_param_i(param_i),
        ylabel="Proportional Persistence", lw=10, color=colors[1], size=(1600, 1000))
    plot!(plt, plottingData.ParameterValue, plottingData.PP_2, label="PP₂", lw=10, color=colors[2])
    plot!(plt, plottingData.ParameterValue, plottingData.PP_3, label="PP₃", lw=10, color=colors[3])
    plot!(plt, plottingData.ParameterValue, plottingData.PP_4, label="PP₄", lw=10, color=colors[4])
    plot!(plt, plottingData.ParameterValue, plottingData.PP_5, label="PP₅", lw=10, color=colors[5])
end

# Plot the "weighted proportional persistence" of all 5 possible coexistance patterns as a function of <param_i>
# Note: it is assumed that <param_i> is the index of a parameter that ONLY varies between 0 and 1
function WeightedProportionalPersistencePlot(param_i, model, baseParams, resourceParams, grain, init, timespan)
    # check input
    if !(param_i in CANDIADATE_PP_PARAMETER_INDICES)
        error("{ProportionalPersistencePlot}: Parameter index chosen is invalid or not yet accounted for.")
    end
    # generate plotting data
    plottingData = DataFrame(ParameterValue=Float64[], PP_1=Float64[],
        PP_2=Float64[], PP_3=Float64[], PP_4=Float64[], PP_5=Float64[])
    for parameter_value in 0.025:PP_GRAIN:0.5
        # alter model parmeter inputs as required
        baseParams1 = copy(baseParams)
        resourceParams1 = copy(resourceParams)
        if param_i < 12
            baseParams1[param_i] = parameter_value
        elseif param_i >= 14 && param_i <= 16
            resourceParams1[param_i-14+1] = parameter_value
        else
            error("{ProportionalPersistencePlot}: Parameter index chosen is invalid.")
        end
        # run model and push PP results to plotting data
        pp = WeightedProportionalPersistencePoint(model, baseParams1, resourceParams1, grain, init, timespan)
        push!(plottingData, (ParameterValue=parameter_value, PP_1=pp[1], PP_2=pp[2],
            PP_3=pp[3], PP_4=pp[4], PP_5=pp[5]))
    end
    # plot
    colors = ["#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"]
    plt = plot(plottingData.ParameterValue, plottingData.PP_1, label="PP₁",
        xlabel=get_label_from_param_i(param_i),
        ylabel="Proportional Persistence", lw=10, color=colors[1], size=(1600, 1000))
    plot!(plt, plottingData.ParameterValue, plottingData.PP_2, label="PP₂", lw=10, color=colors[2])
    plot!(plt, plottingData.ParameterValue, plottingData.PP_3, label="PP₃", lw=10, color=colors[3])
    plot!(plt, plottingData.ParameterValue, plottingData.PP_4, label="PP₄", lw=10, color=colors[4])
    plot!(plt, plottingData.ParameterValue, plottingData.PP_5, label="PP₅", lw=10, color=colors[5])
end

# TODO: replace the data generation logic of these plots with it's own dedicated function and have a switch: weighted = true or false