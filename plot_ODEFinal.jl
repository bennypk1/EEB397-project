using Revise
using DifferentialEquations
using Plots
using DataFrames
using Statistics
include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))
include(joinpath(@__DIR__, "CONSTANTS", "visualConventions.jl"))


# Note: each ODE has 6 variables (none dependent directly on time) and 16 parameters
# Note: of a coexistance pattern, proportional persistance := the area of the Liao Type plot it occupies

# add an early-stopping threshold when P₁ gets very small, to avoid numerical instability.
P1_STOPPING_THRESHOLD = 1e-6
condition(u, t, integrator) = u[1] - P1_STOPPING_THRESHOLD
affect!(integrator) = begin
    terminate!(integrator)
end
cb = ContinuousCallback(condition, affect!)

################################################################################################################
# DATA GENERATING FUNCTIONS
################################################################################################################

# simulate model accross grid and return data (optimized data allocation comapred to LEGACY)
# NOTE: <baseParams> consists of the first 11 nuisance parameters
# NOTE: <resourceParams> are the last 3 parameters that specify the resource life-history
function LiaoTypeGrid(model, baseParams, resourceParams, grain, init, timespan)
    # input check
    if !LiaoTypeGridValidInput(baseParams, resourceParams, grain, init, timespan)
        error("{LiaoTypeGrid}: Invalid input.")
    end
    # create grid, empty dataframe ODE problem
    grid = create_unit_grid(grain)
    gridSize = length(grid)
    data = DataFrame( # fastest way to initialize this kind of dataframe
        Availability=Vector{Float64}(undef, gridSize),
        Connectivity=Vector{Float64}(undef, gridSize),
        DensityR=Vector{Float64}(undef, gridSize),
        Density2=Vector{Float64}(undef, gridSize),
        Density3=Vector{Float64}(undef, gridSize),
        Density23=Vector{Float64}(undef, gridSize),
        Density13=Vector{Float64}(undef, gridSize)
    )
    # preallocate problem and parameter vector
    p = [baseParams; zeros(2); resourceParams]
    p_useful_index = length(baseParams) + 1
    base_problem = ODEProblem(model, init, timespan, p)
    # populate dataframe
    for (i, point) in enumerate(grid)
        # save the current point
        data.Availability[i] = point[1]
        data.Connectivity[i] = point[2]
        # point is invalid, fill with NaN
        if !is_gridPoint_valid(point)
            data.DensityR[i] = NaN
            data.Density2[i] = NaN
            data.Density3[i] = NaN
            data.Density23[i] = NaN
            data.Density13[i] = NaN
            # otherwise, cacluate ODE solution
        else
            # create current parameter vector
            p[p_useful_index:p_useful_index+1] .= transform_goodLandscape_params(point)
            # initializeand solve ODE problem
            curr_problem = remake(base_problem; p=p) # using remake for performance
            curr_solution = solve(curr_problem, Rodas5(); callback=cb) # using Tsit5() for performance
            sol_end = curr_solution[end]
            # if the end solution is mathematically invalid, throw error
            if !is_valid_ODE_end(point, sol_end)
                println(sol_end)
                error("{LiaoTypeGrid}: ODE solution violates model variable constraints.")
            end
            # update dataframe
            data.DensityR[i] = sol_end[1]
            data.Density2[i] = sol_end[2]
            data.Density3[i] = sol_end[3] + sol_end[4]
            data.Density23[i] = sol_end[3]
            data.Density13[i] = sol_end[4]
        end
    end
    return (data)
end

# augment LiaoTypeGrid dataframe with persistence data
function LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan; threshold=0.05)
    # input check
    if length(baseParams) != 11
        error("{ProportionalPersistencePoint}: inputted baseParams has a length different than 11.")
    end
    if length(resourceParams) != 3
        error("{ProportionalPersistencePoint}: inputted resourceParams has a length different than 3.")
    end
    # get raw grid data simGridData
    rawGridData = LiaoTypeGrid(model, baseParams, resourceParams, grain, init, timespan)
    # return a dataframe augmented with binary persistence data
    newData = DataFrame(
        Availability=rawGridData.Availability, Connectivity=rawGridData.Connectivity,
        DensityR=rawGridData.DensityR, Density2=rawGridData.Density2, Density3=rawGridData.Density3,
        Density23=rawGridData.Density23, Density13=rawGridData.Density13,
        PersistenceR=discretizeDensity(rawGridData.DensityR, threshold),
        Persistence2=discretizeDensity(rawGridData.Density2, threshold),
        Persistence3=discretizeDensity(rawGridData.Density3, threshold),
        Persistence23=discretizeDensity(rawGridData.Density23, threshold),
        Persistence13=discretizeDensity(rawGridData.Density13, threshold))
    return newData
end

# generates a single point vector, to be plotted in <ProportionalPersistencePlot>
function ProportionalPersistencePoint(model, baseParams, resourceParams, grain, init, timespan; threshold=0.05)
    # get augmented sim data, add distribution IDs, and save it's number of rows
    simGridData = LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan; threshold)
    simGridData.SpeciesDistributionID = assignSpeciesDistributionID(simGridData)
    n_valid_landscapes = count(!isnan, simGridData.SpeciesDistributionID)
    # input check 2
    if !(n_valid_landscapes > 0)
        error("{ProportionalPersistencePoint}: simGridData somehow has no valid landscapes.")
    end
    # for each unique distribution ID, caclulate its proportional persistence
    pointVector = zeros(length(SPECIES_DISTRIBUTION_IDS))
    for i in eachindex(SPECIES_DISTRIBUTION_IDS)
        curr_ID = SPECIES_DISTRIBUTION_IDS[i]
        pointVector[i] = count(simGridData.SpeciesDistributionID .== curr_ID) / n_valid_landscapes
    end
    # output check
    if abs(sum(pointVector) - 1) > grain^2 # checking threshold chosen arbitrarily
        error("{ProportionalPersistencePoint}: Point vector sum is very different from 1")
    end
    return pointVector
end

# generates a single point vector, to be plotted in <WeightedProportionalPersistencePlot>
# Note: this function weighs the proportion of each coexistance pattern by the mean equilibrium resource density accross it's range
function WeightedProportionalPersistencePoint(model, baseParams, resourceParams, grain, init, timespan; threshold=0.05)
    # get augmented sim data, add distribution IDs, and save it's number of rows
    simGridData = LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan; hreshold)
    simGridData.SpeciesDistributionID = assignSpeciesDistributionID(simGridData)
    n_valid_landscapes = count(!isnan, simGridData.SpeciesDistributionID)
    # input check
    if !(n_valid_landscapes > 0)
        error("{WeightedProportionalPersistencePoint}: simGridData somehow has no valid landscapes.")
    end
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

# generalized function to generate proportional persistence data accross a customizable parameter gradient
function ProportionalPersistenceData(param_i, param_lb, param_ub, model, baseParams, resourceParams, landscapeGrain, init, timespan; threshold=0.05, pointFunction=ProportionalPersistencePoint, PP_grain=0.1)
    # check input
    if param_lb < 0.025 || param_ub > 1
        error("{ProportionalPersistenceData}: Input parameter bounds are invalid.")
    end
    # init dataframe
    plottingData = DataFrame(ParameterValue=Float64[], PP_1=Float64[],
        PP_2=Float64[], PP_3=Float64[], PP_4=Float64[], PP_5=Float64[])

    for parameter_value in param_lb:PP_grain:param_ub
        # alter model parmeter inputs as required
        baseParams1 = copy(baseParams)
        resourceParams1 = copy(resourceParams)
        if param_i < 12
            baseParams1[param_i] = parameter_value
        elseif param_i >= 14 && param_i <= 16
            resourceParams1[param_i-14+1] = parameter_value
        else
            error("{ProportionalPersistenceData}: Parameter index chosen is invalid.")
        end
        # run model and push PP results to plotting data
        pp = pointFunction(model, baseParams1, resourceParams1, landscapeGrain, init, timespan; threshold)
        push!(plottingData, (ParameterValue=parameter_value, PP_1=pp[1], PP_2=pp[2],
            PP_3=pp[3], PP_4=pp[4], PP_5=pp[5]))
    end
    return plottingData
end

################################################################################################################
# GENUINE PLOTTING FUNCTIONS
################################################################################################################

# plots a single run of a fully specified simulation
function plotRun(model, params, init, timespan)
    # input check
    if plotRunValidInput(params, init, timespan)
        problem = ODEProblem(model, init, timespan, params)
        solution = solve(problem; callback=cb)
        sol_end = solution[end]
        # output checks
        # println("Final Pair Sums Less than P1?: ", sol_end[5] + sol_end[6] < sol_end[1])
        # println("Minimum Dynamic Variable is: ", round(minimum(sol_end); digits=4))
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
function LiaoTypeHeatMap(model, baseParams, resourceParams, grain, init, timespan; focalValue=-1.0)
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
    if focalValue == -1.0
        # create heatmaps
        pR = heatmap(grid_length, grid_length, Z_R, title="Resource Density (P₁)",
            aspect_ratio=1, c=cgrad(:jet, scale=:linear))
        p2 = heatmap(grid_length, grid_length, Z_2, title="Consumer Density (P₂)",
            aspect_ratio=1, c=cgrad(:jet, scale=:linear))
        p3 = heatmap(grid_length, grid_length, Z_3, title="Tertiary Density (P₃)",
            aspect_ratio=1, c=cgrad(:jet, scale=:linear))
        p23 = heatmap(grid_length, grid_length, Z_23, title="2-3 Link Density (P₍₂₃₎)",
            aspect_ratio=1, c=cgrad(:jet, scale=:linear))
        p13 = heatmap(grid_length, grid_length, Z_13, title="1-3 Link Density (P₍₁₃₎)",
            aspect_ratio=1, c=cgrad(:jet, scale=:linear))
        # arrange plots
        plot(pR, p2, p3, p23, p13, size=(1600, 1000))
    elseif focalValue == 2.0
        # plot just species 2
        plot(heatmap(grid_length, grid_length, Z_2,
            aspect_ratio=1, c=cgrad(:jet, scale=:linear), size=(1600, 1000), clims=(0, 1))) # temportary clims here
    elseif focalValue == 1.0
        # plot just resource species
        plot(heatmap(grid_length, grid_length, Z_R,
            aspect_ratio=1, c=cgrad(:jet, scale=:linear), size=(1600, 1000), clims=(0, 1)))
    else
        error("Haven't gotten there yet...")
    end
end

# assigns each point on the LiaoType plot a color. Each colors is a unique species combination
function LiaoTypeSpeciesRichnessMap(model, baseParams, resourceParams, grain, init, timespan; threshold=0.05)
    # get augmented sim data
    simGridData = LiaoTypeGridExtra(model, baseParams, resourceParams, grain, init, timespan; threshold)
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
    p = heatmap(grid_length, grid_length, Z_S, aspect_ratio=1, c=cgrad(:viridis, scale=:linear), clims=(0, 1), colorbar=false)
    return p
end

# plot species persistence as a function of a single changing parameter
function SpeciesPersistencePlot(param_i, param_lb, param_ub, model, baseParams, resourceParams; landscapeGrain=0.025, init=CANONICAL_INIT, timespan=CANONICAL_TIMESPAN)
    # first, get data
    ppData = ProportionalPersistenceData(param_i, param_lb, param_ub, model, baseParams, resourceParams,
        landscapeGrain, init, timespan)
    # plot data
    plot1 = plot(ppData.ParameterValue, 1 .- ppData.PP_1) # resource plot
    plot!(ppData.ParameterValue, ppData.PP_3 + ppData.PP_5) # species 2
    plot!(ppData.ParameterValue, ppData.PP_4 + ppData.PP_5) # species 3
    return plot1
end

# plot coexistance pattern persistence as a function of a single changing parameter
function PatternPersistencePlot(param_i, param_lb, param_ub, model, baseParams, resourceParams; landscapeGrain=0.025, init=CANONICAL_INIT, timespan=CANONICAL_TIMESPAN)
    # first, get data
    ppData = ProportionalPersistenceData(param_i, param_lb, param_ub, model, baseParams, resourceParams,
        landscapeGrain, init, timespan)
    # Viridis hex codes
    custom_colors = ["#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"]
    # plot data
    plot1 = plot(ppData.ParameterValue, ppData.PP_2, lw=5, markershape=:circle, ms=12, color=custom_colors[2]) # resource only
    plot!(ppData.ParameterValue, ppData.PP_3, lw=5, markershape=:circle, ms=12, color=custom_colors[3]) # resource + species 2
    plot!(ppData.ParameterValue, ppData.PP_4, lw=5, markershape=:circle, ms=12, color=custom_colors[4]) # resource + species 3
    plot!(ppData.ParameterValue, ppData.PP_5, lw=5, markershape=:circle, ms=12, color=custom_colors[5]) # all species
    return plot1
end

###############################################################################################################
# PLOTS SPECIFICALLY FOR RESULTS SECTION
###############################################################################################################

# plot-prettying function for <LiaoTypeSpeciesRichnessMap>
function FigureGradeLiaoTypePlot(raw_LiaoTypePlot)
    plot!(
        raw_LiaoTypePlot,
        xlims=(0, 1),
        guidefontsize=20,
        tickfontsize=16,
        titlefontsize=22,
        legendfontsize=16,
        grid=false,
        minorgrid=false,
        size=(1600, 1000),
        framestyle=:box
    )
    return raw_LiaoTypePlot
end

function Figure3_2b(ppData_gamma, ppData_e1)
    # input check
    if nrow(ppData_gamma) != nrow(ppData_e1)
        error("{Figure3_2b}: dataframe column lengths differ.")
    end
    # base plot
    plot1 = plot(ppData_gamma.ParameterValue,
        1 .- ppData_gamma.PP_1, # resource PP for gamma
        seriestype=:path,
        marker=:utriangle,
        xlims=(0, max(maximum(ppData_gamma.ParameterValue), maximum(ppData_e1.ParameterValue)) + 0.1),
        ylims=(0, 1),
        lw=5,
        ms=12,
        guidefontsize=20,
        tickfontsize=16,
        size=(2000, 1000),
        grid=false,
        minorgrid=false,
        framestyle=:box,
        legend=false,
        color="#154360" # dark blue
    )
    # add other lines
    plot!(ppData_gamma.ParameterValue,
        ppData_gamma.PP_3 + ppData_gamma.PP_4 + ppData_gamma.PP_5, marker=:utriangle, lw=5, ms=12, color="#7B241C")
    plot!(ppData_e1.ParameterValue, 1 .- ppData_e1.PP_1, marker=:circle, lw=5, ms=12, color="#85C1E9") # light blue
    plot!(ppData_e1.ParameterValue,
        ppData_e1.PP_3 + ppData_e1.PP_4 + ppData_e1.PP_5, marker=:circle, lw=5, ms=12, color="#F1948A")
    return plot1

    # what the colors mean:
    # blue is resource persistence: light blue is for e1 and dark for gamma
    # red is consumer persistence: light red is for e1 and dar for gamma
    # traingle for gamma, circle for e
end

function prettyPatternPlot!(p)
    plot!(p,
        size=(2000, 1000),
        grid=false,
        minorgrid=false,
        guidefontsize=20,
        tickfontsize=16,
        titlefontsize=22,
        legend=false,
        framestyle=:box
    )
    return p
end