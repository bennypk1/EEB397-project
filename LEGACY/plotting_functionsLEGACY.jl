
function LiaoTypeGridLEGACY(model, baseParams, resourceParams, grain, init, timespan)
    # input check
    if !LiaoTypeGridValidInput(baseParams, resourceParams, grain, init, timespan)
        error("{LiaoTypeGrid}: Invalid input.")
    end
    # setup grid and dataframe
    grid = create_unit_grid(grain)
    data = DataFrame(Availability=Float64[], Connectivity=Float64[],
        DensityR=Float64[], Density2=Float64[], Density3=Float64[],
        Density23=Float64[], Density13=Float64[])
    # populate dataframe
    for point in grid
        # if point is invalid, save NaNs
        if !is_gridPoint_valid(point)
            new_row = (point[1], point[2], NaN, NaN, NaN, NaN, NaN)
            push!(data, new_row)
            # otherwise, run ODE and save end points
        else
            UF = transform_goodLandscape_params([point[1], point[2]])
            curr_params = [baseParams; UF; resourceParams]
            curr_prob = ODEProblem(model, init, timespan, curr_params)
            curr_sol = solve(curr_prob, Tsit5(); callback=cb)
            sol_end = curr_sol[end]
            # if the end solution is mathematically invalid, throw error
            if !is_valid_ODE_end(point, sol_end)
                error("{LiaoTypeGrid}: ODE solution violates model variable constraints.")
            end
            new_row = (point[1], point[2], sol_end[1], sol_end[2], sol_end[3] + sol_end[4], sol_end[3], sol_end[4])
            push!(data, new_row)
        end
    end
    return data
end

# Note: only saves data for the 3 higher-order coexistance patterns
function WeightedFragmentationPersistenceDataLEGACY(param_i, model, baseParams, resourceParams, grain, init, timespan)
    # initialize return dataframe
    returnData = DataFrame(ParameterValue=Float64[], Avaliability=Float64[],
        PP_3=Float64[], PP_4=Float64[], PP_5=Float64[])

    # for each fixed parameter value...
    for parameter_value in 0.025:PP_GRAIN:0.5

        # 1) alter model parmeter inputs as required, and get data
        baseParams1 = copy(baseParams)
        resourceParams1 = copy(resourceParams)
        if param_i < 12
            baseParams1[param_i] = parameter_value
        elseif param_i >= 14 && param_i <= 16
            resourceParams1[param_i-14+1] = parameter_value
        else
            error("{WeightedFragmentationPersistenceData}: Parameter index chosen is invalid.")
        end
        simGridData = LiaoTypeGridExtra(model, baseParams1, resourceParams1, grain, init, timespan)
        simGridData.SpeciesDistributionID = assignSpeciesDistributionID(simGridData)

        # 2) for each fixed Availability point, save results
        for availability in 0:grain:1
            # setup
            fragmentation_gradient_data = simGridData[simGridData.Availability.==availability, :]
            n_valid_landscapes = count(!isnan, fragmentation_gradient_data.SpeciesDistributionID)
            pp = zeros(3)
            # for each higher-order coexistance pattern, get it's weighted proportional persistence over the fragmentation gradient
            for i in 3:5
                curr_ID = SPECIES_DISTRIBUTION_IDS[i]
                filtered_data = fragmentation_gradient_data[fragmentation_gradient_data.SpeciesDistributionID.==curr_ID, :]
                # handle if this species coexistance pattern does not exist in the given data
                if nrow(filtered_data) == 0
                    pp[i-2] = 0.0
                else
                    pp[i-2] = sum(filtered_data.DensityR) / n_valid_landscapes
                end
            end
            # save results
            push!(returnData, (ParameterValue=parameter_value, Avaliability=availability,
                PP_3=pp[1], PP_4=pp[2], PP_5=pp[3]))
        end
    end
    return returnData
end

# plots 5 binary persistence maps, 4 for the raw links + 1 for the combined top predator density
function LiaoTypePersistenceMapLEGACY(model, baseParams, resourceParams, grain, init, timespan)
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

# Plot the "proportional persistence" of all 5 possible coexistance patterns as a function of <param_i>
# y-axis: interpreted as the fraction of all possible habitat types in which this coexistance pattern is observed
# Note: it is assumed that <param_i> is the index of a parameter that ONLY varies between 0 and 1
function ProportionalPersistencePlotLEGACY(param_i, model, baseParams, resourceParams, grain, init, timespan)
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
function WeightedProportionalPersistencePlotLEGACY(param_i, model, baseParams, resourceParams, grain, init, timespan)
    # check input and generate plotting data
    if !(param_i in CANDIADATE_PP_PARAMETER_INDICES)
        error("{ProportionalPersistencePlot}: Parameter index chosen is invalid or not yet accounted for.")
    end
    plottingData = WeightedProportionalPersistenceData(param_i, model, baseParams, resourceParams, grain, init, timespan)
    # plot
    colors = COEXISTANCE_PATTERN_COLORS
    plt = plot(plottingData.ParameterValue, plottingData.PP_3, label="Species 2",
        xlabel=get_label_from_param_i(param_i),
        ylabel="Proportional Persistence", lw=10, color=colors[3], size=(1600, 1000))
    plot!(plt, plottingData.ParameterValue, plottingData.PP_4, label="Species 3", lw=10, color=colors[4])
    plot!(plt, plottingData.ParameterValue, plottingData.PP_5, label="Species 2 & 3", lw=10, color=colors[5])
end

# plots the weighted proportional persistence of the 3  coexistance patterns as a function of <param_i>.
# Note: here, persistence is measured for a fixed habitat availability, and therefore measures the proportion of fragmentation extents it can occupy
function WeightedFragmentationPersistencePlotLEGACY(param_i, model, baseParams, resourceParams, grain, init, timespan)
    # check input, then generate plotting data
    if !(param_i in CANDIADATE_PP_PARAMETER_INDICES)
        error("{ProportionalPersistencePlot}: Parameter index chosen is invalid or not yet accounted for.")
    end
    rawData = WeightedFragmentationPersistenceData(param_i, model, baseParams, resourceParams, grain, init, timespan)
    plottingData = flatten_PP_statistics(rawData)

    # plot
    plot(plottingData.ParameterValue, plottingData.PP_3_mean, ribbon=plottingData.PP_3_sd,
        linewidth=2, label="PP_3", xlabel="Parameter Value", ylabel="Mean Persistence (± SD)")
end

# generate a dataframe of weighted proportional persistence
function WeightedProportionalPersistenceDataLEGACY(param_i, model, baseParams, resourceParams, grain, init, timespan)
    # init dataframe
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
            error("{WeightedProportionalPersistenceData}: Parameter index chosen is invalid.")
        end
        # run model and push Weighted PP results to plotting data
        pp = WeightedProportionalPersistencePoint(model, baseParams1, resourceParams1, grain, init, timespan)
        push!(plottingData, (ParameterValue=parameter_value, PP_1=pp[1], PP_2=pp[2],
            PP_3=pp[3], PP_4=pp[4], PP_5=pp[5]))
    end
    return plottingData
end