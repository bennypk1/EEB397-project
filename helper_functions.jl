using Revise
using DifferentialEquations
using DataFrames
include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "CONSTANTS", "visualConventions.jl"))

################################################################################################################
# MAIN HELPERS
################################################################################################################

# Create a unit grid with a set grain
function create_unit_grid(grain)
    result = []
    s_possibilities = 0:grain:1
    q_possibilities = 0:grain:1
    for s_point in s_possibilities
        for q_point in q_possibilities
            push!(result, [s_point, q_point])
        end
    end
    return result
end

# Check if a grid point [S, qₛₛ] is valid
function is_gridPoint_valid(grid_point)
    c1 = grid_point[2] > (2 - 1 / grid_point[1])
    c2 = grid_point[1] > 0.025
    c3 = grid_point[2] > 0.025
    c4 = grid_point[1] < 1
    c5 = grid_point[2] < 1
    return c1 && c2 && c3 && c4 && c5
end

function plotRunValidInput(all_params, all_init, timespan)
    c1 = is_params_valid(all_params)
    c2 = is_proportion(all_init)
    c3 = is_timespan_valid(timespan)
    c4 = length(all_init) == 6
    return c1 && c2 && c3 && c4
end

function LiaoTypeGridValidInput(baseParams, resourceParams, grain, init, timespan)
    c0 = all([resourceParams; baseParams] .>= 0)
    c1 = is_proportion(init)
    c2 = is_timespan_valid(timespan)
    c3 = grain < 0 || grain > 1
    return c1 && c2 && !c3 && c0
end

function mean_nonzero_diff(data; tol=1e-12)
    d = diff(data)
    # Filter for absolute values larger than tolerance
    nonzero_d = d[abs.(d).>tol]
    return mean(nonzero_d)
end

function no_meaningful_slope_violation(resourceDat, consumerDat; tol=0.01)
    # input checks
    if length(resourceDat) != length(consumerDat)
        error("{meaningful_slope_violation}: resource and consumer threshold Vectors have differing lengths.")
    end
    if !(all(diff(resourceDat) .<= 0.0) && all(diff(consumerDat) .<= 0.0))
        error("{TargetCommunity_consumerper...}: some estimated slopes are positive.")
    end
    # if consumer difference from [1] to [2] if very small, disregard this example
    if (consumerDat[2] - consumerDat[1] < tol)
        return true
    end
    # get slopes
    slopeR = mean_nonzero_diff(resourceDat)
    slopeC = mean_nonzero_diff(consumerDat)

    # verify claim
    if (abs(slopeC) - abs(slopeR) > -tol)
        return true
    else
        println([slopeR, slopeC])
        return false
    end
end

# check if the output of solve(ODEProblem) is valid for our model, in a given a landscape defined by <point> (non-exhautive)
function is_valid_ODE_end(point, sol_end)
    if !all(sol_end .>= -1e-4)
        println("Condition 1 failed: At least one final density is significantly negative.")
        return false
    end
    if !all(sol_end .<= point[1])
        println("Condition 2 failed: At least one final density is greater than landscape suitability.")
        return false
    end
    if !(sol_end[2] <= sol_end[1] || sol_end[2] <= SIGNIFICANT_P)
        println("Condition 3 failed: P₂ is significant and larger than P₁.")
        return false
    end
    if !(sol_end[3] <= sol_end[1] || sol_end[3] <= SIGNIFICANT_P)
        println("Condition 4 failed: P₍₁₃₎ is significant and larger than P₁.")
        return false
    end
    if !(sol_end[4] <= sol_end[1] || sol_end[4] <= SIGNIFICANT_P)
        println("Condition 5 failed: P₍₂₃₎ is significant and larger than P₁.")
        return false
    end
    if !(sol_end[3] + sol_end[4] <= sol_end[1] || any(sol_end[3:4] .<= SIGNIFICANT_P))
        println("Condition 3 failed: P₃ is significant and larger than P₁.")
        return false
    end
    return true
end

# returns vector mapping NaN to NaN, and proportions to 0 or 1 depending on <threshold>
# <LiaoTypeVector> must be a vector with only NaN, and Floats from 0-1
function discretizeDensity(LiaoTypeVector, threshold)
    # check that threshold is valid
    if threshold > 1 || threshold < 0
        error("{discretizeDensity} Threshold too large or invalid.")
    end
    # real function body
    returnVector = similar(LiaoTypeVector, Float64)
    for i in eachindex(LiaoTypeVector)
        element = LiaoTypeVector[i]
        if isnan(element)
            returnVector[i] = NaN
        elseif element <= 1 && element >= threshold
            returnVector[i] = 1.0
        elseif element >= -1e-4 && element < threshold
            returnVector[i] = 0.4
        else
            error("{discretizeDensity} Wrong elelment type.")
        end
    end
    return returnVector
end

# returns a vector of discrete fractions, each corresponding to a unique combination of species
# <abundanceGridData> must be output from LiaoTypeGridExtra
function assignSpeciesDistributionID(abundanceGridData)
    n = nrow(abundanceGridData)
    returnVector = zeros(n)
    for i in 1:nrow(abundanceGridData)
        raw_id = [
            abundanceGridData.PersistenceR[i],
            abundanceGridData.Persistence2[i],
            abundanceGridData.Persistence3[i]]
        returnVector[i] = mapRawIDtoSpeciesDistributionID(raw_id)
    end
    return returnVector
end

# flatten a dataframe output from <>
function flatten_PP_statistics(df)
    gdf = groupby(df, :ParameterValue)
    flattened = combine(gdf,
        :PP_3 => mean => :PP_3_mean,
        :PP_4 => mean => :PP_4_mean,
        :PP_5 => mean => :PP_5_mean,
        :PP_3 => std => :PP_3_sd,
        :PP_4 => std => :PP_4_sd,
        :PP_5 => std => :PP_5_sd)
    return flattened
end


################################################################################################################
# SMALL HELPERS
################################################################################################################

# Convert [S, qₛₛ] to [U, F]
function transform_goodLandscape_params(goodLandscape_params)
    U = 1 - goodLandscape_params[1]
    F = 1 - goodLandscape_params[2]
    return ([U, F])
end

# Convert [U, F] to [S, qₛₛ]
function transform_badLandscape_params(badLandscape_params)
    S = 1 - badLandscape_params[1]
    qₛₛ = 1 - badLandscape_params[2]
    return ([S, qₛₛ])
end

#check validity of landscape parameters [U, F]
function is_landscapeUF_valid(UF)
    return (is_gridPoint_valid(transform_badLandscape_params(UF)))
end

# check validity of parameters
function is_params_valid(all_params)
    c0 = length(all_params) == 16
    c1 = is_landscapeUF_valid(all_params[12:13])
    c2 = isinteger(all_params[11])
    return c0 && c1 && c2
end

# check if all elevents of a vector are proportions
function is_proportion(v)
    for i in eachindex(v)
        if v[i] < 0 || v[i] > 1
            return false
        end
    end
    return true
end

# check validity of timespan
function is_timespan_valid(the_timespan)
    c0 = length(the_timespan) == 2
    c1 = the_timespan[1] == 0
    c2 = the_timespan[2] > the_timespan[1]
    return c0 && c1 && c2
end

function mapRawIDtoSpeciesDistributionID(rawID)
    if length(rawID) != 3
        error("{mapRawIDtoSpeciesDistributionID} rawID vector not correct length.")
    end
    no = PERSISTENCE_CODE[1]
    yes = PERSISTENCE_CODE[2]
    if isnan(rawID[1])
        return NaN
    elseif rawID == [no, no, no]
        return SPECIES_DISTRIBUTION_IDS[1]
    elseif rawID == [yes, no, no]
        return SPECIES_DISTRIBUTION_IDS[2]
    elseif rawID == [yes, yes, no]
        return SPECIES_DISTRIBUTION_IDS[3]
    elseif rawID == [yes, no, yes]
        return SPECIES_DISTRIBUTION_IDS[4]
    elseif rawID == [yes, yes, yes]
        return SPECIES_DISTRIBUTION_IDS[5]
    else
        error("{mapRawIDtoSpeciesDistributionID} rawID vector is an invalid ID")
    end
end

function get_label_from_param_i(i)
    if i == 4
        return "Intrinsic Extinction Coefficient (e₁)"
    elseif i == 8
        return "Top-Down Predation Coefficient (μ₂₁)"
    elseif i == 16
        return "Crowding Sensitivity Coefficient (γ)"
    else
        error("{get_label_from_param_i}: Parameter index is invalid.")
    end
end

################################################################################################################
# LEGACY HELPERS
################################################################################################################

# convert [S, qₛₛ] to [pᵤ, pᵤᵤ]
function LEGACYconvert_landscape_params(landscape_params)
    S = 1 - landscape_params[1]
    F = this_pᵤ - landscape_params[1] * (1 - landscape_params[2])
    return ([S, F])
end

# performance function for plot_ODE5
function get_init_from_trophic_config(initial_densities, trophic_configuration)
    P1, P12, P23, P13, P11, Pu1 = initial_densities
    I23, I13 = trophic_configuration
    P23_init = I23 == 1 ? P23 : 0.0
    P13_init = I13 == 1 ? P13 : 0.0
    return [P1, P12, P23_init, P13_init, P11, Pu1]
end

################################################################################################################
# GRID DIAGNOSTICS
################################################################################################################

function LiaoTypeGridOutputDiagnostics(model, baseParams, resourceParams, grain, init, timespan)
    # input check
    if !LiaoTypeGridValidInput(baseParams, resourceParams, grain, init, timespan)
        error("{LiaoTypeGrid}: Invalid input.")
    end
    # create grid and empty dataframe
    grid = create_unit_grid(grain)
    data = DataFrame(
        Availability=Float64[], Connectivity=Float64[],
        DensityR=Float64[], Density2=Float64[], Density3=Float64[],
        Density23=Float64[], Density13=Float64[],
        PairDensity11=Float64[], PairDensityU1=Float64[], PairDensitySum=Float64[]
    )
    # populate dataframe
    for point in grid
        if !is_gridPoint_valid(point)
            push!(data, (point[1], point[2], NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN))
        else
            UF = transform_goodLandscape_params([point[1], point[2]])
            curr_params = [baseParams; UF; resourceParams]
            curr_problem = ODEProblem(model, init, timespan, curr_params)
            curr_solution = solve(curr_problem)
            sol_end = curr_solution[end]
            # check HERE
            if !is_valid_ODE_end(point, sol_end)
                error("{LiaoTypeGrid}: ODE solution violates model variable constraints.")
            end
            push!(data, (point[1], point[2],
                sol_end[1], sol_end[2], sol_end[3] + sol_end[4], sol_end[3], sol_end[4],
                sol_end[5], sol_end[6], sol_end[5] + sol_end[6]))
        end
    end
    return (data)
end