
using DataFrames
include(joinpath(@__DIR__, "sensitivity_constants.jl"))

# This file is a set of helper functions necessary to run the sensitivity (/robustness) analysis of qualitative
# model results.

###############################################################################################################
# CHECKING PARAMETER VALIDITY
###############################################################################################################

# <candidateNuisanceParms> must be of the form: (cₙ, eₙ, μₙ, α)
function is_candidate_nuisance_marginally_valid(candidateNuisanceParms)
    # check input length
    if length(candidateNuisanceParms) != 4
        error("{is_candidate_nuisance_marginally_valid}: input vector does not have length 4.")
    end
    # unpack input vector
    cₙ, eₙ, μₙ, α = candidateNuisanceParms
    # return validity check
    return (
        cₙ > 0.0 &&
        eₙ > 0.0 &&
        cₙ <= SOFT_Cn_MAX &&
        0.0 <= α <= 1.0 &&
        μₙ >= 0.0 &&
        cₙ > eₙ + μₙ &&
        μₙ < eₙ
    )
end

# <mainParms> must be of the form: (e₁, γ, δ)
function is_main_marginally_valid(mainParms)
    # check input length
    if length(mainParms) != 3
        error("{is_main_marginally_valid}: input vector does not have length 3.")
    end
    # unpack input vector
    e₁, γ, δ = mainParms
    # return validity check
    return (e₁ > 0.0 &&
            e₁ >= SOFT_E1_MIN &&
            e₁ < 1.0 &&
            e₁ <= SOFT_E1_MAX &&
            δ > 0.0 &&
            δ >= SOFT_DELTA_MIN &&
            δ <= SOFT_DELTA_MAX &&
            γ >= 0.0 &&
            γ <= SOFT_GAMMA_MAX
    )
end

function is_parms_valid_point(candidateNuisanceParms, mainParms)
    # check marginals
    if !(is_main_marginally_valid(mainParms) && is_candidate_nuisance_marginally_valid(candidateNuisanceParms))
        error("{is_parms_valid_point}: At least one marginal condition of either main or nuisanced parms failed.")
    end
    # check joint constraints
    cₙ, eₙ, μₙ, α = candidateNuisanceParms
    e₁, γ, δ = mainParms
    # TODO: add joint constraints if necessary
    return true
end

function generate_mainParam_grid(; search_grain=SEARCH_GRAIN)
    # define values in valid ranges
    e1_vals = range(SOFT_E1_MIN, stop=SOFT_E1_MAX, length=search_grain)
    gamma_vals = range(0.0, stop=SOFT_GAMMA_MAX, length=search_grain)
    delta_vals = range(SOFT_DELTA_MIN, stop=SOFT_DELTA_MAX, length=search_grain)
    # setup return df
    returnData = DataFrame(e1=Float64[], gamma=Float64[], delta=Float64[])
    # populate and return df
    for e1 in e1_vals, γ in gamma_vals, δ in delta_vals
        push!(returnData, (e1, γ, δ))
    end
    return returnData
end


# this function should test the input vector against a dataframe of main params
function is_parms_valid(candidateNuisanceParms; noJoint=true)
    # check marginal validity of nuisance parameters
    if !is_candidate_nuisance_marginally_valid(candidateNuisanceParms)
        error("{is_parms_valid}: Candidate nuisance parameters failed marginal checks.")
    end
    # if no Joint conditions have been defined yet, return true immediately
    if noJoint
        return true
    end
    # generate main parameter grid and check validity
    mainParamGrid = generate_mainParam_grid()
    for row in eachrow(mainParamGrid)
        if !is_parms_valid_point(candidateNuisanceParms, [row.e1, row.gamma, row.delta])
            return false
        end
    end
    # if all points passed return true
    return true
end

###############################################################################################################
# AUGMENTING PROCESSED (AND VALIDATED) PARAMETER DATA
###############################################################################################################

# map <compactParms> to full parameter list (excluding U, F)
function assemble_parms(compactParms)
    # unpack parameter vectors
    cₙ, eₙ, μₙ, α, e₁, γ = compactParms
    # assign (in order)
    c₂ = cₙ
    c₃₂ = cₙ
    c₃₁ = cₙ
    e₂ = eₙ
    e₃₂ = eₙ
    e₃₁ = eₙ
    μ₂₁ = μₙ
    μ₃₂ = μₙ
    μ₃₁ = μₙ
    z = 4
    β = 1 - α
    return ((c₂, c₃₂, c₃₁, e₁, e₂, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, α, β, γ))
end

# generate new dataframe of full parameters from reduce list of sampled parametrs
function augment_compact_parms(compactParms)
    # input check
    if !(ncol(compactParms) == 6)
        error("{augment_compact_parms}: input dataframe does not have exactly 6 columns.")
    end
    # initialize full dataframe
    fullData = DataFrame(c2=Float64[], c32=Float64[], c31=Float64[],
        e1=Float64[], e2=Float64[], e32=Float64[], e31=Float64[],
        mu21=Float64[], mu32=Float64[], mu31=Float64[],
        z=Float64[], alpha=Float64[], beta=Float64[], gamma=Float64[])
    # iterate over each row of nuisance parameters
    for row in eachrow(compactParms)
        curr_out = assemble_parms(row)
        new_fullData_row = (c2=curr_out[1], c32=curr_out[2], c31=curr_out[3],
            e1=curr_out[4], e2=curr_out[5], e32=curr_out[6], e31=curr_out[7],
            mu21=curr_out[8], mu32=curr_out[9], mu31=curr_out[10],
            z=curr_out[11], alpha=curr_out[12], beta=curr_out[13], gamma=curr_out[14])
        push!(fullData, new_fullData_row)
    end
    return fullData
end