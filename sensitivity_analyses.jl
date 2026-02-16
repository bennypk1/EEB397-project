using Revise
using Distributions
include(joinpath(@__DIR__, "helper_functions.jl"))

###############################################################################################################
# DEFINING ALL HARD MODEL CONSTRAINTS (delete from here when implemented in code below)
###############################################################################################################

# z = 4

# TODO: e₁ + μ₂₁ + μ₃₁ < 1

# 2) Dependencies
# c₂ = c₃₂ = c₃₁ = cₙ
# e₃₂ = e₂ = eₙ
# μ₃₁ = μ₃₂ = μₙ
# α + β = 1
# ψ = e₃₁ / e₃₂
# ω = μ₃₁ / μ₂₁

###############################################################################################################
# Under these constraints, the model effectively has 8 free variables (+ threshold constraints ; 7 if the model is not Omnivory)
###############################################################################################################

# Base Parameters (these may be sampled from a parameter space, adhering to hard AND soft constraints)
# 1) cₙ: higher-order colonization rate
# 2) eₙ: higher-order extinction rate
# 3) μₙ: higher-order top-down predation rate
# 4) α: locality of resource dispersal
# 5) ψ: species feeding preference cost (varies in OMNIVORY ONLY)

# Parameters of interest
# 6) e₁: intrinsic extinction rate
# 7) ω:  variation in top-down extinction rate
# 8) γ:  local crowding sensitivity

# Major Structural Assumptions 
# 1) regional dispersal of species 2
# 2) competitive dominance of species 2
# 3) global dispersal of species 3

###############################################################################################################
# DEFINING ADDITIONAL SEARCH (soft) CONSTRAINTS
###############################################################################################################

PSI_UPPER_BOUND_SOFT = 3.0
RESOURCE_INIT_DENSITY_CAP_SOFT = 0.015
PAIR_DENSITY_SUM_CAP_SOFT = 1.0 # OFF
GAMMA_MIN, GAMMA_MAX = [0, 1] # OFF
OMEGA_MIN, OMEGA_MAX = [0.1, 1] # TODO: this is quite conservative
E1_MIN, E1_MAX = [0.05, 0.5] # TODO: figure this out, currently unstable
CANONICAL_MAIN_BASE = [0.05, 1.0, 0.0]

###############################################################################################################
# HELPERS - SAMPLING FUNCTIONS
###############################################################################################################

# Notation
# fixedParameters : 5-element vector
# mainParameters  : 3-element vector
# init            : 6-element vector

# check if fixed parameters adhere to (marginal) hard AND soft constraints
function is_fixed_valid(fixedParameters)
    if length(fixedParameters) != 5
        error("{is_fixed_valid}: fixedParameters input vector does not have length 5.")
    end
    cₙ, eₙ, μₙ, α, ψ = fixedParameters
    return (
        cₙ > 0.0 &&
        eₙ > 0.0 &&
        0.0 <= α <= 1.0 &&
        μₙ >= 0.0 &&
        ψ >= 1.0 &&
        cₙ > eₙ + μₙ &&
        μₙ < eₙ &&
        μₙ < 1.0
    )
end

# check if main parameters adhere to (marginal) hard AND soft constraints
function is_main_valid(mainParameters)
    if length(mainParameters) != 3
        error("{is_main_valid}: mainParameters input vector does not have length 3.")
    end
    e₁, ω, γ = mainParameters
    return (e₁ > 0.0 &&
            ω > 0.0 &&
            γ >= 0.0 &&
            ω <= 1)
end

# check if initial dynamic variables are valid
function is_init_valid(init)
    if length(init) != 6
        error("{is_init_valid}: init input vector does not have length 6.")
    end
    P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = init
    return (P₁ >= 0.0 && P₍₁₂₎ >= 0.0 && P₍₂₃₎ >= 0.0 && P₍₁₃₎ >= 0.0 && P₁₁ >= 0.0 && Pᵤ₁ >= 0.0 &&
            P₁ < RESOURCE_INIT_DENSITY_CAP_SOFT &&
            P₍₁₂₎ < P₁ &&
            P₍₂₃₎ < P₁ &&
            P₍₂₃₎ + P₍₁₃₎ < P₁ &&
            P₁₁ + Pᵤ₁ < P₁ * PAIR_DENSITY_SUM_CAP_SOFT)
end

# check if all parameters satisfy hard AND soft constraints
function is_model_valid(fixedParameters, mainParameters, omnivoryBool)
    # input checks
    marg1 = is_fixed_valid(fixedParameters)
    marg2 = is_main_valid(mainParameters)
    # function body
    cₙ, eₙ, μₙ, α, ψ = fixedParameters
    e₁, ω, γ = mainParameters
    # next, quickly check omnivory
    omniPsi = omnivoryBool || (ψ == 1.0) # if Omnivory, ψ ≥ 1, else, ψ = 1
    # return validity
    return marg1 && marg2 && omniPsi
end

# HELPER
# returns a 3-column dataframe where each row is a potential mainParameters vector. Used EXCLUSIVELY for checking input parameter validity.
# <length.out> determines how finely we are verifying parameter validity
function generate_mainParameters_data(length_out, e₁_base, ω_base, γ_base)
    # check input validity
    if !(grain < 0.25 && grain > 0.0)
        error("{generate_mainParameters_data}: grain is either too coarse (> 0.25), or invalid.")
    end
    if !all([e₁_base, ω_base, γ_base] .== CANONICAL_MAIN_BASE)
        error("{generate_mainParameters_data}: at least one chosen base main parameter is invalid.")
    end
    # setup
    returnData = DataFrame(e_1=Float64[], omega=Float64[], gamma=Float64[])
    # add each parameter variation list individually
    for e₁_curr in range(E1_MIN, stop=E1_MAX, length=length_out)
        push!(returnData, (e₁_curr, ω_base, γ_base))
    end
    for ω_curr in range(OMEGA_MIN, stop=OMEGA_MAX, length=length_out)
        push!(returnData, (e₁_base, ω_curr, γ_base))
    end
    for γ_curr in range(GAMMA_MIN, stop=GAMMA_MAX, length=length_out)
        push!(returnData, (e₁_base, ω_base, γ_curr))
    end
    # check output marginal validity
    for row in eachrow(returnData)
        if !is_main_valid(row)
            error("{generate_mainParameters_data}: runs but is creating marginally invalid main parameter combinations.")
        end
    end
    return returnData
end

# HELPER
function SAMPLERAW_fixedParameters(omnivoryBool)
    while true
        # create candidate vector
        v = rand(5)
        if omnivoryBool
            v[5] = rand(Uniform(1.0, PSI_UPPER_BOUND_SOFT))
        else
            v[5] = 1.0
        end
        # halt if valid
        if is_fixed_valid(v)
            return v
        end
    end
end

###############################################################################################################
# SAMPLING FUNCTIONS
###############################################################################################################

# generate a random, valid `fixedParameters` vector that is compatable accross a range of `mainParameters`
function SAMPLE_fixedParameters(omnivoryBool=false)
    # main parameters data
    mainParametersData = generate_mainParameters_data(
        0.1, CANONICAL_MAIN_BASE[1], CANONICAL_MAIN_BASE[2], CANONICAL_MAIN_BASE[3])
    # loop
    while true
        # create candidate, marginally valid vector
        sampleFixedParameters = SAMPLERAW_fixedParameters(omnivoryBool)
        # check it against main data
        valid_accross_main = true
        for row in eachrow(mainParametersData)
            if !is_model_valid(sampleFixedParameters, row, omnivoryBool)
                valid_accross_main = false
            end
        end
        # if valid, return
        if valid_accross_main
            return sampleFixedParameters
        end
    end
    return sampleFixedParameters
end

# generate a random, marginally valid `init` vector
function SAMPLE_init()
    while true
        # sample candidate vector
        sampleInit = rand(Uniform(0.0, RESOURCE_INIT_DENSITY_CAP_SOFT), 6)
        if is_init_valid(sampleInit)
            return sampleInit
        end
    end
end

# map <fixedParams> to two parameter lists, assuming CANONICAL_MAIN_BASE
# Note: <fixedParams> must be the output of SAMPLE_fixedParameters()
function assemble_parameters(fixedParams)
    cₙ, eₙ, μₙ, α, ψ = fixedParams
    e₁, ω, γ = CANONICAL_MAIN_BASE
    # base params
    c₂ = cₙ
    c₃₂ = cₙ
    c₃₁ = cₙ
    e₂ = eₙ
    e₃₂ = eₙ
    e₃₁ = e₁ * ψ          # ψ = e₃₁ / e₃₂
    μ₂₁ = μₙ
    μ₃₂ = μₙ
    μ₃₁ = μ₂₁ * ω         # ω = μ₃₁ / μ₂₁
    z = 4
    baseParams = [c₂, c₃₂, c₃₁, e₁, e₂, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z]
    # resource params
    β = 1 - α
    resourceParams = [α, β, γ]
    # return tuple
    return (baseParams, resourceParams)
end

###############################################################################################################
# SENSITIVITY ANALYSIS
###############################################################################################################

# 1) sample a valid fixedParameters
# 2) for each of the three model types, check if results hold

# TODO: need to define all of the results I want to check
# TODO: start by verifying that Resource response to mainParameters is robust


###############################################################################################################
# TARGET FUNCTIONS (each output T or F)
###############################################################################################################



# in the absence of consumers, checks resource response
function TARGET_resourceResponse_Resource(baseParams, resourceParams, grain, init, timespan)
end