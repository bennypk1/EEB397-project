using Revise
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))

function goofLeg!(dP, P, params, t)
    # parameters and variables
    c₂, c₃₂, c₃₁, e₁, e₂, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, U, F, α, β, γ = params
    P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = P
    # derived variables
    fill!(dP, 0.0)
    S = 1 - U
    pᵤ = U
    pᵤᵤ = 1 - 2S + (1 - F) * S # THIS IS A VERY IMPORTANT RELATIONSHIP
    p₀ = 1 - P₁ - pᵤ
    q₁₁ = P₁₁ / P₁
    qᵤ₁ = Pᵤ₁ / P₁
    p₀₁ = P₁ - P₁₁ - Pᵤ₁
    pᵤ₀ = pᵤ - pᵤᵤ - Pᵤ₁
    q₁₀ = p₀₁ / p₀
    # main ODEs
    dP[1] = α * (1 - q₁₁ - qᵤ₁) * P₁ - e₁ * P₁
    # moment closure
    dP[5] = -e₁ * P₁₁ + α * p₀₁ * ((1 / z) + (z - 1) / z * q₁₀)
    dP[6] = -e₁ * Pᵤ₁ + α * pᵤ₀ * ((z - 1) / z * q₁₀)
end

################################################################################################################

params = [CANONICAL_PARAMS; [0.5, 0.5]; [0, 1, 0]]
init = CANONICAL_INIT
timespan = CANONICAL_TIMESPAN
# plotRun(ExploitativeCompetition3!, params, init, [0, 100])

grain = 0.005
CANONICAL_PARAMS[4] = 0.05 # set intrinsic extinction of resource
CANONICAL_PARAMS[8] = 0.025 # set additive mortality rate due to 2 feeding on 1
CANONICAL_PARAMS[10] = 0.025 # set additive mortality rate due to 3 feeding on 1

LiaoTypeSpeciesRichnessMap(SimpleFoodChain!, CANONICAL_PARAMS, [1, 0, 0], grain, init, timespan)
# testGrid = LiaoTypeGridExtra(ExploitativeCompetition!, CANONICAL_PARAMS, [0, 1, 0.5], grain, init, timespan)
# assignSpeciesDistributionID(testGrid)

# TODO: NEED TO WRITE CHECKS AND ERROR THROWS FOR LiaoTypeGrid TO MAKE SURE ALL PERSISTANCES ARE ALWAYS POSITIVE
