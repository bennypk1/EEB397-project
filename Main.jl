using Revise
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "ODEFinal.jl"))


function goofLeg!(dP, P, params, t)
        # parameters and variables
        c₂, c₃₂, c₃₁, e₁, e₂, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, U, F, α, β, γ = params
        P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = P
        # derived variables
        fill!(dP, 0.0)
        S = 1 - U
        pᵤ = U
        pᵤᵤ = F
        p₀ = 1 - P₁ - pᵤ
        q₁₁ = P₁₁ / P₁
        qᵤ₁ = Pᵤ₁ / P₁
        p₀₁ = P₁ - P₁₁ - Pᵤ₁
        pᵤ₀ = pᵤ - pᵤᵤ - Pᵤ₁
        q₁₀ = p₀₁ / p₀
        # main ODEs
        dP[1] = α * (1 - q₁₁ - qᵤ₁) * P₁ +
                β * (S - P₁) * P₁ -
                γ * q₁₁ * P₁ -
                e₁ * P₁
        # moment closure
        dP[5] = -e₁ * P₁₁ -
                γ * P₁₁ * (1 / z + (z - 1) / z * q₁₁) +
                α * p₀₁ * (1 / z + (z - 1) / z * q₁₁) +
                β * p₀₁ * P₁
        dP[6] = -e₁ * Pᵤ₁ -
                γ * Pᵤ₁ * ((z - 1) / z * q₁₁) +
                α * pᵤ₀ * ((z - 1) / z * q₁₀) +
                β * pᵤ₀ * P₁
end

################################################################################################################

params = [CANONICAL_PARAMS; [0.5, 0.4]; [1, 0, 0]]
init = CANONICAL_INIT
timespan = CANONICAL_TIMESPAN
plotRun(goofLeg!, params, init, [0, 100])

# grain = 0.05
# y = LiaoTypeGrid1(goofLeg!, CANONICAL_PARAMS, [0, 1, 0], grain, init, timespan)
# LiaoTypeHeatMap(goofLeg!, CANONICAL_PARAMS, [1, 0, 0], grain, init, timespan)