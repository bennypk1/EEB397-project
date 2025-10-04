
using DifferentialEquations

# Pillai's 2-species predator prey model
function predatorPrey2Pillai!(dP, init, params, t)
    c₁, c₂, e₁, e₂, μ₂₁ = params
    P₁, P₂ = init
    dP[1] = c₁ * (1 - P₁) * P₁ - (e₁ + μ₂₁) * P₁
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₂ + μ₂₁) * P₂
end

# Liao's 2-species predator prey model
function predatorPrey2Liao!(dP, init, params, t)
    c₁, c₂, e₁, e₂, μ₂₁, s, ρₛₛ, z = params
    P₁, P₂, Q₁₁, Qᵤ₁ = init
    dP[1] = c₁ * (1 - Q₁₁ - Qᵤ₁) * P₁ - (e₁ + μ₂₁) * P₁
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₂ + μ₂₁) * P₂
end
