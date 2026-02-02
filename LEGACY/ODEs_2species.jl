
using DifferentialEquations

# All-Global 2-species predator prey model (Pillai's model)
function predatorPrey2AllGlobal!(dP, init, params, t)
    c₁, c₂, e₁, e₂, μ₂₁, s, ρₛₛ = params
    P₁, P₂, Q₁₁, Qᵤ₁ = init
    dP[1] = c₁ * (s - P₁) * P₁ - (e₁ + μ₂₁) * P₁
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₂ + μ₂₁) * P₂
end

# Local-Prey 2-species predator prey model (Liao's Model ; z = 4)
function predatorPrey2LocalPrey4!(dP, init, params, t)
    c₁, c₂, e₁, e₂, μ₂₁, s, ρₛₛ = params
    P₁, P₂, Q₁₁, Qᵤ₁ = init
    dP[1] = c₁ * (1 - Q₁₁ - Qᵤ₁) * P₁ - (e₁ + μ₂₁) * P₁
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₂ + μ₂₁) * P₂
    dP[3] = (c₁ * (1 - Q₁₁ - Qᵤ₁) * (1/2 + (3/2 * P₁ * (1 - Q₁₁ - Qᵤ₁) / (s - P₁)) - Q₁₁)) - (Q₁₁ * (e₁ + (μ₂₁ * P₂ / P₁)))
    dP[4] = c₁ * (3/4 * ((s - ρₛₛ - P₁*Qᵤ₁) / (s - P₁)) - Qᵤ₁) * (1 - Q₁₁ - Qᵤ₁)
end

# Local-Prey 2-species predator prey model (Liao's Model ; z = 8)
function predatorPrey2Liao!(dP, init, params, t)
    # TODO: function body shown here is incorrect
    c₁, c₂, e₁, e₂, μ₂₁, s, ρₛₛ = params
    P₁, P₂, Q₁₁, Qᵤ₁ = init
    dP[1] = c₁ * (1 - Q₁₁ - Qᵤ₁) * P₁ - (e₁ + μ₂₁) * P₁
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₂ + μ₂₁) * P₂
    dP[3] = (c₁ * (1 - Q₁₁ - Qᵤ₁) * (1/2 + (3/2 * P₁ * (1 - Q₁₁ - Qᵤ₁) / (s - P₁)) - Q₁₁)) - (Q₁₁ * (e₁ + (μ₂₁ * P₂ / P₁)))
    dP[4] = c₁ * 3/4 * (((s - s * (ρₛₛ / s) - P₁ * Qᵤ₁) / (s - P₁)) - Qᵤ₁) * (1 - Q₁₁ - Qᵤ₁)
end