
using DifferentialEquations

##########################################################################################
# Common Parameters: c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ
# Common Varaibles: P₁, P₂, P₃, P₁₁, Pᵤ₁
# Common Neighborhood: z = 4 (other simple possibilities are: z = 2, 3, 6)
##########################################################################################

# Globally Dispersing Single Species Model
function singleGlobal4!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] =  c₁ * P₁ * (s - P₁) - e₁ * P₁
end

# Locally Dispersing Single Species Model
function singleLocal4!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - (e₁ * P₁)
    dP[4] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/4 + 3/4 * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁))) - P₁₁ * e₁
    dP[5] = c₁ * (3/4) * (s - ρ - Pᵤ₁) * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁)) - Pᵤ₁ * e₁
end

# All-Global 2-species predator prey model (Pillai's model)
function predatorPrey2AllGlobal!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] = c₁ * (s - P₁) * P₁ - e₁ * P₁ - μ₂₁ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ + e₂ + μ₂₁) * P₂
end

# Local-Prey 2-species predator prey model (modified Liao)
function predatorPrey2LocalPrey4!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ - μ₂₁ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ + e₂ + μ₂₁) * P₂
    dP[4] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/4 + 3/4 * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁))) - P₁₁ * (e₁ + μ₂₁ * (P₂ / P₁))
    dP[5] = c₁ * (3/4) * (s - ρ - Pᵤ₁) * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁)) - Pᵤ₁ * (e₁ + μ₂₁ * (P₂ / P₁))
end

# Local-Prey, Regional-Predator 2-species predator prey model (modified Liao)
function predatorPrey2LocalPrey4!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ - μ₂₁ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ * (ρ / s) - (e₁ + e₂ + μ₂₁) * P₂
    dP[4] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/4 + 3/4 * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁))) - P₁₁ * (e₁ + μ₂₁ * (P₂ / P₁))
    dP[5] = c₁ * (3/4) * (s - ρ - Pᵤ₁) * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁)) - Pᵤ₁ * (e₁ + μ₂₁ * (P₂ / P₁))
end

# All-Global 3-species predator prey model (Pillai's model)
function predatorPrey3AllGlobal!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] = c₁ * (s - P₁) * P₁ - e₁ * P₁ - μ₂₁ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ + e₂ + μ₂₁) * P₂ - μ₃₂ * P₃
    dP[3] = c₃ * (P₂ - P₃) * P₃ - (e₁ + e₂ + e₃ + μ₂₁ + μ₃₂) * P₃
end
# Local-Prey 3-species predator prey model (modified Liao)
function predatorPrey3LocalPrey4!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ - μ₂₁ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ + e₂ + μ₂₁) * P₂ - μ₃₂ * P₃
    dP[3] = c₃ * (P₂ - P₃) * P₃ - (e₁ + e₂ + e₃ + μ₂₁ + μ₃₂) * P₃
    dP[4] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/4 + 3/4 * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁))) - P₁₁ * (e₁ + μ₂₁ * (P₂ / P₁))
    dP[5] = c₁ * (3/4) * (s - ρ - Pᵤ₁) * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁)) - Pᵤ₁ * (e₁ + μ₂₁ * (P₂ / P₁))
end

# Local-Prey, Regional-Consumer 3-species predator prey model (modified Liao)
function predatorPrey3LocalPrey4!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ - μ₂₁ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ * (ρ / s) - (e₁ + e₂ + μ₂₁) * P₂ - μ₃₂ * P₃
    dP[3] = c₃ * (P₂ - P₃) * P₃ - (e₁ + e₂ + e₃ + μ₂₁ + μ₃₂) * P₃
    dP[4] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/4 + 3/4 * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁))) - P₁₁ * (e₁ + μ₂₁ * (P₂ / P₁))
    dP[5] = c₁ * (3/4) * (s - ρ - Pᵤ₁) * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁)) - Pᵤ₁ * (e₁ + μ₂₁ * (P₂ / P₁))
end

# OTHER

# Local-Prey 2-species predator prey model (Liao's Model ; z = 4)
function predatorPrey2LocalPrey4Liao!(dP, init, params, t)
    c₁, c₂, e₁, e₂, μ₂₁, s, ρₛₛ = params
    P₁, P₂, Q₁₁, Qᵤ₁ = init
    dP[1] = c₁ * (1 - Q₁₁ - Qᵤ₁) * P₁ - (e₁ * P₁) - (μ₂₁ * P₂)
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ + e₂ + μ₂₁) * P₂
    dP[3] = c₁ * (1 - Q₁₁ - Qᵤ₁) * (1/2 + (3/2 * P₁ * (1 - Q₁₁ - Qᵤ₁) / (s - P₁)) - Q₁₁) - Q₁₁ * (e₁ + (μ₂₁ * P₂ / P₁))
    dP[4] = c₁ * (3/4 * ((s - ρₛₛ - P₁*Qᵤ₁) / (s - P₁)) - Qᵤ₁) * (1 - Q₁₁ - Qᵤ₁)
end