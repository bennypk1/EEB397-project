# Variations on the basic contact process for single species

##############################################################
# Common Parameters: c₁, c₂, e₁, e₂, μ, pᵤ, pᵤᵤ, z
# Common Varaibles: P₁, P₂, P₃, P₁₁, Pᵤ₁
##############################################################

# Globally Consumer AND Resource (equivalent to Pillai model)
function allGlobal!(dP, P, params, t)
    c₁, c₂, e₁, e₂, μ, pᵤ, pᵤᵤ, z = params
    P₁, P₂, P₁₁, Pᵤ₁ = P
    dP[1] =  c₁ * P₁ * (1 - pᵤ - P₁) - e₁ * P₁ - μ * P₂
    dP[2] =  c₂ * P₂ * (P₁ - P₂)  - (e₁ + e₂ + μ) * P₂
end

# Global Consumer, Locally Dispersing Resource (equivalent to Liao 2016 model)
function partialLocalResource!(dP, P, params, t)
    c₁, c₂, e₁, e₂, μ, pᵤ, pᵤᵤ, z = params
    P₁, P₂, P₁₁, Pᵤ₁ = P
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ - μ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ + e₂ + μ) * P₂
    dP[3] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/z + (z-1)/z * ((P₁ - P₁₁ - Pᵤ₁) / (1 - pᵤ - P₁))) - P₁₁ * (e₁ + μ * P₂ / P₁)
    dP[4] = c₁ * ((z-1)/z) * (pᵤ - Pᵤ₁ - pᵤᵤ) * ((P₁ - P₁₁ - Pᵤ₁) / (1 - pᵤ - P₁)) - Pᵤ₁ * (e₁ + μ * P₂ / P₁)
end

# Global Consumer, Locally Dispersing + Negative Density-Dependent Mortality Resource (inspired by Ellner)
# Density-Dependent Mortality Rule: e = e₁ * (1 + clustering)
function negativeLocalResource!(dP, P, params, t)
    c₁, c₂, e₁, e₂, μ, pᵤ, pᵤᵤ, z = params
    P₁, P₂, P₁₁, Pᵤ₁ = P
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ - P₁₁ * e₁ - μ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ * P₂ * (1 + (P₁₁ / P₁)) + (e₂ + μ) * P₂)
    dP[3] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/z + (z-1)/z * ((P₁ - P₁₁ - Pᵤ₁) / (1 - pᵤ - P₁))) - P₁₁ * (e₁ * (1 + 1/z + ((z-1)/z) * (P₁₁ / P₁)) + μ * P₂ / P₁)
    dP[4] = c₁ * ((z-1)/z) * (pᵤ - Pᵤ₁ - pᵤᵤ) * ((P₁ - P₁₁ - Pᵤ₁) / (1 - pᵤ - P₁)) - Pᵤ₁ * (e₁ * (1 + (P₁₁ / P₁) * ((z-1)/z)) + μ * P₂ / P₁)
end

# Global Consumer, Locally Dispersing + Positive Density-Dependent Mortality Resource (inspired by Ellner)
# Density-Dependent Mortality Rule: e = e₁ * (1 - clustering)
function positiveLocalResource!(dP, P, params, t)
    c₁, c₂, e₁, e₂, μ, pᵤ, pᵤᵤ, z = params
    P₁, P₂, P₁₁, Pᵤ₁ = P
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ + P₁₁ * e₁ - μ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ * P₂ * (1 - (P₁₁ / P₁)) + (e₂ + μ) * P₂)
    dP[3] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/z + (z-1)/z * ((P₁ - P₁₁ - Pᵤ₁) / (1 - pᵤ - P₁))) - P₁₁ * (e₁ * (1 - 1/z - ((z-1)/z) * (P₁₁ / P₁)) + μ * P₂ / P₁)
    dP[4] = c₁ * ((z-1)/z) * (pᵤ - Pᵤ₁ - pᵤᵤ) * ((P₁ - P₁₁ - Pᵤ₁) / (1 - pᵤ - P₁)) - Pᵤ₁ * (e₁ * (1 - (P₁₁ / P₁) * ((z-1)/z)) + μ * P₂ / P₁)
end