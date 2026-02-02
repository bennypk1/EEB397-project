
function everythingCR!(dP, P, params, t)
    c₁, c₂, e₁, e₂, μ, pᵤ, pᵤᵤ, z, a, b, g = params
    P₁, P₂, P₁₁, Pᵤ₁ = P
    # probabilities
    Pₛ = (1 - P₁ - pᵤ)
    P₁ₛ = (P₁ - P₁₁ - Pᵤ₁)
    Q₁ₛ =  P₁ₛ / Pₛ
    Q₁₁ = (P₁₁ / P₁)
    Pₛᵤ = (pᵤ - Pᵤ₁ - pᵤᵤ)
    # species dynamics
    propagation = (a * c₁ * Pₛ * Q₁ₛ)
    seeding = (b * c₁ * Pₛ)
    intrinsic_mortalityR = (e₁ * P₁)
    density_dependent_mortalityR = P₁ * (g * e₁ * Q₁₁) # if g < 0, this term represents decreased mortality as a function of crowding
    predation = (μ * P₂)
    dispersalC = c₂ * (P₁ - P₂) * P₂
    density_independent_mortalityC = (e₂ + μ + e₁) * P₂
    density_dependent_mortalityC = P₂ * (g * e₁ * Q₁₁) # only exists when density_dependent_mortalityR != 0
    # pair dynanmics 1-1
    P₁ₛ_propagated_direct = P₁ₛ * (a * c₁ / z)
    P₁ₛ_propagated_other = P₁ₛ * (a * c₁ * Q₁ₛ * (z-1)/z)
    P₁ₛ_seeded = P₁ₛ * (b * c₁ * P₁)
    P₁₁_predated = P₁₁ * μ * (P₂ / P₁)
    P₁₁_intrinsic_death = P₁₁ * e₁
    P₁₁_density_dependent_death = P₁₁ * (g * e₁) * ((Q₁₁ * (z-1)/z) + (1/z))
    # pair dynamics u-1
    Pₛᵤ_propagated = (a * c₁ * (z-1)/z * Pₛᵤ * Q₁ₛ)
    Pₛᵤ_seeded = Pₛᵤ * (P₁ * b * c₁)
    Pᵤ₁_predated = Pᵤ₁ * μ * (P₂ / P₁)
    Pᵤ₁_intrinsic_death = Pᵤ₁ * e₁
    Pᵤ₁_density_dependent_death = Pᵤ₁ * (g * e₁ * Q₁₁ * (z-1)/z)
    # ODE
    dP[1] = propagation + seeding - intrinsic_mortalityR - density_dependent_mortalityR - predation
    dP[2] = dispersalC - density_independent_mortalityC - density_dependent_mortalityC
    dP[3] = P₁ₛ_propagated_direct +  P₁ₛ_propagated_other + P₁ₛ_seeded - P₁₁_predated - P₁₁_intrinsic_death - P₁₁_density_dependent_death
    dP[4] = Pₛᵤ_propagated + Pₛᵤ_seeded - Pᵤ₁_predated - Pᵤ₁_intrinsic_death - Pᵤ₁_density_dependent_death
end