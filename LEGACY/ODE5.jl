
# CONSUMER DYNAMICS 
# 1) top-down predation:
# 2) explotative competition: 2 outcompetes 3
# 3) omnivory: 2 outcompetes 3 when feeding on 1 ; 3 switches to 1 if 2 goes extinct ; 3 switches to 1 if 2 comes to outcompete it

# Note: species cannot colonize a patch where itself is already there


# cⱼᵢ   : POSITIVE REAL : dispersal rate* of species j when feeding on species i
# eⱼᵢ   : PROPORTION    : intrinsic extinction likelihood of species j when feeding on species i
# μⱼᵢ   : POSITIVE REAL : added extinction likelihood due to top-down predation of species j on species i
# P₍ᵢⱼ₎ : PROPORTION    : proportion of patches where a feeding link from species i to speciesj exists
# I₍ᵢⱼ₎ : BINARY        : parameter inducating whether a feeding link exists between species i and species j

# z     : POSITIVE INT  : number of nearest neighbors a single patch has (z = 4 for entire paper)
# a     : POSITIVE REAL : rate at which resources allocalted to dispersal successfully establish new populations locally
# b     : POSITIVE REAL : rate at which resources allocalted to dispersal successfully establish new populations globally
# g     : POSITIVE REAL : sensitivity to crowding (if negative: survival benefit of crowding)

function everything!(dP, P, params, t)
    # parameters and variables
    c₁, c₂₁ , c₃₂, c₃₁, e₁, e₂₁, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, pᵤ, pᵤᵤ, a, b, g, I₍₂₃₎, I₍₁₃₎ = params
    P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = P
    # dependent dynamic probabilities
    Pₛ = (1 - P₁ - pᵤ)
    P₁ₛ = (P₁ - P₁₁ - Pᵤ₁)
    Q₁ₛ = (P₁ₛ / Pₛ)
    Q₁₁ = (P₁₁ / P₁)
    Pₛᵤ = (pᵤ - Pᵤ₁ - pᵤᵤ)

    # 1) basal resource dynamics
    propagation = (a * c₁ * P₁ₛ)
    seeding = (b * c₁ * Pₛ) * P₁
    intrinsic_extinctionR = (e₁ * P₁)
    density_dependent_extinctionR = P₁ * (g * e₁ * Q₁₁)
    predationR2 = (μ₂₁ * P₍₁₂₎)
    predationR3 = (μ₃₁ * P₍₁₃₎)

    # 2) 1-2 dynamics
    dispersal21 = c₂₁ * (P₁ - P₍₁₂₎) * P₍₁₂₎ # takes into account gain from competitive exclusions of 3 in 3 -> 1 patches by 2
    intrinsic_extinction21 = (e₂₁ + e₁ + μ₂₁) * P₍₁₂₎
    density_dependent_extinction21 = P₍₁₂₎ * (g * e₁ * Q₁₁)
    predation23 = (μ₃₂ * P₍₂₃₎)

    # 3) 2-3 dynamics
    creation23_dispersal32 = c₃₂ * P₍₂₃₎ * (P₍₁₂₎ - P₍₂₃₎)
    creation23_dispersal31 = c₃₁ * P₍₁₃₎ * (P₍₁₂₎ - P₍₂₃₎)
    creation23_2outcompetes3 = c₂₁ * P₍₁₂₎ * P₍₁₃₎ # in 3 -> 1 patches, when 2 outcompetes 3, 3 switches to 2
    loss23_intrinsic_extinction = (e₁ + e₂₁ + e₃₂ + μ₂₁ + μ₃₂) * P₍₂₃₎
    loss23_density_dependent_extinction = P₍₂₃₎ * (g * e₁ * Q₁₁)
    
    # 4) 1-3 dynamics
    creation13_dispersal32 = c₃₂ * P₍₂₃₎ * (P₁ - P₍₁₂₎ - P₍₁₃₎)
    creation13_dispersal31 = c₃₁ * P₍₁₃₎ * (P₁ - P₍₁₂₎ - P₍₁₃₎)
    creation13_extinction2 = (e₂₁ + μ₂₁) * P₍₂₃₎ # in 3 -> 2 -> 1 patches, when just 2 goes extinct and 3 switches to 1
    loss13_intrinsic_extinction = (e₁ + e₃₁ + μ₃₁) * P₍₁₃₎
    loss13_density_dependent_extinction = P₍₁₃₎ * (g * e₁ * Q₁₁)
    
    # 5) [moment closure] 1-1 dynamics
    P₁ₛ_propagated_direct = P₁ₛ * (a * c₁ / z)
    P₁ₛ_propagated_other = P₁ₛ * (a * c₁ * Q₁ₛ * (z-1)/z)
    P₁ₛ_seeded = P₁ₛ * (b * c₁ * P₁)
    P₁₁_predated2 = P₁₁ * μ₂₁ * (P₍₁₂₎ / P₁)
    P₁₁_predated3 = P₁₁ * μ₃₁ * (P₍₁₃₎ / P₁)
    P₁₁_intrinsic_death = P₁₁ * e₁
    P₁₁_density_dependent_death = P₁₁ * (g * e₁) * ((Q₁₁ * (z-1)/z) + (1/z))

    # 6) [moment closure] u-1 dynamics
    Pₛᵤ_propagated = (a * c₁ * (z-1)/z * Pₛᵤ * Q₁ₛ)
    Pₛᵤ_seeded = Pₛᵤ * (P₁ * b * c₁)
    Pᵤ₁_predated2 = Pᵤ₁ * μ₂₁ * (P₍₁₂₎ / P₁)
    Pᵤ₁_predated3 = Pᵤ₁ * μ₃₁ * (P₍₁₃₎ / P₁)
    Pᵤ₁_intrinsic_death = Pᵤ₁ * e₁
    Pᵤ₁_density_dependent_death = Pᵤ₁ * (g * e₁ * Q₁₁ * (z-1)/z)

    # [ODE] main dynamics
    dP[1] = propagation + seeding - intrinsic_extinctionR - density_dependent_extinctionR - predationR2 - (I₍₁₃₎) * predationR3
    dP[2] = dispersal21 - intrinsic_extinction21 - density_dependent_extinction21 - (I₍₂₃₎) * predation23
    dP[3] = creation23_dispersal32 + (I₍₁₃₎) * creation23_dispersal31 + (I₍₁₃₎) * creation23_2outcompetes3 - loss23_intrinsic_extinction - loss23_density_dependent_extinction
    dP[4] = creation13_dispersal31 + (I₍₂₃₎) * creation13_dispersal32 + (I₍₂₃₎) * creation13_extinction2 - loss13_intrinsic_extinction - loss13_density_dependent_extinction - creation23_2outcompetes3
    # [ODE] moment closure dynamics
    dP[5] = P₁ₛ_propagated_direct +  P₁ₛ_propagated_other + P₁ₛ_seeded - P₁₁_predated2 - (I₍₁₃₎) * P₁₁_predated3 - P₁₁_intrinsic_death - P₁₁_density_dependent_death
    dP[6] = Pₛᵤ_propagated + Pₛᵤ_seeded - Pᵤ₁_predated2 - (I₍₁₃₎) * Pᵤ₁_predated3 - Pᵤ₁_intrinsic_death - Pᵤ₁_density_dependent_death
end