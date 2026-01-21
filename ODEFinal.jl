# RESOURCE DYNAMICS


# CONSUMER DYNAMICS 
# 1) top-down predation: species 2 can only eat 1, and 3 can only eat 2
# 2) explotative competition: 2 outcompetes 3 always, but 3 is a better disperser
# 3) omnivory: 2 outcompetes 3 when feeding on 1 ; 3 switches to 1 if 2 goes extinct ; 3 switches to 1 if 2 comes to outcompete it
# 4) species cannot colonize a patch where itself is already there
# 5) species can only exist in patches where their food is present

function Resource!(dP, P, params, t)
    # parameters and variables
    c₂, c₃₂, c₃₁, e₁, e₂, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, U, F, α, β, γ = params
    P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = P
    # derived variables
    S = 1 - U
    pᵤ = U
    pᵤᵤ = F
    p₀ = 1 - P₁ - pᵤ
    q₁₁ = P₁₁ / P₁
    qᵤ₁ = Pᵤ₁ / P₁
    p₀₁ = P₁ - P₁₁ -  Pᵤ₁
    pᵤ₀ = pᵤ - pᵤᵤ - Pᵤ₁
    q₁₀ = p₀₁ / p₀
    # main ODEs
    dP[1] = α * ( 1 - q₁₁ - qᵤ₁ ) * P₁ 
          + β * ( S - P₁ ) * P₁ 
          - γ * q₁₁ * P₁ 
          - e₁ * P₁
    # moment closure
    dP[5] = 
          - 2e₁ * P₁₁ 
          - 2γ * P₁₁ * ( 1/z + (z-1)/z * q₁₁ ) 
          + 2α * p₀₁ * ( 1/z + (z-1)/z * q₁₁ ) 
          + 2β * p₀₁ * P₁
    dP[6] = 
          - e₁ * Pᵤ₁ 
          - γ * Pᵤ₁ * ( (z-1)/z * q₁₁ ) 
          + α * pᵤ₀ * ( (z-1)/z * q₁₀ ) 
          + β * pᵤ₀ * P₁
end

function SimpleFoodChain!(dP, P, params, t)
    # parameters and variables
    c₂, c₃₂, c₃₁, e₁, e₂, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, U, F, α, β, γ = params
    P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = P
    # derived variables
    S = 1 - U
    pᵤ = U
    pᵤᵤ = F
    p₀ = 1 - P₁ - pᵤ
    q₁₁ = P₁₁ / P₁
    qᵤ₁ = Pᵤ₁ / P₁
    p₀₁ = P₁ - P₁₁ -  Pᵤ₁
    pᵤ₀ = pᵤ - pᵤᵤ - Pᵤ₁
    q₁₀ = p₀₁ / p₀
    # main ODEs
    dP[1] = α * ( 1 - q₁₁ - qᵤ₁ ) * P₁ 
          + β * ( S - P₁ ) * P₁ 
          - γ * q₁₁ * P₁ 
          - e₁ * P₁
          - μ₂₁ * P₍₁₂₎
    dP[2] = c₂ * P₍₁₂₎ * (P₁ -  P₍₁₂₎)
          - (γ * q₁₁ + e₁ + e₂ + μ₂₁) * P₍₁₂₎
          - μ₃₂ * P₍₂₃₎
    dP[3] = c₃₂ * P₍₂₃₎ * (P₍₁₂₎ - P₍₂₃₎) 
          - (γ * q₁₁ + e₁ + e₂ + e₃₂ + μ₂₁ + μ₃₂) * P₍₂₃₎
    # moment closure
    dP[5] = - 2e₁ * P₁₁ 
          - 2γ * P₁₁ * ( 1/z + (z-1)/z * q₁₁ ) 
          + 2α * p₀₁ * ( 1/z + (z-1)/z * q₁₁ ) 
          + 2β * p₀₁ * P₁
          - (μ₂₁ * P₍₁₂₎ / P₁) * P₁₁
    dP[6] = - e₁ * Pᵤ₁ 
          - γ * Pᵤ₁ * ( (z-1)/z * q₁₁ ) 
          + α * pᵤ₀ * ( (z-1)/z * q₁₀ ) 
          + β * pᵤ₀ * P₁
          - (μ₂₁ * P₍₁₂₎ / P₁) * P₁₁
end

function CompetitiveExclusion!(dP, P, params, t)
    # parameters and variables
    c₂, c₃₂, c₃₁, e₁, e₂, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, U, F, α, β, γ = params
    P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = P
    # derived variables
    S = 1 - U
    pᵤ = U
    pᵤᵤ = F
    p₀ = 1 - P₁ - pᵤ
    q₁₁ = P₁₁ / P₁
    qᵤ₁ = Pᵤ₁ / P₁
    p₀₁ = P₁ - P₁₁ -  Pᵤ₁
    pᵤ₀ = pᵤ - pᵤᵤ - Pᵤ₁
    q₁₀ = p₀₁ / p₀
    # main ODEs
    dP[1] = α * ( 1 - q₁₁ - qᵤ₁ ) * P₁ 
          + β * ( S - P₁ ) * P₁ 
          - γ * q₁₁ * P₁ 
          - e₁ * P₁
          - μ₂₁ * P₍₁₂₎
          - μ₃₁ * P₍₁₃₎
    dP[2] = c₂ * P₍₁₂₎ * (P₁ -  P₍₁₂₎)             # colonization + competitive exclusion
          - (γ * q₁₁ + e₁ + e₂ + μ₂₁) * P₍₁₂₎      # extibctions + predations
    dP[4] = c₃₁ * P₍₁₃₎ * (P₁ - P₍₁₂₎ - P₍₁₃₎)     # colonization
          - (γ * q₁₁ + e₁ + e₃₁ + μ₃₁) * P₍₁₃₎     # extinctions + predations   
          - c₂ * P₍₁₂₎ * P₍₁₃₎                     # competition
    # moment closure
    dP[5] = - 2e₁ * P₁₁ 
          - 2γ * P₁₁ * ( 1/z + (z-1)/z * q₁₁ ) 
          + 2α * p₀₁ * ( 1/z + (z-1)/z * q₁₁ ) 
          + 2β * p₀₁ * P₁
          - (μ₂₁ * P₍₁₂₎ / P₁) * P₁₁
          - (μ₃₁ * P₍₁₃₎ / P₁) * P₁₁
    dP[6] = - e₁ * Pᵤ₁ 
          - γ * Pᵤ₁ * ( (z-1)/z * q₁₁ ) 
          + α * pᵤ₀ * ( (z-1)/z * q₁₀ ) 
          + β * pᵤ₀ * P₁
          - (μ₂₁ * P₍₁₂₎ / P₁) * P₁₁
          - (μ₃₁ * P₍₁₃₎ / P₁) * P₁₁
end

function Omnivory!(dP, P, params, t)
    # parameters and variables
    c₂, c₃₂, c₃₁, e₁, e₂, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, U, F, α, β, γ = params
    P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = P
    # derived variables
    S = 1 - U
    pᵤ = U
    pᵤᵤ = F
    p₀ = 1 - P₁ - pᵤ
    q₁₁ = P₁₁ / P₁
    qᵤ₁ = Pᵤ₁ / P₁
    p₀₁ = P₁ - P₁₁ -  Pᵤ₁
    pᵤ₀ = pᵤ - pᵤᵤ - Pᵤ₁
    q₁₀ = p₀₁ / p₀
    # main ODEs
    dP[1] = α * ( 1 - q₁₁ - qᵤ₁ ) * P₁ 
          + β * ( S - P₁ ) * P₁ 
          - γ * q₁₁ * P₁ 
          - e₁ * P₁
          - μ₂₁ * P₍₁₂₎
          - μ₃₁ * P₍₁₃₎
    dP[2] = c₂ * P₍₁₂₎ * (P₁ -  P₍₁₂₎)                             # colonization + competitive exclusion
          - (γ * q₁₁ + e₁ + e₂ + μ₂₁) * P₍₁₂₎                      # extibctions + predations
          - μ₃₂ * P₍₂₃₎                                            # predation
    dP[3] = (c₃₁ * P₍₁₃₎ + c₃₂ * P₍₂₃₎) * (P₍₁₂₎ - P₍₁₃₎)          # colonization
          - (γ * q₁₁ + e₁ + e₂ + e₃₂ + μ₂₁ + μ₃₂) * P₍₂₃₎          # extinctions + predations 
    dP[4] = (c₃₁ * P₍₁₃₎ + c₃₂ * P₍₂₃₎) * (P₁ - P₍₁₂₎ - P₍₁₃₎)     # colonization
          + (e₂ + μ₃₂) * P₍₂₃₎                                     # prey switching
          - (γ * q₁₁ + e₁ + e₃₁ + μ₃₁) * P₍₁₃₎                     # extinctions + predations   
          - c₂ * P₍₁₂₎ * P₍₁₃₎                                     # competitive exclusion
    # moment closure
    dP[5] = - 2e₁ * P₁₁ 
          - 2γ * P₁₁ * ( 1/z + (z-1)/z * q₁₁ ) 
          + 2α * p₀₁ * ( 1/z + (z-1)/z * q₁₁ ) 
          + 2β * p₀₁ * P₁
          - (μ₂₁ * P₍₁₂₎ / P₁) * P₁₁
          - (μ₃₁ * P₍₁₃₎ / P₁) * P₁₁
    dP[6] = - e₁ * Pᵤ₁ 
          - γ * Pᵤ₁ * ( (z-1)/z * q₁₁ ) 
          + α * pᵤ₀ * ( (z-1)/z * q₁₀ ) 
          + β * pᵤ₀ * P₁
          - (μ₂₁ * P₍₁₂₎ / P₁) * P₁₁
          - (μ₃₁ * P₍₁₃₎ / P₁) * P₁₁
end