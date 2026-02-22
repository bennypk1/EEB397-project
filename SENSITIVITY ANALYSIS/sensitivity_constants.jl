
# This file holds constants (mainly hard/soft marginal parameter boundaries) require to generate biologically
# plausible parmeters, and test qualitative result robustness to different parametric combinations. Also
# contains descriptive information about the parmeters when relevant.


# INDEPENDENT NUISANCE PARAMETERS
# cₙ  = c₂ = c₃₂ = c₃₁   (consumer colonization rate)
# eₙ  = e₂ = e₃₂ = e₃₁   (consumer extinction rate)
# μₙ  = μ₃₁ = μ₃₂        (species 3 predation rate)
# α                      (resource local dispersal rate)

# MAIN PARAMETERS
# e₁                     (resource intrinsic extinction rate)
# γ                      (resource local crowding sensitivity)
# δ                      (differential predation pressure: larger => μ₂₁ > μ₃₁)  

# DEPENENT PARAMETERS
# β = 1 - α              (resource global dispersal rate)
# z = 4                  (von-Neuman neighborship degree)
# μ₂₁ = δ * μₙ           (species 2 predation rate)


# SOFT CONSTRAINTS

SOFT_Cn_MAX = 1.0
SOFT_E1_MIN = 0.0025
SOFT_E1_MAX = 0.9
SOFT_DELTA_MIN = 0.1
SOFT_DELTA_MAX = 10
SOFT_GAMMA_MAX = 2
SEARCH_GRAIN = 10

SENSITIVITY_MAIN_BASE = [0.05, 0, 1] # Note: these are the Liao ones. In target functions should override e₁ when using Resource-only model

# POTENTIAL OTHER THINGS NOT CODED RIGHT NOW
# 1) increased extinction of species 3 when feeding on species 1 (obviously, only in Omnivory)
# 2) global dispersal of species 2 (Note: this is only interesting if applied in a Simple Food Chain)


# OTHER NOTES
# 1) δ = 1 / ω (where ω is define in Liao 2017)