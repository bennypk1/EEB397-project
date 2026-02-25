using DifferentialEquations
using Plots
using DataFrames

#################################################################################################################################
# IMPORT ODE CODE
#################################################################################################################################

function testODE!(dP, P, params, t)
    # parameters and variables
    c₁, c₂₁, c₃₂, c₃₁, e₁, e₂₁, e₃₂, e₃₁, μ₂₁, μ₃₂, μ₃₁, z, pᵤ, pᵤᵤ, a, b, g, I₍₂₃₎, I₍₁₃₎ = params
    P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁ = P
    # dependent dynamic probabilities
    Pₛ = (1 - P₁ - pᵤ)
    P₁ₛ = (P₁ - P₁₁ - Pᵤ₁)
    Q₁ₛ = (P₁ₛ / Pₛ)
    Q₁₁ = (P₁₁ / P₁)
    Pₛᵤ = (pᵤ - Pᵤ₁ - pᵤᵤ)

    # 1) basal resource dynamics
    propagation = (a * c₁ * P₁ₛ)
    seeding = (b * c₁ * Pₛ)
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
    P₁ₛ_propagated_other = P₁ₛ * (a * c₁ * Q₁ₛ * (z - 1) / z)
    P₁ₛ_seeded = P₁ₛ * (b * c₁ * P₁)
    P₁₁_predated2 = P₁₁ * μ₂₁ * (P₍₁₂₎ / P₁)
    P₁₁_predated3 = P₁₁ * μ₃₁ * (P₍₁₃₎ / P₁)
    P₁₁_intrinsic_death = P₁₁ * e₁
    P₁₁_density_dependent_death = P₁₁ * (g * e₁) * ((Q₁₁ * (z - 1) / z) + (1 / z))

    # 6) [moment closure] u-1 dynamics
    Pₛᵤ_propagated = (a * c₁ * (z - 1) / z * Pₛᵤ * Q₁ₛ)
    Pₛᵤ_seeded = Pₛᵤ * (P₁ * b * c₁)
    Pᵤ₁_predated2 = Pᵤ₁ * μ₂₁ * (P₍₁₂₎ / P₁)
    Pᵤ₁_predated3 = Pᵤ₁ * μ₃₁ * (P₍₁₃₎ / P₁)
    Pᵤ₁_intrinsic_death = Pᵤ₁ * e₁
    Pᵤ₁_density_dependent_death = Pᵤ₁ * (g * e₁ * Q₁₁ * (z - 1) / z)

    # [ODE] main dynamics
    dP[1] = propagation + seeding - intrinsic_extinctionR - density_dependent_extinctionR - predationR2 - (I₍₁₃₎) * predationR3
    dP[2] = dispersal21 - intrinsic_extinction21 - density_dependent_extinction21 - (I₍₂₃₎) * predation23
    dP[3] = (I₍₂₃₎) * (creation23_dispersal32 + (I₍₁₃₎) * creation23_dispersal31 + (I₍₁₃₎) * creation23_2outcompetes3 - loss23_intrinsic_extinction - loss23_density_dependent_extinction)
    dP[4] = (I₍₁₃₎) * (creation13_dispersal31 + (I₍₂₃₎) * creation13_dispersal32 + (I₍₂₃₎) * creation13_extinction2 - loss13_intrinsic_extinction - loss13_density_dependent_extinction - creation23_2outcompetes3)
    # [ODE] moment closure dynamics
    dP[5] = P₁ₛ_propagated_direct + P₁ₛ_propagated_other + P₁ₛ_seeded - P₁₁_predated2 - (I₍₁₃₎) * P₁₁_predated3 - P₁₁_intrinsic_death - P₁₁_density_dependent_death
    dP[6] = Pₛᵤ_propagated + Pₛᵤ_seeded - Pᵤ₁_predated2 - (I₍₁₃₎) * Pᵤ₁_predated3 - Pᵤ₁_intrinsic_death - Pᵤ₁_density_dependent_death
end

#################################################################################################################################
# PARAMETETERS
#################################################################################################################################

# Trophic configurations:
#   [0, 0] : consumer-resource
#   [1, 0] : food chain
#   [0, 1] : exploitative competition
#   [1, 1] : omnivory

# Resource Dynamics:
#   [1, 0, 0] :  locally dispersing + NO density-dependent extinction
#   [1, 0, 1] : locally dispersing + density-dependent extinction
#   [0, 1, 1] : globally dispersing + density-dependent extinction
#   [0.5, 0.5, 1] : mixed dispersing + density-dependent extinction

fixed_parameters = [1, 1, 1, 1, 0.05, 0.05, 0.05, 0.05, 0.0025, 0.0025, 0.0025, 4] # c, e, μ, z
initial_densities = [0.01, 0.005, 0.0025, 0.0025, 0.002, 0.005] # starting P₁, P₍₁₂₎, P₍₂₃₎, P₍₁₃₎, P₁₁, Pᵤ₁
resource_dynamics = [1, 0, 0] # a, b, g
trophic_configuration = [1, 1] # I₍₂₃₎, I₍₁₃₎
time_span = [0, 100]
grain = 0.01

#################################################################################################################################
# SINGLE SOLUTION OVER TIME
#################################################################################################################################

this_s = 0.59
this_q = 0.5
conv_params = convert_landscape_params([this_s, this_q])
this_params = vcat(fixed_parameters,        # c, e, μ, z
    conv_params,             # pᵤ, pᵤᵤ
    resource_dynamics,       # a, b, g
    trophic_configuration)   # I₍₂₃₎, I₍₁₃₎
this_init = get_init_from_trophic_config(initial_densities, trophic_configuration)
this_prob = ODEProblem(testODE!, this_init, time_span, this_params)
this_sol = solve(this_prob, AutoVern7(Rodas5()); reltol=1e-8, abstol=1e-10)
labels = ["P₁", "P₍₁₂₎", "P₍₂₃₎", "P₍₁₃₎", "P₁₁", "Pᵤ₁"]
p1 = plot()
for i in 1:length(this_sol.u[1])
    plot!(p1, this_sol.t, getindex.(this_sol.u, i), label=labels[i])
end

#################################################################################################################################
# LIAO-TYPE HEAT MAPS - DATA GENERATION
#################################################################################################################################

# setup
grid = create_unit_grid(grain)
n_sims = length(grid)
data = DataFrame(Availability=Float64[], Connectivity=Float64[],
    DensityR=Float64[], Density2=Float64[], Density3=Float64[],
    Density23=Float64[], Density13=Float64[])

for point in grid
    if !is_valid(point)
        new_row = (point[1], point[2], NaN, NaN, NaN, NaN, NaN)
        push!(data, new_row)
    else
        curr_init = get_init_from_trophic_config(initial_densities, trophic_configuration)
        conv_point = convert_landscape_params([point[1], point[2]])
        curr_params = vcat(fixed_parameters,        # c, e, μ, z
            conv_point,              # pᵤ, pᵤᵤ
            resource_dynamics,       # a, b, g
            trophic_configuration)   # I₍₂₃₎, I₍₁₃₎
        curr_prob = ODEProblem(testODE!, curr_init, time_span, curr_params)
        curr_sol = solve(curr_prob, AutoVern7(Rodas5()); reltol=1e-8, abstol=1e-10)
        sol_end = curr_sol[end]
        new_row = (point[1], point[2], sol_end[1], sol_end[2], sol_end[3] + sol_end[4], sol_end[3], sol_end[4])
        push!(data, new_row)
    end
end

#################################################################################################################################
# LIAO-TYPE HEAT MAPS - DATA VISUALIZATION
#################################################################################################################################

grid_length = 0:grain:1
n = length(grid_length)

Z_R = reshape(data.DensityR, n, n)
Z_2 = reshape(data.Density2, n, n)
Z_3 = reshape(data.Density3, n, n)
Z_23 = reshape(data.Density23, n, n)
Z_13 = reshape(data.Density13, n, n)

pR = heatmap(grid_length, grid_length, Z_R, title="Resource Density (P₁)",
    aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
p2 = heatmap(grid_length, grid_length, Z_2, title="Consumer Density (P₂)",
    aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
p3 = heatmap(grid_length, grid_length, Z_3, title="Tertiary Density (P₃)",
    aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
p23 = heatmap(grid_length, grid_length, Z_23, title="2-3 Link Density (P₍₂₃₎)",
    aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))
p13 = heatmap(grid_length, grid_length, Z_13, title="1-3 Link Density (P₍₁₃₎)",
    aspect_ratio=1, c=cgrad(:jet, scale=:linear), clims=(0, 1))

plot(pR, p2, p3, p1, p23, p13,
    layout=@layout([a b c; d e f]),
    size=(1200, 900))