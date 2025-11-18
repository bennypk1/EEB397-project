using DifferentialEquations
using Plots
using DataFrames
include("helper_functions.jl")

# Replace function body with one of the ODEs
function testODE!(dP, P, params, t)
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

# Single solution over time
this_s = 0.59
this_q = 0.5
conv_params = convert_landscape_params([this_s, this_q])
fixed_params = [1, 1, 0.05, 0.05, 0.0025, conv_params[1], conv_params[2], 4, 1, 0.1, 0]
tm_spn = [0, 100]
initial_density = [0.01, 0.005, 0.002, 0.005]
this_prob = ODEProblem(testODE!, initial_density, tm_spn, fixed_params)
this_sol = solve(this_prob, AutoVern7(Rodas5()); reltol=1e-8, abstol=1e-10)
p1 = plot(this_sol)

# Full Simulation
grain = 0.01
grid = create_unit_grid(grain)
n_sims = length(grid)
data = DataFrame(Availability = Float64[], Connectivity = Float64[], ResDensity = Float64[], ConsDensity = Float64[])
wavy_data = DataFrame(Availability = Float64[], Connectivity = Float64[], ResRange = Float64[])

for point in grid
    if !is_valid(point)
        new_row = (point[1], point[2], NaN, NaN)
        push!(data, new_row)
        new_wavy_row = (point[1], point[2], NaN)
        push!(wavy_data, new_wavy_row)
    else
        # save data for colorful Liao-type plots
        conv_point = convert_landscape_params([point[1], point[2]])
        curr_params = [1, 1, 0.05, 0.05, 0.0025, conv_point[1], conv_point[2], 4, 1, 0, 0]
        curr_prob = ODEProblem(testODE!, initial_density, tm_spn, curr_params)
        curr_sol = solve(curr_prob, AutoVern7(Rodas5()); reltol=1e-8, abstol=1e-10)
        sol_end = curr_sol[end]
        new_row = (point[1], point[2], sol_end[1], sol_end[2])
        push!(data, new_row)
        
        # save data to analyze oscillations
        last_times = [t for t in curr_sol.t if t ≥ tm_spn[end] - 10]
        last_vals = [curr_sol(t)[1] for t in last_times]
        new_wavy_row = (point[1], point[2], maximum(last_vals) - minimum(last_vals))
        push!(wavy_data, new_wavy_row)
    end
end

# Heat Map Plots
grid_length = 0:grain:1
z_matrix2 = reshape(data.ResDensity, length(grid_length), length(grid_length))
z_matrix3 = reshape(data.ConsDensity, length(grid_length), length(grid_length))
z_matrix4 = reshape(wavy_data.ResRange, length(grid_length), length(grid_length))
p2 = heatmap(grid_length, grid_length, z_matrix2,
    xlabel="Availability", ylabel="Connectivity", colorbar_title=" Resource Density",
    aspect_ratio = 0.85, c = cgrad(:jet, scale = :linear), clims = (0, 1)
)
p3 = heatmap(grid_length, grid_length, z_matrix3,
    xlabel="Availability", ylabel="Connectivity", colorbar_title="Consumer Density",
    aspect_ratio = 0.85, c = cgrad(:jet, scale = :linear), clims = (0, 1)
)
p4 = heatmap(grid_length, grid_length, z_matrix4,
    xlabel="Availability", ylabel="Connectivity", colorbar_title="Resource Wavy-ness",
    aspect_ratio = 0.85, c = cgrad(:jet, scale = :linear)
)
plot(p1, p2, p3, p4, size = (800, 800))