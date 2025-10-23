using DifferentialEquations
using Plots
using DataFrames
include("helper_functions.jl")

# Replace function body with one of the ODEs in plot_ODE1.jl
function testODE!(dP, init, params, t)
    c₁, c₂, c₃, e₁, e₂, e₃, μ₂₁, μ₃₂, s, ρ = params
    P₁, P₂, P₃, P₁₁, Pᵤ₁ = init
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ - μ₂₁ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ * (ρ / s) - (e₁ + e₂ + μ₂₁) * P₂
    dP[4] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/4 + 3/4 * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁))) - P₁₁ * (e₁ + μ₂₁ * (P₂ / P₁))
    dP[5] = c₁ * (3/4) * (s - ρ - Pᵤ₁) * ((P₁ - P₁₁ - Pᵤ₁) / (s - P₁)) - Pᵤ₁ * (e₁ + μ₂₁ * (P₂ / P₁))
end

# Single solution over time
this_s = 0.6
this_q = 0.75
fixed_params = [1, 1, 1, 0.05, 0.05, 0.05, 0.0025, 0.0025, this_s, this_s * this_q]
tm_spn = [0, 100]
initial_density = [0.01, 0.005, 0.002, 0.002, 0.005]
this_prob = ODEProblem(testODE!, initial_density, tm_spn, fixed_params)
this_sol = solve(this_prob)
p1 = plot(this_sol)

# Full Simulation
grain = 0.015
grid = create_unit_grid(grain)
n_sims = length(grid)
data = DataFrame(Availability = Float64[], Connectivity = Float64[], PreyDensity = Float64[], PredDensity = Float64[], Pred2Density = Float64[])

for point in grid
    if !is_valid(point)
        new_row = (point[1], point[2], NaN, NaN, NaN)
        push!(data, new_row)
    else
    curr_params = [1, 1, 1, 0.05, 0.05, 0.05, 0.0025, 0.0025, point[1], point[1] * point[2]]
    curr_prob = ODEProblem(testODE!, initial_density, tm_spn, curr_params)
    curr_sol = solve(curr_prob)
    sol_end = curr_sol[end]
    new_row = (point[1], point[2], sol_end[1], sol_end[2], sol_end[3])
    push!(data, new_row)
    end
end

# Heat Map Plots
grid_length = 0:grain:1
z_matrix2 = reshape(data.PreyDensity, length(grid_length), length(grid_length))
z_matrix3 = reshape(data.PredDensity, length(grid_length), length(grid_length))
z_matrix4 = reshape(data.Pred2Density, length(grid_length), length(grid_length))
p2 = heatmap(grid_length,grid_length, z_matrix2,
    xlabel="Availability", ylabel="Connectivity", colorbar_title="Density",
    aspect_ratio = 0.85, c = cgrad(:jet, scale = :linear), clims = (0, 1)
)
p3 = heatmap(grid_length,grid_length, z_matrix3,
    xlabel="Availability", ylabel="Connectivity", colorbar_title="Density",
    aspect_ratio = 0.85, c = cgrad(:jet, scale = :linear), clims = (0, 1)
)
p4 = heatmap(grid_length,grid_length, z_matrix4,
    xlabel="Availability", ylabel="Connectivity", colorbar_title="Density",
    aspect_ratio = 0.85, c = cgrad(:jet, scale = :linear), clims = (0, 1)
)
plot(p1, p2, p3, p4, size = (800, 800))
