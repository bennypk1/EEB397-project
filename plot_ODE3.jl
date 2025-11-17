using DifferentialEquations
using Plots
using DataFrames
include("helper_functions.jl")

# Replace function body with one of the ODEs
function testODE!(dP, P, params, t)
    c₁, c₂, e₁, e₂, μ, pᵤ, pᵤᵤ, z = params
    P₁, P₂, P₁₁, Pᵤ₁ = P
    dP[1] = c₁ * (P₁ - P₁₁ - Pᵤ₁) - e₁ * P₁ - μ * P₂
    dP[2] = c₂ * (P₁ - P₂) * P₂ - (e₁ + e₂ + μ) * P₂
    dP[3] = c₁ * (P₁ - P₁₁ - Pᵤ₁) * (1/z + (z-1)/z * ((P₁ - P₁₁ - Pᵤ₁) / (1 - pᵤ - P₁))) - P₁₁ * (e₁ + μ * P₂)
    dP[4] = c₁ * ((z-1)/z) * (pᵤ - Pᵤ₁ - pᵤᵤ) * ((P₁ - P₁₁ - Pᵤ₁) / (1 - pᵤ - P₁)) - Pᵤ₁ * (e₁ + μ * P₂)
end

# Single solution over time
this_s = 0.59
this_q = 0.5
conv_params = convert_landscape_params([this_s, this_q])
fixed_params = [1, 1, 0.05, 0.05, 0.0025, conv_params[1], conv_params[2], 4]
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
        curr_params = [1, 1, 0.05, 0.05, 0.0025, conv_point[1], conv_point[2], 4]
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
