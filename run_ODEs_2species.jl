using DifferentialEquations
using Plots
using DataFrames
include("ODEs_2species.jl")

# Helper functions
function is_valid(grid_point)
    return grid_point[2] > (2 - 1/grid_point[1])
end

# RUN WITH DIFFERENT LANDSCAPE PARAMETERS
# 1 Create Parameter Grid
grain = 0.01
s_possibilities = 0:grain:1
q_possibilities = 0:grain:1
grid = []
for s_point in s_possibilities
    for q_point in q_possibilities
        push!(grid, [s_point, q_point])
    end
end

# 2 Initialize Simulation
fixed_params = [1, 1, 0.05, 0.05, 0.025]
time_span = [0, 10000]
initial_density = [0.01, 0.05, 0.08, 0.0005]
n_sims = length(grid)
data = DataFrame(Availability = Float64[], Connectivity = Float64[], PreyDensity = Float64[], PredDensity = Float64[])

# Simulation
for point in grid
    if !is_valid(point) # check if point is in invalid region
        new_row = (point[1], point[2], NaN, NaN) # TODO: replace 1 with something else that changes color
        push!(data, new_row)
    else
    curr_point = [point[1], point[2] * point[1]]
    curr_params = [fixed_params; curr_point] # since rho = qss * s
    # curr_prob = ODEProblem(predatorPrey2AllGlobal!, initial_density, time_span, curr_params)
    curr_prob = ODEProblem(predatorPrey2LocalPrey4!, initial_density, time_span, curr_params)
    curr_sol = solve(curr_prob)
    sol_end = curr_sol[end]
    new_row = (point[1], point[2], sol_end[1], sol_end[2])
    push!(data, new_row)
    end
end

# Visualization
n_avail = length(unique(data.Availability))        # number of unique Avail values
n_pair = length(unique(data.Connectivity))      # number of unique PairProb values

z_matrix = reshape(data.PreyDensity, n_pair, n_avail)
heatmap(
    unique(data.Availability),           # x-axis
    unique(data.Connectivity),        # y-axis
    z_matrix,
    xlabel="Availability",
    ylabel="Connectivity",
    colorbar_title="Density (Prey)",
    aspect_ratio = 1,
    c = cgrad(:jet, scale = :linear),
)
