
# Create a unit grid with a set grain
function create_unit_grid(grain)
    result = []
    s_possibilities = 0:grain:1
    q_possibilities = 0:grain:1
    for s_point in s_possibilities
        for q_point in q_possibilities
            push!(result, [s_point, q_point])
        end
    end
    return result
end

# Check if a grid point is valid
function is_valid(grid_point)
    c1 = grid_point[2] > (2 - 1/grid_point[1])
    c2 = grid_point[1] > 0.025
    c3 = grid_point[2] > 0.025
    c4 = grid_point[1] < 1
    c5 = grid_point[2] < 1
    return c1 && c2 && c3 && c4 && c5
end

# Convert [s, qₛₛ] to [pᵤ, pᵤᵤ]
function convert_landscape_params(landscape_params)
    this_pᵤ = 1 - landscape_params[1]
    this_pᵤᵤ = this_pᵤ - landscape_params[1] * (1 - landscape_params[2])
    return([this_pᵤ, this_pᵤᵤ])
end

# performance function for plot_ODE5
function get_init_from_trophic_config(initial_densities, trophic_configuration)
    P1, P12, P23, P13, P11, Pu1 = initial_densities
    I23, I13 = trophic_configuration
    P23_init = I23 == 1 ? P23 : 0.0
    P13_init = I13 == 1 ? P13 : 0.0
    return [P1, P12, P23_init, P13_init, P11, Pu1]
end