
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
    return c1 && c2 && c3
end