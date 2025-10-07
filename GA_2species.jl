using DifferentialEquations
using Plots
include("ODEs_2species.jl")

# 1) VISUALIZE SOLUTION DIRECTLY

# Note: ρₛₛ must be less than or equal to s, so that qₛₛ is a valid probability measure
parameters = (2, 2, 1, 1, 0.5, 1, 0.5, 0.5)
time_span = (0, 10)
initial_2SAllGlobal = [0.5, 0.5, 0.5, 0.5]
initial_2SLocalPrey4 = [0.5, 0.5, 0.5, 0.5]

prob_2SAllGlobal = ODEProblem(predatorPrey2AllGlobal!, initial_2SAllGlobal, time_span, parameters)
prob_2SLocalPrey4 = ODEProblem(predatorPrey2LocalPrey4!, initial_2SLocalPrey4, time_span, parameters)

sol_2SAllGlobal = solve(prob_2SAllGlobal)
sol_2SLocalPrey4 = solve(prob_2SLocalPrey4)

plot(sol_2SAllGlobal)
# plot(sol_2SLocalPrey4)



