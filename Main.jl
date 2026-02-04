using Revise
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "ODEFinal.jl"))




################################################################################################################

params = [CANONICAL_PARAMS; [0.1, 0.01]; [0, 1, 2]]
init = CANONICAL_INIT
timespan = CANONICAL_TIMESPAN
plotRun(Resource!, params, init, [0, 100])

# grain = 0.05
# y = LiaoTypeGrid1(goofLeg!, CANONICAL_PARAMS, [0, 1, 0], grain, init, timespan)
# LiaoTypeHeatMap(goofLeg!, CANONICAL_PARAMS, [1, 0, 0], grain, init, timespan)