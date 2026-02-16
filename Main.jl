using Revise
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))

################################################################################################################

params = [CANONICAL_PARAMS; [0.5, 0.25]; [0, 1, 0]]
init = CANONICAL_INIT
init[1] = 0.96
timespan = CANONICAL_TIMESPAN
# plotRun(Resource!, params, init, [0, 100])

grain = 0.01
CANONICAL_PARAMS[4] = 0.05 # set intrinsic extinction of resource
# CANONICAL_PARAMS[7] = CANONICAL_PARAMS[6] * 3 # apply penalty to species 3 when feeding on 1 instead of 2 (TO GET EXTRA COOL BEHAVIOIR IN OMNIVORY MODEL)
CANONICAL_PARAMS[8] = 0.025 # set additive mortality rate due to 2 feeding on 1
CANONICAL_PARAMS[10] = 0 # set additive mortality rate due to 3 feeding on 1

# WeightedProportionalPersistencePlot3(16, Omnivory!, CANONICAL_PARAMS, [1.0, 0, 0.0], grain, init, timespan)
# x = ProportionalPersistencePoint8(ExploitativeCompetition!, CANONICAL_PARAMS, [0, 1, 0.025], grain, init, timespan)
# ProportionalPersistencePlot(8, ExploitativeCompetition!, CANONICAL_PARAMS, [0., 1.0, 0.0], grain, init, timespan)
# WeightedProportionalPersistencePlot(8, ExploitativeCompetition!, CANONICAL_PARAMS, [0., 1.0, 0.0], grain, init, timespan)
# LiaoTypeSpeciesRichnessMap(ExploitativeCompetition!, CANONICAL_PARAMS, [1, 0, 0], grain, init, timespan)
# testGrid = LiaoTypeGrid(ExploitativeCompetition!, CANONICAL_PARAMS, [1, 0, 0], grain, init, timespan)
