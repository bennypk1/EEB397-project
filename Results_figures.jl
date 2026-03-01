using Revise
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))

# base parameters used in all plots
basePlotParams = [CANONICAL_PARMS; [0.5, 0.25]; [0.5, 0.5, 0]]
basePlotInit = CANONICAL_INIT
basePlotTimespan = [0, 100]
# basePlotGrain = 0.0025
basePlotGrain = 0.01
basePlotParams_altered = copy(basePlotParams)

functionalGrain = 0.05





# change alpha and beta HERE
basePlotParams_altered[14] = 0.5
basePlotParams_altered[15] = 0.5

# change e and gamma HERE
basePlotParams_altered[4] = 0.15
basePlotParams_altered[16] = 0

# run the plot
# p = LiaoTypeSpeciesRichnessMap(
#     # SimpleFoodChain!,
#     SimpleFoodChain_GlobalConsumers!,
#     # Omnivory!,
#     # ExploitativeCompetition!,
#     basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)
# FigureGradeLiaoTypePlot(p)

# run Figure 3.2b generator
# theDataE1 = ProportionalPersistenceData(4, 0.025, 0.5, SimpleFoodChain_GlobalConsumers!,
#     basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     functionalGrain, basePlotInit, basePlotTimespan)
# theDataGamma = ProportionalPersistenceData(16, 0.025, 0.5, SimpleFoodChain_GlobalConsumers!,
#     basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     functionalGrain, basePlotInit, basePlotTimespan)
# p = Figure3_2b(theDataGamma, theDataE1)

# run Figure 3.2b on problematic parms
problematicParms = [0.9480874493022413, 0.9480874493022413, 0.9480874493022413, 0.05, 0.2998348678681639, 0.2998348678681639, 0.2998348678681639, 0.24208056136531367, 0.24208056136531367, 0.24208056136531367, 4.0, 0.24916129448582738, 0.7508387055141726, 0.0]
theDataE1 = ProportionalPersistenceData(4, 0.025, 0.5, SimpleFoodChain!,
    problematicParms[1:11], problematicParms[12:14],
    functionalGrain, basePlotInit, basePlotTimespan)
theDataGamma = ProportionalPersistenceData(16, 0.025, 0.5, SimpleFoodChain!,
    problematicParms[1:11], problematicParms[12:14],
    functionalGrain, basePlotInit, basePlotTimespan)
p = Figure3_2b(theDataGamma, theDataE1)


# testing
# x = LiaoTypeGrid(Resource!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)

# x = LiaoTypeGridLEGACYY(Resource!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)

# x = ProportionalPersistencePoint(ExploitativeCompetition!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)