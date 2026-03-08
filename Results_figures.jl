using Revise
using CSV
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
basePlotParams_altered[4] = 0.025
basePlotParams_altered[16] = 0

# LANDSCAPE DATA ANALYSIS
focalLandscapeData = AllLandscapeData(4, 0.025, 0.5, ExploitativeCompetition!, basePlotParams_altered[1:11], basePlotParams_altered[14:16])
CSV.write("landscapeDataE1EC.csv", focalLandscapeData)



# run Figure 3.2b generator
# theDataE1 = ProportionalPersistenceData(4, 0.025, 0.5, SimpleFoodChain_GlobalConsumers!,
#     basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     functionalGrain, basePlotInit, basePlotTimespan)
# theDataGamma = ProportionalPersistenceData(16, 0.025, 0.5, SimpleFoodChain_GlobalConsumers!,
#     basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     functionalGrain, basePlotInit, basePlotTimespan)
# p = Figure3_2b(theDataGamma, theDataE1)

# run Figure 3.2b on problematic parms
# problematicParms = [0.9480874493022413, 0.9480874493022413, 0.9480874493022413, 0.05, 0.2998348678681639, 0.2998348678681639, 0.2998348678681639, 0.24208056136531367, 0.24208056136531367, 0.24208056136531367, 4.0, 0.24916129448582738, 0.7508387055141726, 0.0]
# theDataE1 = ProportionalPersistenceData(4, 0.025, 0.5, SimpleFoodChain!,
#     problematicParms[1:11], problematicParms[12:14],
#     functionalGrain, basePlotInit, basePlotTimespan)
# theDataGamma = ProportionalPersistenceData(16, 0.025, 0.5, SimpleFoodChain!,
#     problematicParms[1:11], problematicParms[12:14],
#     functionalGrain, basePlotInit, basePlotTimespan)
# p = Figure3_2b(theDataGamma, theDataE1)


# testing
# x = LiaoTypeGrid(Resource!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)
# x = LiaoTypeHeatMap(Resource!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)



# p = LiaoTypeSpeciesRichnessMap(ExploitativeCompetition!,
#     basePlotParams_altered[1:11], basePlotParams_altered[14:16], 0.02, init, timespan)
# FigureGradeLiaoTypePlot(p)

# available point functions:
# ProportionalPersistencePoint
# WeightedByResourcePersistencePoint
# LinkDensityPersistencePoint

# p = PatternPersistencePlot2(16, 0.025, 1, ExploitativeCompetition!, basePlotParams_altered[1:11], basePlotParams_altered[14:16];
#     pointFunction=LinkDensityPersistencePoint3)
# prettyPatternPlot!(p)

# p = ConsumerPersistencePlot2(16, 0.025, 1, ExploitativeCompetition!, basePlotParams_altered[1:11], basePlotParams_altered[14:16] ;
# pointFunction = WeightedByResourcePersistencePoint)
# prettyPatternPlot!(p)

# ppData = ProportionalPersistenceData(16, 0.025, 1, ExploitativeCompetition!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     0.025, CANONICAL_INIT, CANONICAL_TIMESPAN; pointFunction=LinkDensityPersistencePoint3)