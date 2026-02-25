using Revise
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))

# base parameters used in all plots
basePlotParams = [CANONICAL_PARAMS; [0.5, 0.25]; [0.5, 0.5, 0]]
basePlotInit = CANONICAL_INIT
basePlotTimespan = [0, 100]
basePlotGrain = 0.0025
# basePlotGrain = 0.01
basePlotParams_altered = copy(basePlotParams)





# change alpha and beta HERE
basePlotParams_altered[14] = 0.2
basePlotParams_altered[15] = 0.8

# change e and gamma HERE
basePlotParams_altered[4] = 0.15
basePlotParams_altered[16] = 0

# run the plot
p = LiaoTypeSpeciesRichnessMap(
    # SimpleFoodChain!,
    # SimpleFoodChain_GlobalConsumers!,
    # Omnivory!,
    ExploitativeCompetition!,
    basePlotParams_altered[1:11], basePlotParams_altered[14:16],
    basePlotGrain, basePlotInit, basePlotTimespan)
FigureGradeLiaoTypePlot(p)

# run Figure 3.2b generator
# theDataE1 = ProportionalPersistenceData(4, 0.025, 0.5, SimpleFoodChain_GlobalConsumers!,
#     basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)
# theDataGamma = ProportionalPersistenceData(16, 0.025, 0.5, SimpleFoodChain_GlobalConsumers!,
#     basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)
# p = Figure3_2b(theDataE1, theDataGamma)


# testing
# x = LiaoTypeGrid(Resource!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)

# x = LiaoTypeGridLEGACYY(Resource!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)

# x = ProportionalPersistencePoint(ExploitativeCompetition!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)