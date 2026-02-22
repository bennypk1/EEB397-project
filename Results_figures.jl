using Revise
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))

# base parameters used in all plots
basePlotParams = [CANONICAL_PARAMS; [0.5, 0.25]; [0.5, 0.5, 0]]
basePlotParams[4] = 0.15 # this is for more dramatic plots in the Resource plots
basePlotInit = CANONICAL_INIT
basePlotTimespan = [0, 100]
# basePlotGrain = 0.0025
basePlotGrain = 0.015
basePlotParams_altered = copy(basePlotParams)





# change alpha and beta HERE
basePlotParams_altered[14] = 0.5
basePlotParams_altered[15] = 0.5

# change e and gamma HERE
basePlotParams_altered[4] = 0.05
basePlotParams_altered[16] = 0

# run the plot
p = LiaoTypeSpeciesRichnessMap(
    Omnivory!,
    basePlotParams_altered[1:11], basePlotParams_altered[14:16],
    basePlotGrain, basePlotInit, basePlotTimespan)
FigureGradeLiaoTypePlot(p)

# testing
# x = LiaoTypeGrid(Resource!, basePlotParams_altered[1:11], basePlotParams_altered[14:16],
#     basePlotGrain, basePlotInit, basePlotTimespan)