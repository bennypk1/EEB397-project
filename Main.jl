using Revise
using CSV
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))

################################################################################################################

# SIMPLE SINGLE RUN
parms_full = CANONICAL_FULL
init = CANONICAL_INIT
timespan = CANONICAL_TIMESPAN
# plotRun(Resource!, parms_full, init, [0, 100])

problematicParms = [0.871882111605091, 0.871882111605091, 0.871882111605091, 0.05, 0.36766583968624617, 0.36766583968624617, 0.36766583968624617, 0.1620643721349284, 0.1620643721349284, 0.1620643721349284, 4.0, 0.30862956543136105, 0.6913704345686389, 0.0]

# RUN COEXISTANCE PARTTERN PLOTS
# p = SpeciesPersistencePlot2(16, 0.025, 0.5, Omnivory!, CANONICAL_PARMS, [1.0, 0.0, 0.0])
# p = PatternPersistencePlot(4, 0.025, 1.0, Omnivory!, CANONICAL_PARMS, [1.0, 0.0, 0.0])
# prettyPatternPlot2!(p)
# LiaoTypeSpeciesRichnessMap(ExploitativeCompetition!,
#     Vector(parmSample[3,1:11]), Vector(parmSample[3,12:14]), 0.0025, init, timespan)
# LiaoTypeHeatMap6(Omnivory!,
#     Vector(parmSample[3, 1:11]), Vector(parmSample[3, 12:14]), 0.01, init, timespan; focalValue=2.0)
LiaoTypeHeatMap7(Omnivory!,
    CANONICAL_PARMS, [1.0, 0, 0], 0.01, init, timespan; focalValue=2.0)