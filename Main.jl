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

# problematicParms = [0.8390240473409691, 0.8390240473409691, 0.8390240473409691, 0.03682354898685283, 0.4840526218684192, 0.4840526218684192, 0.4840526218684192, 0.2973819645287637, 0.2973819645287637, 0.2973819645287637, 4.0, 0.5658720172120596, 0.4341279827879404, 0.22771263432574465]
# p = LiaoTypeSpeciesRichnessMap(Omnivory!,
#     problematicParms[1:11], problematicParms[12:14], 0.01, init, timespan)
# FigureGradeLiaoTypePlot(p)







# RUN COEXISTANCE PARTTERN PLOTS
# p = SpeciesPersistencePlot2(16, 0.025, 0.5, Omnivory!, CANONICAL_PARMS, [1.0, 0.0, 0.0])
# p = PatternPersistencePlot(4, 0.025, 1.0, Omnivory!, CANONICAL_PARMS, [1.0, 0.0, 0.0])
# prettyPatternPlot2!(p)
# LiaoTypeSpeciesRichnessMap(ExploitativeCompetition!,
#     Vector(parmSample[3,1:11]), Vector(parmSample[3,12:14]), 0.0025, init, timespan)
# LiaoTypeHeatMap6(Omnivory!,
#     Vector(parmSample[3, 1:11]), Vector(parmSample[3, 12:14]), 0.01, init, timespan; focalValue=2.0)
# LiaoTypeHeatMap7(Omnivory!,
#     CANONICAL_PARMS, [1.0, 0, 0], 0.01, init, timespan; focalValue=2.0)