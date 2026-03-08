using Revise
using DataFrames
using CSV
using Random
include(joinpath(@__DIR__, "plot_ODEFinal.jl"))
include(joinpath(@__DIR__, "helper_functions.jl"))
include(joinpath(@__DIR__, "CONSTANTS", "modelsInfo.jl"))


# This script contains all target functions used for sensitivity analysis. Each target function verifies a
# specific, qualitative result, for a single, arbitrary (biologically valid). As such, each function must ONLY
# output true or false (functions may also print diagnostic information). The final function in this script runs
# any function of choice against a pre-defined dataframe of biologically valid parameter combinations.

# Notes:
# 1) currently not checking againts different input parameters, timespans, grid extents or persistence thresholds

# Read in OG parm data
# ParmSet1 = CSV.read("ParmSet1.csv", DataFrame)
# parmSample = ParmSet1[sample(1:nrow(ParmSet1), 300, replace=false), :]
# parmExample = Vector(parmSample[1, :])







# set seed for sampling
Random.seed!(416)

# Read in better parm data
ParmSet2 = CSV.read("ParmSet2.csv", DataFrame)
parmSample = ParmSet2[randperm(nrow(ParmSet2))[1:50], :]
parmExample = Vector(parmSample[1, :])

###############################################################################################################
# TARGET FUNCTIONS
###############################################################################################################

# Some things to check
# 1) in Simple Food Chain, is species 2 persistence ever affected by species 3


function TargetResource_dispersalpersistencetradeoff()
end

# VERIFIED
function TargetCommunity_resourcepersistencenotafectedbyconsumerpresence(parmRow; tol=0.01)
    # input check
    if length(parmRow) != 14
        error("{TargetCommunity_resourceper...}: input parameter list does not have a length of exactly 14.")
    end
    # respurce-only
    dataR = ProportionalPersistencePoint(Resource!, parmRow[1:11], parmRow[12:14],
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
    # 4 distinct food web motifs checked
    dataSFC = ProportionalPersistencePoint(SimpleFoodChain!, parmRow[1:11], parmRow[12:14],
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
    dataSFCGC = ProportionalPersistencePoint(SimpleFoodChain_GlobalConsumers!, parmRow[1:11], parmRow[12:14],
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
    dataEC = ProportionalPersistencePoint(ExploitativeCompetition!, parmRow[1:11], parmRow[12:14],
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
    dataO = ProportionalPersistencePoint(Omnivory!, parmRow[1:11], parmRow[12:14],
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
    # verify target claim
    if !equal_within_tolerance([dataR[1], dataSFC[1], dataSFCGC[1], dataEC[1], dataO[1]]; tol)
        println([dataR[1], dataSFC[1], dataSFCGC[1], dataEC[1], dataO[1]])
        return false
    end
    return true
end

# VERIFIED
function TargetFoodChains_species2persistencenotaffectedbyspecies3presence(parmRow; tol=0.01)
    # input check
    if length(parmRow) != 14
        error("{TargetFoodChains_species2per...}: input parameter list does not have a length of exactly 14.")
    end
    # species 3 absent
    control_init = copy(CANONICAL_INIT)
    control_init[3] = 0.0
    control_init[4] = 0.0 # just for good measure
    dataSFC_control = ProportionalPersistencePoint(SimpleFoodChain!, parmRow[1:11], parmRow[12:14],
        0.05, control_init, CANONICAL_TIMESPAN)
    dataSFCGC_control = ProportionalPersistencePoint(SimpleFoodChain_GlobalConsumers!, parmRow[1:11], parmRow[12:14],
        0.05, control_init, CANONICAL_TIMESPAN)
    # species 3 present
    dataSFC = ProportionalPersistencePoint(SimpleFoodChain!, parmRow[1:11], parmRow[12:14],
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
    dataSFCGC = ProportionalPersistencePoint(SimpleFoodChain_GlobalConsumers!, parmRow[1:11], parmRow[12:14],
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
    # verify target claim
    cSFC = dataSFC_control[3] + dataSFC_control[5] - (dataSFC[3] + dataSFC[5]) < tol
    cSFCGC = dataSFCGC_control[3] + dataSFCGC_control[5] - (dataSFCGC[3] + dataSFCGC[5]) < tol
    if !(cSFC && cSFCGC)
        println([dataSFC_control[3] + dataSFC_control[5], dataSFC[3] + dataSFC[5]])
        println("SFCGC:::", [dataSFCGC_control[3] + dataSFCGC_control[5], dataSFCGC[3] + dataSFCGC[5]])
        return false
    end
    return true
end

# VERIFIED
function TargetResource_resourcepersistanceaffectedmorebyethangamma(parmRow)
    # input check
    if length(parmRow) != 14
        error("{TargetResource_resourceper...}: input parameter list does not have a length of exactly 14.")
    end
    # resource-only model: changing E1
    dataR_E1 = ProportionalPersistenceData(
        4, 0.025, 0.5,                             # specify E1 + searching range
        Resource!, parmRow[1:11], parmRow[12:14], # specify model + parameters
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN  # extra stuff
    )
    # resource-only model: changing Gamma
    dataR_Gamma = ProportionalPersistenceData(
        16, 0.025, 0.5,                             # specify Gamma + searching range
        Resource!, parmRow[1:11], parmRow[12:14],  # specify model + parameters
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN   # extra stuff
    )
    # output check 1: consumers shouldn't be present
    if !all(dataR_Gamma.PP_5 + dataR_Gamma.PP_4 + dataR_Gamma.PP_3 .== 0)
        error("{TargetResource_resourceper...}: some outputted PP data for Resource-only model had non-zero Consumer PPs.")
    end
    if !all(dataR_E1.PP_5 + dataR_E1.PP_4 + dataR_E1.PP_3 .== 0)
        error("{TargetResource_resourceper...}: some outputted PP data for Resource-only model had non-zero Consumer PPs.")
    end
    # ouput check 2: resource persistence should be non-increasing
    if any(diff(dataR_Gamma.PP_2) .> 0)
        error("{TargetResource_resourceper...}: some outputted PP data showed that persistance increased with increasing gamma.")
    end
    if any(diff(dataR_E1.PP_2) .> 0)
        error("{TargetResource_resourceper...}: some outputted PP data showed that persistance increased with increasing e₁.")
    end
    # verify claim
    slope_E1 = mean(diff(dataR_E1.PP_2))
    slope_Gamma = mean(diff(dataR_Gamma.PP_2))
    if !(slope_E1 < slope_Gamma) # Note: checking that E1 slope is "more negative"
        println([slope_E1, slope_Gamma])
        return false
    end
    return true
end

# VERIFIED
function TargetCommunity_consumerpersistenceaffectedmorebygammathanresource(parmRow)
    # input check
    if length(parmRow) != 14
        error("{TargetCommunity_consumerper...}: input parameter list does not have a length of exactly 14.")
    end
    # generate data
    dataSFC = ProportionalPersistenceData(
        16, 0.025, 0.5,                             # specify Gamma + searching range
        SimpleFoodChain!, parmRow[1:11], parmRow[12:14],  # specify model + parameters
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN   # extra stuff
    )
    dataSFCGC = ProportionalPersistenceData(
        16, 0.025, 0.5,                             # specify Gamma + searching range
        SimpleFoodChain_GlobalConsumers!, parmRow[1:11], parmRow[12:14],  # specify model + parameters
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN   # extra stuff
    )
    dataEC = ProportionalPersistenceData(
        16, 0.025, 0.5,                             # specify Gamma + searching range
        ExploitativeCompetition!, parmRow[1:11], parmRow[12:14],  # specify model + parameters
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN   # extra stuff
    )
    dataO = ProportionalPersistenceData(
        16, 0.025, 0.5,                             # specify Gamma + searching range
        Omnivory!, parmRow[1:11], parmRow[12:14],  # specify model + parameters
        0.05, CANONICAL_INIT, CANONICAL_TIMESPAN   # extra stuff
    )
    # get thresholds and checks
    thresholdResourceSFC = 1 .- dataSFC.PP_1
    thresholdConsumerSFC = dataSFC.PP_3 + dataSFC.PP_4 + dataSFC.PP_5
    cSFC = no_meaningful_slope_violation(thresholdResourceSFC, thresholdConsumerSFC)
    thresholdResourceSFCGC = 1 .- dataSFCGC.PP_1
    thresholdConsumerSFCGC = dataSFCGC.PP_3 + dataSFCGC.PP_4 + dataSFCGC.PP_5
    cSFCGC = no_meaningful_slope_violation(thresholdResourceSFCGC, thresholdConsumerSFCGC)
    thresholdResourceEC = 1 .- dataEC.PP_1
    thresholdConsumerEC = dataEC.PP_3 + dataEC.PP_4 + dataEC.PP_5
    cEC = no_meaningful_slope_violation(thresholdResourceEC, thresholdConsumerEC)
    thresholdResourceO = 1 .- dataO.PP_1
    thresholdConsumerO = dataO.PP_3 + dataO.PP_4 + dataO.PP_5
    cO = no_meaningful_slope_violation(thresholdResourceO, thresholdConsumerO)
    if cSFC && cSFCGC && cEC && cO
        return true
    else
        println(parmRow)
        return false
    end
end


# function TargetCommunity_regionalspecies2persistenceinSFCandECandOissimilar4(parmRow)
#     # input check
#     if length(parmRow) != 14
#         error("{TargetCommunity_regionalspe...}: input parameter list does not have a length of exactly 14.")
#     end
#     # get data
#     dataSFC = ProportionalPersistencePoint(SimpleFoodChain!, parmRow[1:11], parmRow[12:14],
#         0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
#     dataEC = ProportionalPersistencePoint(ExploitativeCompetition!, parmRow[1:11], parmRow[12:14],
#         0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
#     dataO = ProportionalPersistencePoint(Omnivory!, parmRow[1:11], parmRow[12:14],
#         0.05, CANONICAL_INIT, CANONICAL_TIMESPAN)
#     # verify claim
#     println([dataSFC[3] + dataSFC[5], dataEC[3] + dataEC[5], dataO[3] + dataO[5]])
#     epsilon = 0.02 # have dome small error tolerance here
#     return abs(dataSFC[3] + dataSFC[5] - (dataEC[3] + dataEC[5])) < epsilon
# end

###############################################################################################################
# RESULTS THAT DON'T REQUIRE SENSITIVITY ANALYSES (usually because can be deduced logically)
###############################################################################################################

# 1) Consumers always persist less than resources


###############################################################################################################
# TESTING FUNCTIONS
###############################################################################################################

# <targetFunction> must output a boolean value ; <parmSet> is a 14-column dataframe (no U, F)
function BinaryRobustnessAnalysis(targetFunction, parmSet)
    # input checks
    if nrow(parmSet) == 0
        error("{BinaryRobustnessAnalysis}: inputted parameter set is empty.")
    end
    if ncol(parmSet) != 14
        error("{BinaryRobustnessAnalysis}: inputted parameter set does not have exactly 14 columns.")
    end
    # check each row in parm set
    for row in eachrow(parmSet)
        if !targetFunction(Vector(row))
            return false
        end
    end
    return true
end