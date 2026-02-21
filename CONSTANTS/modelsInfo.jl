
using Revise
using DataFrames

# MODEL PARAMETERS

# "safe" parameter combination inspired by Liao Supplimentary Material. Defines c, e, μ, and z
CANONICAL_PARAMS = [1, 1, 1, 0.05, 0.05, 0.05, 0.05, 0.0025, 0.0025, 0.0025, 4]

# starting densities used by Liao
CANONICAL_INIT = [0.01, 0.005, 0.0025, 0.0025, 0.002, 0.005]

# specifying different landscape types
LANDSCAPE_FULL = [0, 0]
LANDSCAPE_RANDOM = [0.32, 0.11]
LANDSCAPE_CONNECTED = [0.5, 0]

# specifying different resource life histories: cominations of α, β, γ
RD_GLOBAL = [0, 1, 0]
RD_LOCAL = [1, 0, 0]

# OTHER SPECIFICATIONS

CANONICAL_TIMESPAN = [0, 100]



# Sensitivity Analysis Stuff

# cube of main parameter values my analyses will be conducted on 
MAIN_CUBE = DataFrame()


# KEY SOFT CONSTRASINTS
SIGNIFICANT_P = 0.01