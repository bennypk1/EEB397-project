
using Revise
using RCall
using DataFrames
using CSV
include(joinpath(@__DIR__, "sensitivity_helpers.jl"))
include(joinpath(@__DIR__, "sensitivity_constants.jl"))

# This script uses R's hitandrun package to uniformly sample a large set of biologically realistic parameter
# constraints. Once generated and validated, the parameter combinations are written to a .csv file for later
# use. Therefore, unless constraint changes or additional parameters are needed, ONLY RUN THIS SCRIPT ONCE.

# R code to define nuisance parameter constraints and generate valid combinatons
@rput SOFT_Cn_MAX
R"""
library(hitandrun)
set.seed(416)
# Order of parameters: (c_n, e_n, mu_n, alpha,      e_1, gamma)

# define constriant objects
A <- rbind(
  c(-1,  0,  0,  0,      0,  0),   # check if all 6 are greater than 0
  c( 0, -1,  0,  0,      0,  0),
  c( 0,  0, -1,  0,      0,  0),
  c( 0,  0,  0, -1,      0,  0),
  c( 0,  0,  0,  0,     -1,  0),
  c( 0,  0,  0,  0,      0, -1),
  c( 0,  0,  0,  1,      0,  0),   # alpha and c_n less than 1
  c( 1,  0,  0,  0,      0,  0),
  c(-1,  1,  1,  0,      0,  0),   # c_n >= e_n + mu_n
  c(-1,  0,  0,  0,      2,  1),   # c_n >= 2e_1 + gamma (2 to avoid e_1 getting too large, because otherwise uninteresting)
  c( 0, -1,  1,  0,      0,  0))   # mu_n <= e_n
b <- c(
  0,  # check if all 6 are greater than 0
  -0.025, # need a specialized constraint for e_n because e_n too close to 0 destroys the ODE
  0,  
  0, 
  -0.025, # need a specialized constraint for e_1 because e_1 too close to 0 destroys the ODE
  0, 
  1,  # alpha and c_n less than 1
  SOFT_Cn_MAX,
  0,  # c_n >= e_n + mu_n
  0,  #  c_n >= e_1 + gamma
  0) # mu_n <= e_n
n <- 1000

# implement hit and run
constr <- list(
  constr = A,
  rhs    = b,
  dir    = rep("<=", length(b))
)
state <- har.init(constr)
result <- har.run(state, n.samples=n)
samples <- result$samples

# convert to dataframe with readable names
candidateParmData <- as.data.frame(samples)
colnames(candidateParmData) <- c("c_n", "e_n", "mu_n", "alpha",      "e_1", "gamma")
"""
@rget candidateParmData
candidateParmData = DataFrame(candidateParmData)

# LEGACY
# # Julia function to filter for nuisance parm combinations that also satisfy joint constraints with main parms
# function filter_candidate_nuisance_params_by_validity!(candidateNuisanceParmData)
#   candidateNuisanceParmData = filter(row -> begin
#       is_parms_valid([row.c_n, row.e_n, row.mu_n, row.alpha])
#     end, candidateNuisanceParmData)
#   return candidateNuisanceParmData
# end

# TODO: NEED TO REWRITE GENERATED PARAMETER VALIDATION CODE HERE


# # check diagnostics and write valid nuisance parm combinations to a .csv file with a custom name
# # Note: outputted dataframe should be a full list of 16 parameter combinations
filename = "ParmSet2.csv"
# filteredCandidates = filter_candidate_nuisance_params_by_validity!(candidateParmData)
augmentedParms = augment_compact_parms(candidateParmData)
CSV.write(filename, augmentedParms)
println("Parameter set written to CSV: $filename")