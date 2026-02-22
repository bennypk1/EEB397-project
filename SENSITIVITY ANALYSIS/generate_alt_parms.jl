
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
# Order of parameters: (c_n, e_n, mu_n, alpha)

# define constriant objects
A <- rbind(
  c(-1,  0,  0,  0),   # check if all 4 are greater than 0
  c( 0, -1,  0,  0),
  c( 0,  0, -1,  0),
  c( 0,  0,  0, -1),
  c( 0,  0,  0,  1),   # alpha and c_n less than 1
  c( 1,  0,  0,  0),
  c(-1,  1,  1,  0),   # c_n >= e_n + mu_n
  c( 0, -1,  1,  0))   # mu_n <= e_n
b <- c(
  0,  # check if all 4 are greater than 0
  0,
  0,  
  0,  
  1,  # alpha and c_n less than 1
  SOFT_Cn_MAX,
  0,  # c_n >= e_n + mu_n
  0) # mu_n <= e_n
n <- 1000

# implement hit and run
constr <- list(
  constr = A,
  rhs    = b,
  dir    = rep("<=", 8)
)
state <- har.init(constr)
result <- har.run(state, n.samples=n)
samples <- result$samples

# convert to dataframe with readable names
candidateParmData <- as.data.frame(samples)
colnames(candidateParmData) <- c("c_n", "e_n", "alpha", "mu_n")
"""
@rget candidateParmData
candidateParmData = DataFrame(candidateParmData)

# Julia function to filter for nuisance parm combinations that also satisfy joint constraints with main parms
function filter_candidate_nuisance_params_by_validity!(candidateNuisanceParmData)
  candidateNuisanceParmData = filter(row -> begin
      is_parms_valid([row.c_n, row.e_n, row.alpha, row.mu_n])
    end, candidateNuisanceParmData)
  return candidateNuisanceParmData
end

# check diagnostics and write valid nuisance parm combinations to a .csv file with a custom name
# Note: outputted dataframe should be a full list of 16 parameter combinations
filename = "ParmSet1.csv"
filteredCandidates = filter_candidate_nuisance_params_by_validity!(candidateParmData)
augmentedCandidates = augment_filtered_nuisance_parms(filteredCandidates)
CSV.write(filename, augmentedCandidates)
println("Filtered nuisance parameters written to CSV: $filename")