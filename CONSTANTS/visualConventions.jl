
using Revise

# plotRun
LABELS_plotRun = ["P₁", "P₂", "P₍₂₃₎", "P₍₁₃₎", "p₁₁", "pᵤ₁"]
LINEWIDTHS_plotRun = [4, 4, 4, 4, 2, 2]
LINEOPACITY_plotRun = [1.0, 1.0, 1.0, 1.0, 0.75, 0.75]


# Liao-type Persistence Thresholds
PERSISTENCE_THRESHOLD = 0.05

PERSISTENCE_CODE = [0.4, 1.0]

SPECIES_DISTRIBUTION_IDS = [0, 0.25, 0.5, 0.75, 1]
# length = 5 = number of possible species distributions in the model
# our convention is (by index, the actual values represent colors):

# 1 = no species present (purple)
# 2 = only resource      (blue)
# 3 = resource + 2       (turquoise)
# 4 = resource + 3       (green)     
# 5 = all species        (yellow)

COEXISTANCE_PATTERN_COLORS = ["#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"]


# Proportional Persistance-related plots

PP_GRAIN = 0.1
CANDIADATE_PP_PARAMETER_INDICES = [4, 8, 16] # can extend functionality to incorporate more in the future if necessary