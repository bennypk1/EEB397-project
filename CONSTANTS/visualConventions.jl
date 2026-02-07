
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
# 1 = no species present
# 2 = only resource
# 3 = resource + 2
# 4 = resource + 3
# 5 = all 3
