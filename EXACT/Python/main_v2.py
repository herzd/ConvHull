"""
Computes the convex hull for production envelopes of metabolic network. Solution is
the list of hyperplanes and set of extreme points of the Convex hull. Inputs are:
  - REACTION_FILE : absolute path to file containing one reaction name string in each line, the names must be aligned with the columns of the Stoichiometric matrix.
  - STOICHIOMETRIC_FILE : absolute path to file containing Stoichiometric matrix
  - DOMAIN_FILE : absolute path to file containing lower and upper bounds for each reaction, one pair in a line
  - INPUT_REACTIONS: list of indices for the dimensions onto which the CH should be computed
"""

import chm_exact_v2
# given information
REACTION_FILE = "/home/dherzig/ConvHull/DATA/toy/toy_reactions.txt"
STOICHIOMETRIC_FILE = "/home/dherzig/ConvHull/DATA/toy/toy_stoichs.txt"
DOMAIN_FILE = "/home/dherzig/ConvHull/DATA/toy/toy_domains.txt"
INPUT_REACTIONS = [0, 1]

# create dictionary from the above files for further processing
PROBLEM_READ = chm_exact_v2.read_problem(REACTION_FILE, STOICHIOMETRIC_FILE, DOMAIN_FILE)
# create some information for qsopt_ex (this makes a list of 3 zeros)
OBJ = [0] * PROBLEM_READ["Aeq"].shape[1]
# this sets the first zero to 1
OBJ[INPUT_REACTIONS[0]] = 1
# create linear problem from the dictionary to get rid off global statements within functions
PROBLEM_CREATED = chm_exact_v2.create_lp(PROBLEM_READ,OBJ)
# extract reaction ids, to get rid off the global statements within functions
REACTION_IDS = PROBLEM_READ["rids"]
CHULL, EPTS = chm_exact_v2.compute_CH(PROBLEM_READ, INPUT_REACTIONS, REACTION_IDS, PROBLEM_CREATED)
# refine results
REFINED_CHULL, REFINED_EPTS = chm_exact_v2.incremental_refinement(CHULL, EPTS, INPUT_REACTIONS, REACTION_IDS, PROBLEM_CREATED)
print("refined convex hull after refinement:")
print(REFINED_CHULL)
print("refined set of points:")
print(REFINED_EPTS)
