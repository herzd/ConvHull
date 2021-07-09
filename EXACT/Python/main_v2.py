"""
Computes the convex hull for production envelopes of metabolic network. Solution is
the list of hyperplanes and set of extreme points of the Convex hull. Inputs are:
  - REACTION_FILE : absolute path to file containing one reaction name string in each line, the names must be aligned with the columns of the Stoichiometric matrix.
  - STOICHIOMETRIC_FILE : absolute path to file containing Stoichiometric matrix
  - DOMAIN_FILE : absolute path to file containing lower and upper bounds for each reaction, one pair in a line
  - INPUT_REACTIONS: list of indices for the dimensions onto which the CH should be computed
"""

import chm_exact
REACTION_FILE = "/home/daniel/ConvHull/DATA/toy/toy_reactions.txt"
STOICHIOMETRIC_FILE = "/home/daniel/ConvHull/DATA/toy/toy_stoichs.txt"
DOMAIN_FILE = "/home/daniel/ConvHull/DATA/toy/toy_domains.txt"
INPUT_REACTIONS = [0, 1]

chm_exact.compute_CH(REACTION_FILE, STOICHIOMETRIC_FILE, DOMAIN_FILE, INPUT_REACTIONS)
