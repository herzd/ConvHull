"""
Computes the convex hull for production envelopes of metabolic network. Solution is
the list of hyperplanes and set of extreme points of the Convex hull. Inputs are:
* fname: name of file without extension (must be the same for all files
  - fname_r.txt: list of reaction names - order must follow that of S columns
  - fname_S.txt: Stoichiometric matrix
  - fname_d.txt : lb ub for each reaction
* impt_reactions: list of indices for the dimensions onto which the CH should be computed
"""

import chm_exact
REACTION_FILE = "/home/daniel/ConvHull/DATA/toy/toy_reactions.txt"
STOICHIOMETRIC_FILE = "/home/daniel/ConvHull/DATA/toy/toy_stoichs.txt"
DOMAIN_FILE = "/home/daniel/ConvHull/DATA/toy/toy_domains.txt"
REACTIONS = [0, 1]

chm_exact.compute_CH(REACTION_FILE, STOICHIOMETRIC_FILE, DOMAIN_FILE, REACTIONS)
