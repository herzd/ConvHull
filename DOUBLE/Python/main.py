
from chm_double import CHM
import time

reactions = [0, 2, 71] 

chm = CHM(reactions)

chm.set_stoichiometric_matrix("../../DATA/Aeq.txt")
chm.set_reaction_domains("../../DATA/domains.txt")

# It is slightly faster to store the stoichiometric constraints
# as a gruboi model, rather than reinitializing each time
chm.set_model()

start_time = time.time()
print("Initial Points...")
init_points = chm.initial_points()

print("Initial Hull...")
hull = chm.initial_hull(init_points)

print("Incremental Refinement...")
(final_points, final_hyperplanes) = chm.incremental_refinement(hull, init_points)

print("Final Extreme Points...")
print(final_points)
#print(final_hyperplanes)

print("--- {} seconds ---".format(time.time() - start_time))
print("--- " + str(chm.solver_count) + " models were build. ---")
