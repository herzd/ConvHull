
import gurobipy as gp
import numpy as np

from itertools import combinations
from sympy import Matrix


class CHM:

    def __init__(self, impt_reactions, tol=pow(10, -5), dec_cmp=5, dec_pnt=9, adjust_tol=True, bug_report=False):
        """
        Keyword arguments:
        impt_reactions -- list of reactions of interest
        tol -- acceptable tolerance between points (default 10^-5)
        dec_cmp -- number of decimal places points are compared to (default 5)
        dec_pnt -- number of decimal places points are rounded to (default 9)
        adjust_tol -- boolean indicating whether adjustments can be made to the tolerance
        """
        self.impt_reactions = impt_reactions
        self.impt_reactions_len = len(self.impt_reactions)

        # Optional parameters
        self.tol = tol
        self.dec_cmp = dec_cmp
        self.dec_pnt = dec_pnt
        self.adjust_tol = adjust_tol
        self.bug_report = bug_report

        # Set "Output" flag to 0 to prevent gurobi from printing license information and model statistics.
        self.env = gp.Env(empty=True)
        self.env.setParam("OutputFlag", 0)
        self.env.start()

        self.OPTIMAL = gp.GRB.OPTIMAL
        self.rebuild_model = True
        self.hull_id = 1
        self.solver_count = 0 # Keep track of the number of models built

    def set_stoichiometric_matrix(self, path):
        """
        Read stoichiometric matrix from specified file.
        """
        input_file = open(path, "r")
        rows = []

        for line in input_file.readlines():
            rows.append(line.strip().split())

        input_file.close()

        self.variable_count = len(rows)
        self.reaction_count = len(rows[0])

        self.matrix_Aeq = np.array(rows)
        self.matrix_Beq = np.zeros(self.variable_count)

    def set_reaction_domains(self, path):
        """
        Read upper and lower bounds (domain) for reactions from specified file.
        """
        input_file = open(path, "r")

        self.lower_bounds = []
        self.upper_bounds = []

        for line in input_file.readlines():
            bounds = line.strip().split()
            self.lower_bounds.append(float(bounds[0]))
            self.upper_bounds.append(float(bounds[1]))

        input_file.close()   

    def init_model(self, hyp, matrix_A, matrix_B):
        """
        Initialize a standard Gurobi model. 
        """
        self.model = gp.Model(env=self.env) # Set environment

        self.x = self.model.addMVar((self.reaction_count,),self.lower_bounds, self.upper_bounds)

        self.model.setObjective(hyp @ self.x, gp.GRB.MAXIMIZE)
        self.model.addConstr(self.matrix_Aeq @ self.x == self.matrix_Beq)
        self.model.addConstr(matrix_A @ self.x <= matrix_B)

        self.model.optimize()

    def set_model(self):
        """
        Set the parameters of the model once.
        If this method is called the model will never be rebuilt only adjusted.
        This may lead to possible performance improvements.
        """
        hyp = np.zeros(self.reaction_count)
        matrix_A = np.zeros((self.variable_count, self.reaction_count))
        matrix_B = np.zeros(self.variable_count)

        self.rebuild_model = False
        self.init_model(hyp, matrix_A, matrix_B)

    def build_model(self, hyp, matrix_A, matrix_B):
        """
        Solve model using Gurobi.
        The model can be rebuild each time or only the parameters are adjusted.
        This does not effect the outcome, only performance.
        """
        if not self.rebuild_model:
            self.model.reset() # Delete remaining data from last model build
            self.model.NumObj = 0
            # Remove the second constrait (multiple smaller constraits)
            all_const = self.model.getConstrs()
            del all_const[:self.variable_count]
            self.model.remove(all_const)

            self.model.update()
        
            self.model.setObjective(hyp @ self.x, gp.GRB.MAXIMIZE)
            self.model.addConstr(matrix_A @ self.x <= matrix_B)

            self.model.optimize()
        else:
            self.init_model(hyp, matrix_A, matrix_B)

        sol = np.round(self.x.X, self.dec_pnt) if self.model.status == self.OPTIMAL else 0.0

        self.solver_count += 1

        return (sol, self.model.status)

    def get_hyperplane(self, points):
        """
        Compute the hessian normal form of a set of points.
        """
        points_of_interest = points[self.impt_reactions, :].T # Index important rows then transpose
        distance = np.ones((len(points[0]), 1)) * -1
        C = Matrix(np.hstack((points_of_interest, distance)))
        hess = C.nullspace()
        # First column of the nullspace matrix is the hyperplane
        hyperplane = hess[0]  
        hyp = np.zeros(self.reaction_count)

        for i in range(self.impt_reactions_len):
            hyp[self.impt_reactions[i]] = hyperplane[i] 

        hyp_null = hyperplane[-1]

        return (hyp, hyp_null) 

    def extreme_points(self, hyp, hyp_null, opt):
        """
        Computes an extreme points for a given projection.
        """
        matrix_A = np.array([hyp, -hyp])
        matrix_B = np.array([hyp_null, -hyp_null])
        new_hyp = np.zeros(self.reaction_count)
        new_hyp[self.impt_reactions] = 1.0 * opt
        (x_optimum, sol_flag) = self.build_model(new_hyp, matrix_A, matrix_B)

        if not sol_flag == self.OPTIMAL:
            # If no feasible solution is found, then
            # include tolerance for computing the linear program
            hyp_tol = np.round(hyp,self.dec_cmp)
            matrix_A = np.array([hyp_tol, -hyp_tol])
            matrix_B = np.array([hyp_null + self.tol, -hyp_null + self.tol])
            (x_optimum, sol_flag) = self.build_model(new_hyp, matrix_A, matrix_B)

            if self.adjust_tol:
                new_tol = self.tol
                while sol_flag is not self.OPTIMAL:
                    # Keeps adjusting the tolerance until a solution is found 
                    new_tol = new_tol * 10.0
                    if self.bug_report:
                        print("Tolerance adjusted to {} in extreme point calculation.".format(new_tol))
                    matrix_B = np.array([hyp_null + new_tol, -hyp_null + new_tol])
                    (x_optimum, sol_flag) = self.build_model(new_hyp, matrix_A, matrix_B)
            else:
                # If tolerance is not adjusted, 
                # check that the solution is feasible 
                self.check_sol(sol_flag)

        return x_optimum

    def add_point(self, point_list, point):
        """
        Check if a point is in a list, if not add it. 
        Returns a new reference to the list.
        """
        for i in range(len(point_list[0])):
            if np.array_equal(point_list[:, i], point):
                return point_list
        return np.hstack((point_list, point.reshape(-1, 1)))

    def check_sol(self,sol_flag):
        """
        Checks if the solution is optimal,
        if not raise custom InfeasibleSolution exception.
        """
        if not sol_flag == self.OPTIMAL:
            raise self.InfeasibleSolution

    def initial_points(self):
        """
        Calculate initial points.
        """
        matrix_A = np.zeros((self.variable_count, self.reaction_count))
        matrix_B = np.zeros(self.variable_count)
        hyp = np.zeros(self.reaction_count)

        hyp[self.impt_reactions[0]] = 1.0
        (x_optimum, sol_flag) = self.build_model(hyp, matrix_A, matrix_B)
        self.check_sol(sol_flag)

        hyp_null = x_optimum.dot(hyp)
        all_points = self.extreme_points(hyp, hyp_null, -1.0)

        (x_optimum, sol_flag) = self.build_model(-hyp, matrix_A, matrix_B)
        self.check_sol(sol_flag)
        hyp_null = x_optimum.dot(hyp)
        points = self.extreme_points(hyp, hyp_null, -1.0)

        all_points.shape = (-1, 1) # Change from row vector to column vector
        all_points = self.add_point(all_points, points)

        while len(all_points[0]) <= self.impt_reactions_len:

            (hyp, hyp_null) = self.get_hyperplane(all_points)
            (x_optimum, sol_flag) = self.build_model(-hyp, matrix_A, matrix_B)
            self.check_sol(sol_flag)
            x_hyp = hyp.dot(x_optimum)

            if abs(x_hyp - hyp_null) > self.tol:
                # check if the extreme point has already been calculated
                points = self.extreme_points(hyp, x_hyp, 1.0)
                all_points = self.add_point(all_points, points)
            else:
                (x_optimum, sol_flag) = self.build_model(hyp, matrix_A, matrix_B)
                self.check_sol(sol_flag)
                hyp_null = x_optimum.dot(hyp)
                points = self.extreme_points(hyp, hyp_null, -1.0)
                all_points = self.add_point(all_points, points)

        return all_points

    def list_tolerance(self, ls):
        """
        Check if all elements in list are >= than acceptable tolerance.
        """
        return np.all(abs(ls) < self.tol)
            
    def hyperplane_in_chull(self, hyp, hyp_null, points, hull):
        """
        Check if hyperplane is in convex hull.
        """
        for key, (h, null_h, _, _) in hull.items():
            if (self.list_tolerance(hyp_null - null_h) and
                self.list_tolerance(hyp[self.impt_reactions] - h[self.impt_reactions])):
                return True
        return False

    def initial_hull(self, all_points):
        """
        Generate hull.
        """
        hull = dict() # List of tuples (hyp, hyp_null, points)
        points_range = range(0, len(all_points[0]))
        combs = np.array(list(combinations(points_range, self.impt_reactions_len)))

        for columns in combs:
            # Select points from combination tuple
            points = all_points[:, columns]
            (hyp, hyp_null) = self.get_hyperplane(points)

            if (self.hyperplane_in_chull(hyp, hyp_null, points, hull) or 
                self.hyperplane_in_chull(-hyp, -hyp_null, points, hull)):
                continue

            set_difference = np.setdiff1d(points_range, columns)
            extra_points = hyp.dot(all_points[:, set_difference])

            if round(min(extra_points), self.dec_cmp) >= round(hyp_null, self.dec_cmp):
                hull[self.hull_id] = (-hyp, -hyp_null, points, True)
            else:
                hull[self.hull_id] = (hyp, hyp_null, points, True)
            self.hull_id += 1
            # Possible optimization: store only positive but additional -1 or 1 in tuple

        return hull

    def remove_hyperplane_in_hull(self, hull, all_points):
        """
        Remove hyperplanes laying inside the hull.
        """
        del_keys = []

        for key in list(hull.keys()):
            (hyp, hyp_null, _, _) = hull[key]

            new_points = hyp.dot(all_points)
            new_hyp_null = round(hyp_null, self.dec_cmp)
       
            if (round(min(new_points), self.dec_cmp) < new_hyp_null and 
                round(max(new_points), self.dec_cmp) > new_hyp_null):
                del hull[key]
                del_keys.append(key)

        return(del_keys)
                
    def update_chull(self, hull, all_points, new_point):
        """
        Given a new extreme point, compute all possible hyperplanes with the new point.
        """
        for key in list(hull.keys()):
            (hyp, hyp_null, points, _) = hull[key]
        
            if (np.array_equal(new_point[self.impt_reactions], points[self.impt_reactions])
                or round(hyp.dot(new_point), self.dec_cmp) <= round(hyp_null, self.dec_cmp)):
                continue
           
            for j in range(np.size(points, 1)):
                points_copy = np.copy(points)
                
                points_copy[:, j] = new_point
                (hyp, hyp_null) = self.get_hyperplane(points_copy)

                if (self.hyperplane_in_chull(hyp, hyp_null, points, hull) or 
                    self.hyperplane_in_chull(-hyp, -hyp_null, points, hull)):
                    continue

                new_points = hyp.dot(all_points)
              
                if round(max(new_points), self.dec_cmp) <= round(hyp_null, self.dec_cmp):
                    hull[self.hull_id] = (hyp, hyp_null, points_copy, True)
                    self.hull_id += 1
                    
                    if round(min(new_points), self.dec_cmp) >= round(hyp_null, self.dec_cmp):
                        hull[self.hull_id] = (-hyp, -hyp_null, points_copy, True)
                        self.hull_id += 1 

        del_keys = self.remove_hyperplane_in_hull(hull, all_points) 

        return(del_keys)  

    def get_list_hyperplanes(self, hull):
        """
        Gathers a list of all the hyperplanes in the hull dictionary. 
        """
        hyp_list = []

        for key in hull.keys():
            (hyp, hyp_null, _, _) = hull[key]
            # Concatenate hyp and hyp_null
            hyperplane = list(np.append(np.round(hyp[self.impt_reactions], self.dec_pnt), round(hyp_null, self.dec_pnt)))
            hyp_list.append(hyperplane)
        # Return only unique hyperplanes as some might be repeated
        unique_hyp = [list(x) for x in set(tuple(x) for x in hyp_list)]

        return unique_hyp      

    def incremental_refinement(self, hull, all_points):
        """
        Initial convex hull is refined by maximizing/minimizing the hyperplanes
        containing the extreme points until all the facets of the projection are terminal.

        Returns the final points ad hyperplanes for the specified dimensions. 
        """
        matrix_A = np.zeros((self.variable_count, self.reaction_count))
        matrix_B = np.zeros(self.variable_count)

        keys = [] # Stores the keys of non-terminal hyperplanes
        previous_keys = [] # Stores the keys of the previous non-terminal hyperplanes
        tol_new = self.tol # Can be adjusted if hyperplanes can not be terminated
       
        while True:

            keys.clear()
            counter = 0

            for key in list(hull.keys()):
                (hyp, hyp_null, points, non_terminal) = hull[key]
                if non_terminal:
                    keys.append(key)
            
            if np.array_equal(previous_keys,np.array(keys)):
                # cannot find any new points from the non-terminal
                # hyperplanes, so adjust the tolerance
                if self.adjust_tol:
                    tol_new = self.tol * 10 
                    counter += 1
                    if self.bug_report:
                        print("Tolerance adjusted to {} in incremental refinement.".format(tol_new))
                else:
                    raise SystemExit("Cannot find all terminal hyperplanes. Please adjust tolerance.")
                break

            if len(keys) == 0:
                # break from loop only once all hyperplanes are 
                # marked as terminal 
                break
            
            del_keys = np.empty(0)
            for key in keys:
                if key in del_keys:
                    pass
                else:
                    (hyp, hyp_null, points, _) = hull[key] 
                    # changed hyp sign
                    (x_opt, sol_flag) = self.build_model(hyp, matrix_A, matrix_B)
                    self.check_sol(sol_flag)
                    hyp_x = hyp.dot(x_opt)

                    if abs(hyp_x - hyp_null) < tol_new:
                        hull[key] = (hyp, hyp_null, points, False)
                    else:
                        point = self.extreme_points(hyp, hyp_x, 1.0)
                        all_points = self.add_point(all_points, point)
                        # Keep track of keys that have been removed
                        del_keys = np.append(del_keys,self.update_chull(hull, all_points, point))
            # Check whether the tolerance should be reset
            if counter >= 1 and not np.array_equal(previous_keys,np.array(keys)):
                tol_new = self.tol
                counter = 0
                if self.bug_report:
                    print("Tolerance is reset.")
            # Store keys for the next iteration
            previous_keys = np.array(keys)

        del_keys = self.remove_hyperplane_in_hull(hull, all_points)  

        # Return final points and hyperplanes for the specified dimensions
        final_points = all_points[self.impt_reactions]
        final_planes = self.get_list_hyperplanes(hull)

        return final_points, final_planes

    class InfeasibleSolution(Exception):

        def __init__(self, message="No feasible solution found by solver."):
            """
            Custom exception which can be called when an infeasible solution,
            was found by the solver.
            """
            self.message = message
            super().__init__(self.message)
