import sys
import qsoptex
from sympy import Matrix, sympify
from fractions import Fraction

def solve_lp_exact(obj_inds, opt, h_add, h0_add):
    """
    Solves LP exactly
    """
    global RIDS
    global lp_prob
    flag_a = 0
    lp = lp_prob
    # change objective
    new_obj = {}
    # set integers when possible to speed up computation
    for i in range(len(obj_inds)):
        if sympify(obj_inds[i]).is_integer or obj_inds[i] == 0:
            new_obj[RIDS[i]] = int(obj_inds[i])
        elif sympify(obj_inds[i]).is_rational:
            new_obj[RIDS[i]] = Fraction(str(obj_inds[i]))
    lp.set_linear_objective(new_obj)
    # additional constraints other than stoichiometric, if any
    if h_add and h0_add:
        flag_a = 1
        constr = {}
        for j in range(len(h_add)):
            if h_add[j] != 0:
                if sympify(h_add[j]).is_integer:
                    constr[RIDS[j]] = int(h_add[j])
                elif sympify(h_add[i]).is_rational:
                    constr[RIDS[j]] = Fraction(str(h_add[j]))
        lp.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, constr, rhs=Fraction(str(h0_add[0])))
    if opt == -1:
        lp.set_objective_sense(qsoptex.ObjectiveSense.MAXIMIZE)
    elif opt == 1:
        lp.set_objective_sense(qsoptex.ObjectiveSense.MINIMIZE)
    else:
        sys.exit("opt takes 2 possible values: -1 or 1")
    lp.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = lp.solve()
    # remove last constraint (if added) from the basis LP
    if flag_a:
        lp.delete_linear_constraint(lp.get_constraint_count() - 1)
    if status == qsoptex.SolutionStatus.OPTIMAL:
        return Matrix(lp.get_values())
    else:
        sys.exit("Solver status is not optimal. Status:" + str(status))

def get_hyperplane(pts, dims):
    """
    Compute the Hessian Normal form of a set of points
    """
    h = Matrix.zeros(1, pts.shape[0])
    dis = -Matrix.ones(pts.shape[1], 1)
    pnts_dims = pts[dims, :].T
    C = pnts_dims.col_insert(pnts_dims.shape[1], dis)
    hess = C.nullspace()
    for i in range(len(dims)):
        h[dims[i]] = hess[0][i]
    h0 = hess[0][-1]
    return [h, h0]

def create_lp(polyt, obj_inds):
    """ Creates core LP problem with the Stoichiometric Matrix and list of constraints"""
    # create problem
    p = qsoptex.ExactProblem()
    [Aeq, beq, rids, domain] = [polyt["Aeq"], polyt["beq"], polyt["rids"], polyt["domain"]]
    [lbs, ubs] = domain
    # add variables to lp
    for i in range(len(rids)):
        p.add_variable(name=rids[i], objective=Fraction(str(obj_inds[i])), lower=lbs[i], upper=ubs[i])
    # constraints
    # for each row in S (metabolite) = for each constraint
    for i in range(Aeq.shape[0]):
        constr = {}
        # for each column in S = for each reaction
        for j in range(Aeq.shape[1]):
            if Aeq[i, j] != 0:
                constr[rids[j]] = int(Aeq[i, j])
        p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, constr, rhs=int(beq[i]))
    return p

def read_problem(reactions_path, s_matrix_path, domains_path):
    """
    Read LP problem from 3 files: reactions, Stoichiometric matrix, and constraints
    """
    probl = {}
    # read reaction names
    reac_names = []
    with open(reactions_path, "r") as filetoberead:
        for line in filetoberead.readlines():
            line = line.strip()
            reac_names.append(line)
    probl["rids"] = reac_names
    # read upper and lower bounds of reactions (domain)
    lbs = []
    ubs = []
    with open(domains_path, "r") as filetoberead:
        for line in filetoberead.readlines():
            info = line.strip()
            info = line.split()
            lbs.append(int(info[0]))
            ubs.append(int(info[1]))
    probl["domain"] = [lbs, ubs]
    # read stoichiometric matrix. Rows=metabolites, columns=reactions
    S = []
    with open(s_matrix_path, "r") as filetoberead:
        for line in filetoberead.readlines():
            line = line.strip()
            row = []
            for column in line.split():
                row.append(int(column))
            S.append(row)
    beq = [0] * len(S)
    probl["Aeq"] = Matrix(S)
    probl["beq"] = Matrix(beq)
    return probl

def hp_in_CH(h, h0, v, chull):
    """this function checks if hyperplane and points are already in the CH"""
    flag = 0
    if any([[[h, h0], v] == chull[i][:-1] for i in range(len(chull))]):
        flag = 1
    return flag

def update_CH(new_p, epts, chull, dims):
    """
    Given a new extreme point, compute all possible HP with the new EP
    """
    for i in range(len(chull)):
        pts = chull[i][1]
        if any([pts[dims, p] == new_p[dims, :] for p in range(pts.shape[1])]):
            continue
        bla = chull[i][0][0] * new_p
        if bla[0] <= chull[i][0][1]:
            continue
        for j in range(pts.shape[1]):
            v = pts[:, :]
            v[:, j] = new_p
            [h, h0] = get_hyperplane(v, dims)
            if hp_in_CH(h, h0, v, chull) or hp_in_CH(-h, -h0, v, chull):
                continue
            eh = h * epts
            if max(eh) <= h0:
                chull.append([[h, h0], v, 1])
            else:
                if min(eh) >= h0:
                    chull.append([[-h, -h0], v, 1])
    to_remove = []
    for i in range(len(chull)):
        ec = chull[i][0][0] * epts
        h0 = chull[i][0][1]
        if min(ec) < h0 and max(ec) > h0:
            to_remove.append(i)
    chull = [i for j, i in enumerate(chull) if j not in to_remove]
    return chull

def compute_CH(lp_data, impt_reactions, RIDS): # depends on read_problem, create_lp, initial_points, initial_hull, and incremental_refinement
    """
    Computes the convex hull for production envelopes of metabolic network. Solution is 
    the list of hyperplanes and set of extreme points of the Convex hull. Inputs are:
    * fname: name of file without extension (must be the same for all files
      - fname_r.txt: list of reaction names - order must follow that of S columns
      - fname_S.txt: Stoichiometric matrix
      - fname_d.txt : lb ub for each reaction
    * impt_reactions: list of indices for the dimensions onto which the CH should be computed
    """
    # INITIAL POINTS
    epts = initial_points(impt_reactions, RIDS)
    # INITIAL HULL
    chull = initial_hull(epts, impt_reactions)
    return [chull,epts]

def incremental_refinement(chull, eps, dims): # depends on solve_lp_exact, extreme_point, update_CH
    """
    Refine initial convex hull is refined by maximizing/minimizing the \hps
    containing the \eps until all the facets of the projection are terminal.
    """
    while sum([chull[k][2] for k in range(len(chull))]) != 0:
        for i in range(len(chull)):
            if i >= len(chull):
                break
            h = chull[i][0][0]
            h0 = chull[i][0][1]
            opt = solve_lp_exact(h, -1, [], [])
            hx = h * opt
            if hx[0] == h0:
                chull[i][2] = 0
            else:
                ep = extreme_point(h, hx, -1, dims)
                if not any([eps[dims, j] == ep[dims, :] for j in range(eps.shape[1])]):
                    eps = eps.col_insert(eps.shape[1], ep)
                    chull = update_CH(ep, eps, chull, dims)
        to_remove = []
        for i in range(len(chull)):
            ec = chull[i][0][0] * eps
            h0 = chull[i][0][1]
            if min(ec) < h0 and max(ec) > h0:
                to_remove.append(i)
        chull = [i for j, i in enumerate(chull) if j not in to_remove]
    return [chull, eps]

def initial_points(dims,RIDS): # depends on solve_lp_exact and extreme_point
    """
    Computes Initial set of Extreme Points
    """
    
    num_vars = len(RIDS)
    h = [0] * num_vars
    h[dims[0]] = 1
    h = Matrix([h])
    # max
    opt = solve_lp_exact(h, -1, [], [])
    hx = h * opt
    eps = extreme_point(h, hx, -1, dims)
    # min
    opt = solve_lp_exact(h, 1, [], [])
    hx = h * opt
    ep = extreme_point(h, hx, 1, dims)
    # if extreme point already in the list of EPs
    if not any([eps[dims, j] == ep[dims, :] for j in range(eps.shape[1])]):
        eps = eps.col_insert(eps.shape[1], ep)
    while eps.shape[1] <= len(dims):
        [h, h0] = get_hyperplane(eps, dims)
        opt = solve_lp_exact(h, 1, [], [])
        hx = h * opt
        if hx[0] != h0:
            ep = extreme_point(h, hx, 1, dims)
            if not any([eps[dims, j] == ep[dims, :] for j in range(eps.shape[1])]):
                eps = eps.col_insert(eps.shape[1], ep)
        else:
            opt = solve_lp_exact(h, -1, [], [])
            hx = h * opt
            ep = extreme_point(h, hx, -1, dims)
            if not any([eps[dims, j] == ep[dims, :] for j in range(eps.shape[1])]):
                eps = eps.col_insert(eps.shape[1], ep)
    return eps

def extreme_point(h, h0, optim, dims): # depends on solve_lp_exact
    """
    Computes the extreme point of the projection
    """
    obj = [0] * len(h)
    for i in range(len(dims)):
        obj[dims[i]] = 1
    opt = solve_lp_exact(obj, optim, h, h0)
    return opt

def initial_hull(pnts, dims): # depends on get_hyperplane
    """
    Computes initial hull for the initial set of extreme points
    """
    hull = []
    for i in range(pnts.shape[1]):
        v = pnts[:, :]
        v.col_del(i)
        [h, h0] = get_hyperplane(v, dims)
        if (h * pnts[:, i])[0] >= h0:
            hull.append([[-h, -h0], v, 1])
        else:
            hull.append([[h, h0], v, 1])
    return hull
