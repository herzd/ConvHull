#!/usr/bin/env python3
""" library of functions for the exact convhull python algorithm"""

import argparse
from fractions import Fraction
import gzip
import os
import sys
import urllib
import xml.etree.ElementTree
import qsoptex
import sympy

def get_plane_equation(p1,p2,p3):
    '''takes three points given as lists or tuples and
    returns the parameters of the equation of the corresponding
    plane as tuple of four Fraction objects (the first three values make up the normal 
    vector onto the plane)'''
    point1 = numpy.array([Fraction(str(number)) for number in p1])
    point2 = numpy.array([Fraction(str(number)) for number in p2])
    point3 = numpy.array([Fraction(str(number)) for number in p3])
    vector1,vector2 = point3-point1,point2-point1
    normal_vector = numpy.cross(vector1,vector2)
    a,b,c = normal_vector
    d = numpy.dot(normal_vector,point3)
    return a,b,c,d

def intersect_planes(plane1,plane2):
    '''takes two planes, described as equation through its 4 parameters as tuples
    or lists of fraction objects, and returns a point on and the direction vector of 
    their intersection as a tuple of two lists with 3 fraction objects each.'''
    normal_vector1 = numpy.array(plane1[:3])
    normal_vector2 = numpy.array(plane2[:3])
    direction_vector = list(numpy.cross(normal_vector1,
                                        normal_vector2))
    a1,b1,d1 = plane1[0],plane1[1],plane1[3]
    a2,b2,d2 = plane2[0],plane2[1],plane2[3]
    x,y = sympy.symbols('x,y')
    equation_plane1 = sympy.Eq(a1*x+b1*y,d1)
    equation_plane2 = sympy.Eq(a2*x+b2*y,d2)
    point_on_both_planes_sympy = list(sympy.solve([equation_plane1,
                                                   equation_plane2],
                                                  (x,y)).values())
    point_on_both_planes_sympy.append(0)
    point_on_both_planes_fractions = [Fraction(str(number))
                                      for number in point_on_both_planes_sympy]
    return point_on_both_planes_fractions, direction_vector

def intersect_planes_v2(plane1,plane2):
    '''takes two planes, described as equation through its 4 parameters as tuples
    or lists of fraction objects, and returns a point on and the direction vector of 
    their intersection as a tuple of two lists with 3 fraction objects each.'''
    normal_vector1 = numpy.array(plane1[:3])
    normal_vector2 = numpy.array(plane2[:3])
    direction_vector = list(numpy.cross(normal_vector1,
                                        normal_vector2))
    a1,b1,d1 = plane1[0],plane1[1],plane1[3]
    a2,b2,d2 = plane2[0],plane2[1],plane2[3]
    x,y = sympy.symbols('x,y')
    equation_plane1 = sympy.Eq(a1*x+b1*y,d1)
    equation_plane2 = sympy.Eq(a2*x+b2*y,d2)
    point_on_both_planes_sympy = list(sympy.solve([equation_plane1,
                                                   equation_plane2],
                                                  (x,y)).values())
    point_on_both_planes_sympy.append(0)
    point_on_both_planes_fractions = [Fraction(str(number))
                                      for number in point_on_both_planes_sympy]
    return point_on_both_planes_fractions, direction_vector

def download_model(url):
    '''downloads existing model from the web, url has to be a string.
    returns absolute path of downloaded file as string.'''
    filename = url.split("/")[-1:][0]
    filepath = os.path.abspath(filename)
    with urllib.request.urlopen(url) as traveller:
        with open(filepath,"wb") as destination:
            destination.write(traveller.read())
    return filepath

def parse_xml_model(model):
    '''parses existing xml model, given path as string.
    if gunzipped and ends with the extension .gz, unzips the file.
    returns xml.etree.ElementTree object root.'''
    if model.endswith("gz"):
        with gzip.open(model) as unpacked:
            tree = xml.etree.ElementTree.parse(unpacked)
            root = tree.getroot()
    else:
        tree = xml.etree.ElementTree.parse(model)
        root = tree.getroot()
    return root

def extract_metabolites(model):
    '''takes a parsed model in xml.etree.ElementTree.parse-getroot format and returns
    a list of the ids of the metabolites stored in the model.'''
    metabolites_model = model.findall("{http://www.sbml.org/sbml/level3/version1/core}model/"
                                      "{http://www.sbml.org/sbml/level3/version1/core}listOfSpecies/"
                                      "{http://www.sbml.org/sbml/level3/version1/core}species")
    metabolite_list = [metabolite.attrib["id"]
                       for metabolite in metabolites_model]
    return metabolite_list

def extract_stoichiometry(model):
    '''takes a parsed model in xml.etree.ElementTree.parse-getroot format and returns
    a list of tuples containing reaction-name and the negative stoichiometric value as int'''
    # make sure that we match the item we want to process
    stoichiometric_list = []
    reaction_list = model.findall("{http://www.sbml.org/sbml/level3/version1/core}model/"
                                  "{http://www.sbml.org/sbml/level3/version1/core}listOfReactions/")
    for reaction in reaction_list:
        for child in reaction[0]:
            if "listOfProducts" in str(reaction[0]):
                stoichiometric_list.append((reaction.attrib['name'],
                                           child.attrib['species'],
                                           float(child.attrib['stoichiometry'])))
            elif "listOfReactants" in str(reaction[0]):
                stoichiometric_list.append((reaction.attrib['name'],
                                           child.attrib['species'],
                                           -1*float(child.attrib['stoichiometry'])))
    return stoichiometric_list

def extract_parameters(model):
    '''takes a parsed model in xml.etree.ElementTree.parse-getroot format and returns
    a list of tuples containing id and value of the given parameters.'''
    parameters_model = model.findall("{http://www.sbml.org/sbml/level3/version1/core}model/"
                                     "{http://www.sbml.org/sbml/level3/version1/core}listOfParameters/"
                                     "{http://www.sbml.org/sbml/level3/version1/core}parameter")
    parameter_list = [(parameter.attrib["id"],
                       parameter.attrib["value"])
                      for parameter in parameters_model]
    return parameter_list

def extract_reactions(model):
    '''takes a parsed model in xml.etree.ElementTree.parse-getroot format and returns
    a list of tuples containing reaction-name, lower bound, and upper bound'''
    # the following two definitions are just here to keep the linelength in range
    fbc_lb_string = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}lowerFluxBound"
    fbc_ub_string = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}upperFluxBound"
    reactions_model = model.findall("{http://www.sbml.org/sbml/level3/version1/core}model/"
                                    "{http://www.sbml.org/sbml/level3/version1/core}listOfReactions/"
                                    "{http://www.sbml.org/sbml/level3/version1/core}reaction")
    reaction_list = [(reaction.attrib["name"],
                      reaction.attrib[fbc_lb_string],
                      reaction.attrib[fbc_ub_string])
                     for reaction in reactions_model]
    return reaction_list

def resolve_parameters(reaction_list, parameters):
    '''takes the list of reactions containing the unresolved (simply named) parameters
    and replaces them by the actual values from the parameters list (second argument).
    returns the reaction list with filled integer values for bounds .'''
    updated_reaction_list = []
    for reaction in reaction_list:
        for parameter in parameters:
            if reaction[1] == parameter[0]:
                lower_bound_int = int(parameter[1])
            elif reaction[2] == parameter[0]:
                upper_bound_int = int(parameter[1])
        updated_reaction_list.append((reaction[0],
                                      lower_bound_int,
                                      upper_bound_int))
    return updated_reaction_list

def solve_lp_exact(obj_inds, opt, h_add, h0_add, reaction_ids, lp_prob):
    """
    Solves LP exactly
    """
    flag_a = 0
    # change objective
    new_obj = {}
    # set integers when possible to speed up computation
    for pos, value in enumerate(obj_inds):
        if sympy.sympify(obj_inds[pos]).is_integer or obj_inds[pos] == 0:
            new_obj[reaction_ids[pos]] = int(obj_inds[pos])
        elif sympy.sympify(obj_inds[pos]).is_rational:
            new_obj[reaction_ids[pos]] = fractions.Fraction(str(obj_inds[pos]))
    lp_prob.set_linear_objective(new_obj)
    # additional constraints other than stoichiometric, if any
    if h_add and h0_add:
        flag_a = 1
        constr = {}
        for pos,value in enumerate(h_add):
            if h_add[pos] != 0:
                if sympy.sympify(h_add[pos]).is_integer:
                    constr[reaction_ids[pos]] = int(h_add[pos])
                elif sympy.sympify(h_add[pos]).is_rational:
                    constr[reaction_ids[pos]] = fractions.Fraction(str(h_add[pos]))
        lp_prob.add_linear_constraint(qsoptex.ConstraintSense.EQUAL,
                                      constr,
                                      rhs=fractions.Fraction(str(h0_add[0])))
    if opt == -1:
        lp_prob.set_objective_sense(qsoptex.ObjectiveSense.MAXIMIZE)
    elif opt == 1:
        lp_prob.set_objective_sense(qsoptex.ObjectiveSense.MINIMIZE)
    else:
        sys.exit("opt takes 2 possible values: -1 or 1")
    lp_prob.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = lp_prob.solve()
    # remove last constraint (if added) from the basis LP_PROB
    if flag_a:
        lp_prob.delete_linear_constraint(lp_prob.get_constraint_count() - 1)
    if status == qsoptex.SolutionStatus.OPTIMAL:
        return sympy.Matrix(lp_prob.get_values())
    else:
        sys.exit("Solver status is not optimal. Status:" + str(status))

def get_hyperplane(pts, dims):
    """
    Compute the Hessian Normal form of a set of points
    """
    h = sympy.Matrix.zeros(1, pts.shape[0])
    dis = -sympy.Matrix.ones(pts.shape[1], 1)
    pnts_dims = pts[dims, :].T
    C = pnts_dims.col_insert(pnts_dims.shape[1], dis)
    hess = C.nullspace()
    for i in range(len(dims)):
        h[dims[i]] = hess[0][i]
    h0 = hess[0][-1]
    return [h, h0]

def read_problem(reactions_file, stoichiometrics_file, domains_file):
    """
    Read LP problem from 3 files: reactions, Stoichiometric matrix, and constraints
    """
    problem = {}
    # read reaction names
    reac_names = []
    with open(reactions_file, "r") as file_to_read:
        for line in file_to_read.readlines():
            line = line.strip()
            reac_names.append(line)
    problem["rids"] = reac_names
    # read upper and lower bounds of reactions (domain)
    lower_bounds = []
    upper_bounds = []
    with open(domains_file, "r") as file_to_read:
        for line in file_to_read.readlines():
            line = line.strip()
            info = line.split()
            lower_bounds.append(int(info[0]))
            upper_bounds.append(int(info[1]))
    problem["domain"] = [lower_bounds, upper_bounds]
    # read stoichiometric matrix. Rows=metabolites, columns=reactions
    S = []
    with open(stoichiometrics_file, "r") as file_to_read:
        for line in file_to_read.readlines():
            line = line.strip()
            row = []
            for column in line.split():
                row.append(int(column))
            S.append(row)
    beq = [0] * len(S)
    problem["Aeq"] = sympy.Matrix(S)
    problem["beq"] = sympy.Matrix(beq)
    return problem

def create_lp(polyt, obj_inds):
    """ Creates core LP problem with the Stoichiometric Matrix and list of constraints"""
    # create problem
    problem = qsoptex.ExactProblem()
    Aeq,beq,rids = polyt["Aeq"], polyt["beq"], polyt["rids"]
    lower_bound,upper_bound = polyt["domain"]
    # add variables to linear problem
    for i in range(len(rids)):
        problem.add_variable(name=rids[i],
                             objective=fractions.Fraction(str(obj_inds[i])),
                             lower=lower_bound[i],
                             upper=upper_bound[i])
    # constraints
    # for each row in S (metabolite) = for each constraint
    for i in range(Aeq.shape[0]):
        constr = {}
        # for each column in S = for each reaction
        for j in range(Aeq.shape[1]):
            if Aeq[i, j] != 0:
                constr[rids[j]] = int(Aeq[i, j])
        problem.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, constr, rhs=int(beq[i]))
    return problem

def hp_in_chull(h, h0, v, chull):
    """this function checks if hyperplane and points are already in the CH"""
    flag = 0
    if any([[[h, h0], v] == chull[i][:-1] for i in range(len(chull))]):
        flag = 1
    return flag

def update_chull(new_p, epts, chull, dims):
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
            if hp_in_chull(h, h0, v, chull) or hp_in_chull(-h, -h0, v, chull):
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

def compute_chull(lp_data, impt_reactions, reaction_ids, lp_prob):
    """
    Computes the convex hull for production envelopes of metabolic network. Solution is
    the list of hyperplanes and set of extreme points of the Convex hull. Inputs are:
      - lp_data: dictionary containing the keys 'rids', 'Aeq', 'beq'
      - fname_r.txt: list of reaction names - order must follow that of S columns
      - fname_S.txt: Stoichiometric matrix
      - fname_d.txt : lb ub for each reaction
    * impt_reactions: list of indices for the dimensions onto which the CH should be computed
    """
    # INITIAL POINTS
    epts = initial_points(impt_reactions, reaction_ids,lp_prob)
    # INITIAL HULL
    chull = initial_hull(epts, impt_reactions)
    return chull,epts

def incremental_refinement(chull, eps, dims, reaction_ids, lp_prob):
    """
    Refine initial convex hull is refined by maximizing/minimizing the hps
    containing the eps until all the facets of the projection are terminal.
    """
    while sum([chull[k][2] for k in range(len(chull))]) != 0:
        for i in range(len(chull)):
            if i >= len(chull):
                break
            h = chull[i][0][0]
            h0 = chull[i][0][1]
            opt = solve_lp_exact(h, -1, [], [], reaction_ids, lp_prob)
            hx = h * opt
            if hx[0] == h0:
                chull[i][2] = 0
            else:
                ep = extreme_point(h, hx, -1, dims, reaction_ids, lp_prob)
                if not any([eps[dims, j] == ep[dims, :] for j in range(eps.shape[1])]):
                    eps = eps.col_insert(eps.shape[1], ep)
                    chull = update_chull(ep, eps, chull, dims)
        to_remove = []
        for i in range(len(chull)):
            ec = chull[i][0][0] * eps
            h0 = chull[i][0][1]
            if min(ec) < h0 and max(ec) > h0:
                to_remove.append(i)
        chull = [i for j, i in enumerate(chull) if j not in to_remove]
    return chull, eps

def initial_points(dims,reaction_ids, lp_prob): # depends on solve_lp_exact and extreme_point
    """
    Computes Initial set of Extreme Points
    """
    h = [0] * len(reaction_ids)
    h[dims[0]] = 1
    h = sympy.Matrix([h])
    # max
    opt = solve_lp_exact(h, -1, [], [], reaction_ids, lp_prob)
    hx = h * opt
    eps = extreme_point(h, hx, -1, dims, reaction_ids, lp_prob)
    # min
    opt = solve_lp_exact(h, 1, [], [], reaction_ids, lp_prob)
    hx = h * opt
    ep = extreme_point(h, hx, 1, dims, reaction_ids, lp_prob)
    # if extreme point already in the list of EPs
    if not any([eps[dims, j] == ep[dims, :] for j in range(eps.shape[1])]):
        eps = eps.col_insert(eps.shape[1], ep)
    while eps.shape[1] <= len(dims):
        [h, h0] = get_hyperplane(eps, dims)
        opt = solve_lp_exact(h, 1, [], [], reaction_ids, lp_prob)
        hx = h * opt
        if hx[0] != h0:
            ep = extreme_point(h, hx, 1, dims, reaction_ids, lp_prob)
            if not any([eps[dims, j] == ep[dims, :] for j in range(eps.shape[1])]):
                eps = eps.col_insert(eps.shape[1], ep)
        else:
            opt = solve_lp_exact(h, -1, [], [])
            hx = h * opt
            ep = extreme_point(h, hx, -1, dims, reaction_ids, lp_prob)
            if not any([eps[dims, j] == ep[dims, :] for j in range(eps.shape[1])]):
                eps = eps.col_insert(eps.shape[1], ep)
    return eps

def extreme_point(h, h0, optim, dims, reaction_ids, lp_prob): # depends on solve_lp_exact
    """
    Computes the extreme point of the projection
    """
    obj = [0] * len(h)
    for i in range(len(dims)):
        obj[dims[i]] = 1
    opt = solve_lp_exact(obj, optim, h, h0, reaction_ids, lp_prob)
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

def main():
    """
    main procedure using functions from above. the program either takes an URL to a
    downloadable xml(sbml) model, or a filepath to one.  Models may be gzipped.
    If no whole model is available neither online or offline, one can also feed
    the program with three separate files:
    - reaction_file: absolute path to file with one reaction name string in each line.
    - stoichiometric_file: absolute path to file with stoichiometric values for the reactions.
    - boundary_file: absolute path to file with lower and upper bounds for each reaction.
    A dimensionality has to be given as int, to define onto how many dimensions the result
    should be projected.
    This function does not return any object.
    The list of hyperplanes and set of the hull's extreme points of the hull are printed to stdout.
    """
    # argument parser for easier operations with the program
    parser = argparse.ArgumentParser(description='Calculate projected convex set of given network. \
Needs either an xml/sbml-model or 3 files containing the relevant data.')
    parser.add_argument("dimensionality", type = int,
                        help="dimension N onto which convhull shall be projected.")
    parser.add_argument("-m", "--model_location", type=str,
                        help="weburl pointing to sbml(xml)-model, can be a filepath too.")
    parser.add_argument("-r", "--reaction_file", type=str,
                        help="path to file containing reaction ids line by line.")
    parser.add_argument("-s", "--stoichiometric_file", type=str,
                        help="path to file containing stoichiometrics of reactions.")
    parser.add_argument("-b", "--boundary_file", type=str,
                        help="path to file containing upper and lower bounds.")
    arguments = parser.parse_args()
    if arguments.model_location:
        if arguments.model_location.startswith("http"):
            model_path = download_model(arguments.model_location)
            print("downloaded model to:\n\{}\n".format(model_path))
        else:
            model_path = arguments.model_location
            print("loaded model:\n{}\n".format(model_path))
        parsed_model = parse_xml_model(model_path)
        stoichiometry = extract_stoichiometry(parsed_model)
        parameters = extract_parameters(parsed_model)
        reactions = extract_reactions(parsed_model)
        reactions_with_bounds = resolve_parameters(reactions,parameters)
        input_reactions = list(range(arguments.dimensionality))
#        this one here would already deliver all needed information for the calculation
#        print(stoichiometry,reactions_with_bounds)
        # what's happening here dirtily is done to keep structural compatibility with the
        # initial version of the algorithm that still sits in this code. this
        # will be rearranged, once that one is standardized. this will have to be done,
        # as for now it forces to destroy the node-hierarchy of sbml. we are actively
	# losing structural information from the model. The goal is to rearrange
	# the 'polytope = dictionary = problem_read' to hold all the information
	# by reactions. This is technically easy, but the create_lp-function
	# relies on this intermingled dictionary for now.
        read_problem_init = dict.fromkeys(["rids","domain","Aeq","beq"])
        read_problem_init["rids"] = [entry[0] for entry in stoichiometry]
        lower_bounds = [entry[1] for entry in reactions_with_bounds]
        upper_bounds = [entry[2] for entry in reactions_with_bounds]
        read_problem_init["domain"] = [lower_bounds,upper_bounds]
        read_problem_init["Aeq"] = sympy.Matrix([[entry[2] for entry in stoichiometry]])
        read_problem_init["beq"] = sympy.Matrix([0]*len([[entry[2] for entry in stoichiometry]]))
        problem_read = read_problem_init
    else:
#       this here reads the information from files (that most likely derive from
#       a given model).
        reaction_file = arguments.reaction_file
        stoichiometric_file = arguments.stoichiometric_file
        boundary_file = arguments.boundary_file
        problem_read = read_problem(reaction_file,
                                    stoichiometric_file,
                                    boundary_file)
    print("\ninitital problem:\n{}\n".format(problem_read))
    print("projecting hull onto {} dimensions.\n".format(arguments.dimensionality))
    input_reactions = list(range(arguments.dimensionality))
    # create for qsopt_ex (with the toy data this makes a list of 3 zeros)
    objective = [0] * problem_read["Aeq"].shape[1]
    # this sets the first zero to 1
    objective[input_reactions[0]] = 1
    # create linear problem from the dictionary to get rid off global statements within functions
    problem_created = create_lp(problem_read,
                                objective)
    print("created linear problem with qsopt_ex.\n")
    # extract reaction ids, to get rid off the global statements within functions
    reaction_ids = problem_read["rids"]
    chull, epts = compute_chull(problem_read,
                                input_reactions,
                                reaction_ids,
                                problem_created)
    # refine results
    refined_chull, refined_epts = incremental_refinement(chull,
                                                         epts,
                                                         input_reactions,
                                                         reaction_ids,
                                                         problem_created)
    print("refined convex hull:\n{}\n".format(refined_chull))
    print("refined set of points:\n{}\n".format(refined_epts))

if __name__ == "__main__":
    main()
