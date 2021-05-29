# Convex Hull Method (CHM) 

_________________________________________________________________

## Background 

Using concepts from computational geometry, these CHM implementations allow for the calculation of multi-dimensional production envelopes in large metabolic models. To evluate the impact of numerical precision, this repository contains three different CHM implementations:
- an exact arthimetic implementation in Python
- a floating point, double precision arithmetic implementation in Python 
- a floating point, double precision arithmetic implementation in MATLAB 

The floating point arithmetic implementations find less extreme points overall, due to rounding, but scale easily to higher dimensions on genome-scale metabolic models (tested for up to 6 reactions on interest). To compare precision and run time, we use both implementations to calculate production envelopes of the genome-scale metabolic model of *E. coli* [iJO1366](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3261703/) and its core model [EColiCore2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5206746/). To calculate multi-dimensional production envelopes of the syntropic exchange reactions of a microbial communities, we apply the floating point arithmetic implementation to a [biogas-producing community model](https://biotechnologyforbiofuels.biomedcentral.com/articles/10.1186/s13068-016-0429-x).

Please see the Methods in our associated publication (currently in submission) for further details on the implementations. To use the above implementations, please follow the below outline installation requirements and example runs. 

The double precision Python implementation contains an additional feature, compared to the MATLAB double precision implementation, which stores the stoichiometric constraints of the gurobi model such that there is no need to initialize the gurobi model from scratch each time that a linear program is computed. 

_________________________________________________________________

## Installation  

__Exact Implementation in Python__

All Python requirements can be installed as follows:
```
pip3 install -r req_exact.txt
```

*Notes*
- Python implementation tested on Python 3.6.9 using `QSopt` (version 2.5.10) for exact linear optimization.
- The [python-qsoptex](https://github.com/jonls/python-qsoptex) module requires the [GMP](https://gmplib.org/) and [QSopt_ex](https://github.com/jonls/qsopt-ex) libraries to be installed.
- Debian:
    ```
    apt-get install libgmp3-dev
    apt-get install libqsopt-ex-dev
    ```

__Double Precision Implementation in Python__

All Python requirements can be installed as follows:
```
pip3 install -r req_double.txt
```

*Notes*
- The Python implementation was tested on Python 3.7.6 using [`gurobi`](https://www.gurobi.com/) (version 9.0.2) for double-precision linear optimization.
- A [`gurobi`](https://www.gurobi.com/) license must be downloaded pior to running this implementation

__Double Precision Implementation in MATLAB__

The MATLAB implementation was tested on MATLAB R2020a using [`gurobi`](https://www.gurobi.com/) (version 9.0.2) for double-precision linear optimization. A [`gurobi`](https://www.gurobi.com/) license must be downloaded pior to running this implementation. No additional packages, other than the functions outlined in the `.m` files in the [`DOUBLE/MATLAB/chm`](DOUBLE/MATLAB/chm) folder are required. 

_________________________________________________________________

## Usage

Please see below for individual use cases of the three different implementations. The following input files should be used to run any of the three implementations: 
- `domains.txt`: contains two columns, the 1st with the lower bounds, the 2nd with the upper bounds
- `stoichs.txt`: contains the stoichiometric matrix (reactions x metabolites) of the model
- `reactions.txt`: (optional) contains the name annotations of the reactions as listed in the other two files

Please note that input **values for the exact implementation must be integers** and that all values should be scaled accordingly. 

__Exact Implementation in Python__

Run the [`EXACT/Python/main.py`](EXACT/Python/main.py) file for an example or a variation of the following code:
```python
import chm_exact
reactions = [0, 1] # list of reactions for which to calculate the PE
data_path = "YOUR_DATA_SET_FOLDER"
chm_exact.compute_CH(data_path + "reactions.txt", data_path + "stoichs.txt", \
    data_path + "tdomains.txt", reactions)
```


__Double Precision Implementation in Python__

Run the [`DOUBLE/Python/main.py`](DOUBLE/Python/main.py) file for an example or a variation of the following code:

```python
from chm_double import CHM
reactions = [0, 2, 71] 
chm = CHM(reactions)
data_path = "YOUR_DATA_SET_FOLDER"
chm.set_stoichiometric_matrix(data_path + "stoichs.txt")
chm.set_reaction_domains(data_path + "domains.txt")
chm.set_model()
init_points = chm.initial_points()
hull = chm.initial_hull(init_points)
(final_points, final_hyperplanes) = chm.incremental_refinement(hull, init_points)
print("Final Extreme Points...")
print(final_points)
``` 


__Double Precision Implementation in MATLAB__

Run the [`DOUBLE/MATLAB/examples/EColi/main.m`](DOUBLE/MATLAB/examples/EColi/main.m) file for an example or a variation of the following code:

```MATLAB
data_path = "YOUR_DATA_SET_FOLDER";
domain = load('domain.txt');
reactions = [1, 3]; % list of reactions for which to calculate the PE
CH=computeCH(load(data_path + "stoichs.txt",), domain(:,1), domain(:,2), reactions);
```