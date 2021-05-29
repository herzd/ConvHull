%This script requires gurobi and the cobratoolbox
% Load the Model
model = readCbModel('KochModel_Constrained');

%% Construct the new Stoichiometric Matrix

Aeq = model.S;
or = size(model.S,1);
% get the reaction IDs for constraining f_DV, f_MB, f_MM
mumax = 0.05213;
DV{1} = findRxnIDs(model,[{'ATPdrain_DV'},{'Ac_DV__Ac_EX'},{'BM_DV'}]);
DV{2} = [4.3 50 mumax];
MB{1} = findRxnIDs(model,[{'ATPdrain_MB'},{'Meth_MB__Meth_EX'},{'BM_MB'}]);
MB{2} = [2.5 15 mumax];
MM{1} = findRxnIDs(model,[{'ATPdrain_MM'},{'Meth_MM__Meth_EX'},{'BM_MM'}]);
MM{2} = [0.9 15 mumax];
Constr = {DV,MB,MM};
outR = findRxnIDs(model,[{'BM_DV'},{'BM_MB'},{'BM_MM'}]);

% Constrain out-going hydrogen from BM and MM to zero
constr = findRxnIDs(model,[{'H2_MM__H2_EX'},{'H2_MB__H2_EX'}]);
model.ub(constr) = 0.0;

% add additional constraints to S
for i=1:size(Constr,2)
    new_col = zeros(size(Aeq,1),1);
    Aeq = [Aeq new_col]; % add a column of zeros to S
    for j=1:size(Constr{i}{1},2)
        new_row = zeros(1,size(Aeq,2));
        idx = Constr{i}{1}(j);
        new_row(idx) = 1;
        new_row(end) = -Constr{i}{2}(j);
        Aeq = [Aeq; new_row]; % add row with fj constraints 
    end
end
%Set constraint that sum(fj)=1
new_col = zeros(size(Aeq,1),1);
Aeq = [Aeq new_col]; % add a column of zeros to S
new_row = zeros(1,size(Aeq,2));
new_row([(length(new_row)-3):(length(new_row)-1)]) = [1 1 1];
new_row(end) = -1;
Aeq = [Aeq; new_row];

%Set constraint that sum(BMj)=mumax
new_col = zeros(size(Aeq,1),1);
Aeq = [Aeq new_col]; % add a column of zeros to S
new_row = zeros(1,size(Aeq,2));
new_row(outR) = [1 1 1];
new_row(end) = -1;
Aeq = [Aeq; new_row];

% Bound the f_DV, f_MB, f_MM and sum(fj)
lbs = [model.lb; zeros(3,1); 1; mumax*0.9995];
ubs = [model.ub; ones(4,1); mumax];

%% Run the CHM
p1 = 16;
p2 = 12;

biomass = [332,330];
nutrients = findRxnIDs(model,[{'H2_EX__H2_MB'},{'H2_EX__H2_MM'},{'CO2_Medium'},{'Eth_Medium_in'},{'Meth_Medium'}]);
dims = [biomass,nutrients];
%dims = nutrients;
sg = [repmat('=',or,1);repmat(['>';'<';'<'],3,1);'=';'='];

tic
CH=computeCH(Aeq,lbs,ubs,dims,p1,p2,1,1,sg);
final_time = toc;
disp(toc)
