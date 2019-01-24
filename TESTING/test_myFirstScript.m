clear; close all; clc;
addpath('../UTILITIES','-end');  % Add utilities folder at the end of search path
addpath('../code-systems','-end');  % Add system's folder at the end of search path
addpath('../code-beamforming','-end');  % Add beamforming folder at the end of search path
addpath('../code-wirelessEmulation','-end');  % Add channel folder at the end of search path
addpath('../data','-end');  % Add data folder at the end of search path

%% Basic comms parameters
problem = o_read_input_problem('metaproblem_test.dat');
conf = o_read_config('config_test.dat');

%% Input parameters
problem.N_Antennas   = 8.^2; % Number of antennas
problem.nUsers       = 4;  % Number of users
candSet              = (1:1:problem.nUsers);  % Users to be considered for antenna alloc.
conf.detLocation     = false;  % true to use preconfigured locations, false for random
conf.useCaseLocation = false;  % if detLocation is true, retrieve locations
conf.useCaseLocation = 1;  % Use-case location ID
conf.verbosity       = 1;  % Use-case location ID

%% Adjust antenna parameters and create Antenna Array
problem.NxPatch = floor(sqrt(problem.N_Antennas));
problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
problem.NzPatch = problem.NxPatch;
problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                        [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
                        'CosinePower',[1.5 2.5]);
problem.dz = problem.dy;
handle_ConformalArray = phased.URA([problem.NyPatch,problem.NzPatch],...
                        'Lattice','Rectangular','Element',problem.handle_Ant,...
                        'ElementSpacing',[problem.dy,problem.dz]);
problem.possible_locations = handle_ConformalArray.getElementPosition;
elementPos = problem.possible_locations./problem.lambda;

% Basic configuration - Location of users
[problem,~,~] = f_configuration(conf,problem);

% Intermediate variables
W_LCMV = zeros(length(candSet),1);
W_CBF = zeros(length(candSet),1);
availableAnt = (1:1:problem.N_Antennas);

%% First, compute the beamforming weights (CBF and LCMV) with fixed antenna allocation
for id = 1:problem.nUsers
    % Randomly and equally assign antennas to users
    antennaSelected = randsample(availableAnt,problem.N_Antennas/length(candSet));
	availableAnt = setdiff(availableAnt,antennaSelected);
    
    % Retrieve antenna locations
    elementPos1 = elementPos(:,antennaSelected);
    problem.IDUserAssigned = id;
    
    % Retrieve users locations in space
    angles = [-problem.phiUsers ; -problem.thetaUsers];
    
    % Preconfigure beamformers
    sv = steervec(elementPos1,angles);
    Sn = eye(length(antennaSelected));
    resp = zeros(problem.nUsers,1);
    resp(problem.IDUserAssigned) = db2pow(33);
    
    % Call LCMV beamformer
    W_LCMV1 = lcmvweights(sv,resp,Sn);
    W_LCMV(id,antennaSelected) = W_LCMV1.';  % Create global arrray for weights
    
    % Call Conventional (CBF) beamformer
    W2 = cbfweights(elementPos1,angles(:,problem.IDUserAssigned));  % conventional beamformer
    W_CBF(id,antennaSelected) = W2.';  % Create global arrray for weights
    
    % Normalize weigths
    W_LCMV(id,:) = (1/sqrt(W_LCMV(id,:)*W_LCMV(id,:)'))*W_LCMV(id,:);
    W_CBF(id,:) = (1/sqrt(W_CBF(id,:)*W_CBF(id,:)'))*W_CBF(id,:);
end

%% Second, compute the weights using the proposed HELB mechanism
problem.MaxObjF = ones(1,length(candSet));
problem.MinObjF = ones(1,length(candSet));
if conf.MinObjFIsSNR;     problem.MinObjF = 2.^problem.MinObjF - 1;
end
problem.initialW = W_LCMV;
[sol_found,W_HEU,handle_ConformalArray,PRx,I,bestScores] = CBG_solveit(problem,conf,candSet);

% Normalize weights
for id = 1:problem.nUsers; W_HEU(id,:) = (1/sqrt(W_HEU(id,:)*W_HEU(id,:)'))*W_HEU(id,:); end

%% Results

% Parse output results
[~,~,Cap_LCMV_lin,~]  = f_BF_results(W_LCMV,handle_ConformalArray,candSet,problem,conf,false);
[~,~,Cap_CBF_lin,~]  = f_BF_results(W_CBF,handle_ConformalArray,candSet,problem,conf,false);
[~,~,Cap_HEU_lin,~]  = f_BF_results(W_HEU,handle_ConformalArray,candSet,problem,conf,false);

% Parse output
Cap_LCMV_tot = mean(Cap_LCMV_lin);
Cap_CBF_tot = mean(Cap_CBF_lin);
Cap_HEU_tot = mean(Cap_HEU_lin);
fprintf('LCMV: %.3f\tCBF: %.3f\tHELB: %.3f\n',Cap_LCMV_tot,Cap_CBF_tot,Cap_HEU_tot);

% Plot assignation
px = problem.possible_locations(3,:);  % Antenna allocation on x-axis
py = problem.possible_locations(2,:);  % Antenna allocation on y-axis
pz = problem.possible_locations(1,:);  % Antenna allocation on z-axis
patch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
% Assignation by LCMV
arrays = o_getArrays(problem.nUsers,W_LCMV,px,py,pz);
o_plot_feasible_comb(problem,conf,patch,arrays);
% Assignation by LCMV+HEU
arrays = o_getArrays(problem.nUsers,W_HEU,px,py,pz);
o_plot_feasible_comb(problem,conf,patch,arrays);



% EOF
