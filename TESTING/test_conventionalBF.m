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
maxIter              = 1;  % Maximum number of iterations
configList           = [1 2 3 4 -1];  % antenna configuration (distribution across users) (-1 for random asignation)
problem.N_Antennas   = 4.^2; % Number of antennas
problem.nUsers       = 4;  % Number of users
candSet              = (1:1:problem.nUsers);  % Subset to be allocated antennas
conf.detLocation     = false;
conf.useCaseLocation = false;
conf.useCaseLocation = 1;  % Use-case location ID
conf.verbosity       = 1;  % Use-case location ID
%% Output parameters
Cap_LCMV_tot = zeros(length(configList),1);
Cap_CBF_tot = zeros(length(configList),1);
Cap_HEU_tot = zeros(length(configList),1);
%% Antenna Configurations
% Configuration 1
antennaConf(1,:,1) = [1,2,3,4]; antennaConf(1,:,2) = [5,6,7,8]; antennaConf(1,:,3) = [9,10,11,12]; antennaConf(1,:,4) = [13,14,15,16];
% Configuration 2
antennaConf(2,:,1) = [1,3,9,11];  antennaConf(2,:,2) = [2,4,10,12];  antennaConf(2,:,3) = [5,7,13,15];  antennaConf(2,:,4) = [6,8,14,16];
% Configuration 3
antennaConf(3,:,1) = [1,4,13,16];  antennaConf(3,:,2) = [6,7,10,11];  antennaConf(3,:,3) = [5,9,8,12];  antennaConf(3,:,4) = [2,3,14,15];
% Configuration 4
antennaConf(4,:,1) = [1,5,9,13];  antennaConf(4,:,2) = [2,6,10,14];  antennaConf(4,:,3) = [3,7,11,15];  antennaConf(4,:,4) = [4,8,12,16];
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

%%
for idxConf = 1:length(configList)
    configuration = configList(idxConf);
    Cap_LCMV_tot1 = zeros(maxIter,1);
    Cap_CBF_tot1 = zeros(maxIter,1);
    Cap_HEU_tot1 = zeros(maxIter,1);
    for iter = 1:maxIter
        % Intermediate variables
        W_LCMV = zeros(length(candSet),length(problem.possible_locations));
        W_CBF = zeros(length(candSet),length(problem.possible_locations));
        availableAnt = (1:1:problem.N_Antennas);
        for id = 1:problem.nUsers
            if configuration == -1
                antennaSelected = randsample(availableAnt,problem.N_Antennas/length(candSet));
                availableAnt = setdiff(availableAnt,antennaSelected);
            else
                antennaSelected = antennaConf(configuration,:,id);
            end
            elementPos1 = elementPos(:,antennaSelected);
            problem.IDUserAssigned = id;
            angles = [-problem.phiUsers ; -problem.thetaUsers];
            sv = steervec(elementPos1,angles);
            Sn = eye(length(antennaSelected));
            resp = zeros(problem.nUsers,1);
            resp(problem.IDUserAssigned) = db2pow(33);
            % Call LCMV
            W_LCMV1 = lcmvweights(sv,resp,Sn);
            W_LCMV(id,antennaSelected) = W_LCMV1.';  % Create global arrray for weights
            % Call Conventional (CBF)
            W2 = cbfweights(elementPos1,angles(:,problem.IDUserAssigned));  % conventional beamformer
            W_CBF(id,antennaSelected) = W2.';  % Create global arrray for weights
            % Normalize weigths
            W_LCMV(id,:) = (1/sqrt(W_LCMV(id,:)*W_LCMV(id,:)'))*W_LCMV(id,:);
            W_CBF(id,:) = (1/sqrt(W_CBF(id,:)*W_CBF(id,:)'))*W_CBF(id,:);
        end
        % Call LCMV + HEU
        problem.MaxObjF = ones(1,length(candSet));
        problem.MinObjF = ones(1,length(candSet));
        if conf.MinObjFIsSNR;     problem.MinObjF = 2.^problem.MinObjF - 1;
        end
        problem.initialW = W_LCMV;
        [sol_found,W_HEU,handle_ConformalArray,PRx,I,bestScores] = CBG_solveit(problem,conf,candSet);
        % Normalize weights
        for id = 1:problem.nUsers; W_HEU(id,:) = (1/sqrt(W_HEU(id,:)*W_HEU(id,:)'))*W_HEU(id,:); end
        % Parse output results
        [~,~,Cap_LCMV_lin,~]  = f_BF_results(W_LCMV,handle_ConformalArray,candSet,problem,conf,false);
        [~,~,Cap_CBF_lin,~]  = f_BF_results(W_CBF,handle_ConformalArray,candSet,problem,conf,false);
        [~,~,Cap_HEU_lin,~]  = f_BF_results(W_HEU,handle_ConformalArray,candSet,problem,conf,false);
        % Parse output
        Cap_LCMV_tot1(iter) = mean(Cap_LCMV_lin);
        Cap_CBF_tot1(iter) = mean(Cap_CBF_lin);
        Cap_HEU_tot1(iter) = mean(Cap_HEU_lin);
    end
    Cap_LCMV_tot(idxConf) = mean(Cap_LCMV_tot1);
    Cap_CBF_tot(idxConf) = mean(Cap_CBF_tot1);
    Cap_HEU_tot(idxConf) = mean(Cap_HEU_tot1);
    fprintf('Conf %d\tLCMV: %.3f\tCBF: %.3f\tHEU: %.3f\n',configuration,...
          Cap_LCMV_tot(idxConf),Cap_CBF_tot(idxConf),Cap_HEU_tot(idxConf));
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
end

data(:,1,1) = Cap_CBF_tot;
data(:,2,1) = Cap_LCMV_tot;
data(:,3,1) = Cap_HEU_tot;
groupLabels = {'Alloc. 1','Alloc. 2','Alloc. 3','Alloc. 4','Random'};
figure; hold on;
xlim([0.5 length(configList)+0.5]);
grid minor;
plotBarStackGroups(data, groupLabels);
leg = legend('Conventional','LCMV','LCMV + HEU');
set(leg,'FontSize',12);
title('Average capacity achieved per user','FontSize',14);
ylabel('Capacity in bits/Hz/s','FontSize',12);


%% 
problem = o_read_input_problem('metaproblem_test.dat');
conf = o_read_config('config_test.dat');
%% Input parameters
maxIter              = 1;  % Maximum number of iterations
configList           = -1;  % antenna configuration (distribution across users) (-1 for random asignation)
problem.N_Antennas   = 4.^2; % Number of antennas
problem.nUsers       = 4;  % Number of users
candSet              = (1:1:problem.nUsers);  % Subset to be allocated antennas
conf.detLocation     = false;
conf.useCaseLocation = false;
conf.useCaseLocation = 1;  % Use-case location ID
conf.verbosity       = 1;  % Use-case location ID
%% Output parameters
Cap_LCMV_tot = zeros(length(configList),1);
Cap_CBF_tot = zeros(length(configList),1);
Cap_HEU_tot = zeros(length(configList),1);
nUserList = [2 4 6 8];
configList
for idxConf = 1:length(configList)
    for usrIdx = 1:length(nUserList)
        problem.nUsers = nUserList(usrIdx);
        % Intermediate variables
        W_LCMV = zeros(length(candSet),length(problem.possible_locations));
        W_CBF = zeros(length(candSet),length(problem.possible_locations));
        availableAnt = (1:1:problem.N_Antennas);
        for id = 1:problem.nUsers
            % Assign antennas
            antennaSelected = randsample(availableAnt,problem.N_Antennas/length(candSet));
            availableAnt = setdiff(availableAnt,antennaSelected);
            % Start computing beamforming algorithm
            elementPos1 = elementPos(:,antennaSelected);
            problem.IDUserAssigned = id;
            angles = [-problem.phiUsers ; -problem.thetaUsers];
            sv = steervec(elementPos1,angles);
            Sn = eye(length(antennaSelected));
            resp = zeros(problem.nUsers,1);
            resp(problem.IDUserAssigned) = db2pow(33);
            % Call LCMV
            W_LCMV1 = lcmvweights(sv,resp,Sn);
            W_LCMV(id,antennaSelected) = W_LCMV1.';  % Create global arrray for weights
            % Call Conventional (CBF)
            W2 = cbfweights(elementPos1,angles(:,problem.IDUserAssigned));  % conventional beamformer
            W_CBF(id,antennaSelected) = W2.';  % Create global arrray for weights
            % Normalize weigths
            W_LCMV(id,:) = (1/sqrt(W_LCMV(id,:)*W_LCMV(id,:)'))*W_LCMV(id,:);
            W_CBF(id,:) = (1/sqrt(W_CBF(id,:)*W_CBF(id,:)'))*W_CBF(id,:);
        end
        % Call LCMV + HEU
        problem.MaxObjF = ones(1,length(candSet));
        problem.MinObjF = ones(1,length(candSet));
        if conf.MinObjFIsSNR;     problem.MinObjF = 2.^problem.MinObjF - 1;
        end
        problem.initialW = W_LCMV;
        [sol_found,W_HEU,handle_ConformalArray,PRx,I,bestScores] = CBG_solveit(problem,conf,candSet);
        % Normalize weights
        for id = 1:problem.nUsers; W_HEU(id,:) = (1/sqrt(W_HEU(id,:)*W_HEU(id,:)'))*W_HEU(id,:); end
        % Parse output results
        [~,~,Cap_LCMV_lin,~]  = f_BF_results(W_LCMV,handle_ConformalArray,candSet,problem,conf,false);
        [~,~,Cap_CBF_lin,~]  = f_BF_results(W_CBF,handle_ConformalArray,candSet,problem,conf,false);
        [~,~,Cap_HEU_lin,~]  = f_BF_results(W_HEU,handle_ConformalArray,candSet,problem,conf,false);
        % Parse output
        Cap_LCMV_tot1(iter) = mean(Cap_LCMV_lin);
        Cap_CBF_tot1(iter) = mean(Cap_CBF_lin);
        Cap_HEU_tot1(iter) = mean(Cap_HEU_lin);
    end
end


% EOF
