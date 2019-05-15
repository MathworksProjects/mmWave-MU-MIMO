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
maxIter               = 1;  % Maximum number of iterations
configList            = -1;  % antenna configuration (distribution across users) (-1 for random asignation)
problem.N_Antennas    = 8.^2; % Number of antennas
nUserList             = [1 2 4 8];  % Number of users to iterate over
candSet               = (1:1:problem.nUsers);  % Subset to be allocated antennas
conf.detLocation      = true;
conf.useCasesLocation = true;
conf.useCaseLocation  = 3;  % Use-case location ID
conf.verbosity        = 0;  % Use-case location ID
%% Output parameters
Cap_LCMV_tot = zeros(length(configList),1);
Cap_CBF_tot = zeros(length(configList),1);
Cap_HEU_tot = zeros(length(configList),1);
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
%% Override config parameters
conf.Maxgenerations_Data      = 20;
conf.MaxStallgenerations_Data = 10;

for usrIdx = 1:length(nUserList)
    problem.nUsers = nUserList(usrIdx);
    candSet = (1:1:problem.nUsers);
    problem.dUsers = 5.*ones(1,problem.nUsers);
    % Basic configuration - Location of users
    [problem,~,~] = f_configuration(conf,problem);
    % Intermediate variables
    W_LCMV = zeros(length(candSet),length(problem.possible_locations));
    W_CBF = zeros(length(candSet),length(problem.possible_locations));
    availableAnt = (1:1:problem.N_Antennas);
    for id = 1:problem.nUsers
        % Assign antennas
        antennaSelected = randsample(availableAnt,floor(problem.N_Antennas/length(candSet)));
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
    Cap_LCMV_tot(usrIdx) = mean(Cap_LCMV_lin);
    Cap_CBF_tot(usrIdx) = mean(Cap_CBF_lin);
    Cap_HEU_tot(usrIdx) = mean(Cap_HEU_lin);
    fprintf('LCMV: %.3f\tCBF: %.3f\tHEU: %.3f\n',...
          Cap_LCMV_tot(usrIdx),Cap_CBF_tot(usrIdx),Cap_HEU_tot(usrIdx));
% 	% Plot assignation
%     px = problem.possible_locations(3,:);  % Antenna allocation on x-axis
%     py = problem.possible_locations(2,:);  % Antenna allocation on y-axis
%     pz = problem.possible_locations(1,:);  % Antenna allocation on z-axis
%     patch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
%     % Assignation by LCMV
%     arrays = o_getArrays(problem.nUsers,W_LCMV,px,py,pz);
%     o_plot_feasible_comb(problem,conf,patch,arrays);
%     % Assignation by LCMV+HEU
%     arrays = o_getArrays(problem.nUsers,W_HEU,px,py,pz);
%     o_plot_feasible_comb(problem,conf,patch,arrays);
end

% data = [];
% data(:,1,1) = Cap_CBF_tot;
% data(:,2,1) = Cap_LCMV_tot;
% data(:,3,1) = Cap_HEU_tot;
% groupLabels = {'1','2','4','8'};
% figure; hold on;
% xlim([nUserList(1)-0.5 nUserList(end)+0.5]);
% grid minor;
% plotBarStackGroups(data, groupLabels);
% leg = legend('Conventional','LCMV','LCMV + HEU');
% set(leg,'FontSize',12);
% title('Average capacity achieved per user','FontSize',14);
% ylabel('Capacity in bits/Hz/s','FontSize',12);

idx = [1 3 2 4];
xq = (1:1:nUserList(end));
data = [];
data(:,1,1) = interp1(nUserList,Cap_CBF_tot,xq,'linear');
data(:,2,1) = interp1(nUserList,Cap_LCMV_tot(idx),xq,'linear');
data(:,3,1) = interp1(nUserList,Cap_HEU_tot(idx),xq,'linear');
groupLabels = {'1','2','3','4','5','6','7','8'};
figure; hold on;
xlim([nUserList(1)-0.5 nUserList(end)+0.5]);
grid minor;
plotBarStackGroups(data, groupLabels);
leg = legend('Conventional','LCMV','LCMV + HEU');
set(leg,'FontSize',12);
title('Average capacity achieved per user','FontSize',14);
ylabel('Capacity in bits/Hz/s','FontSize',12);

figure; hold on;
Cap_CBF_tot_interp = interp1(nUserList,Cap_CBF_tot,xq,'pchip');
Cap_LCMV_tot_interp = interp1(nUserList,Cap_LCMV_tot,xq,'pchip');
Cap_HEU_tot_interp = interp1(nUserList,Cap_HEU_tot,xq,'pchip');
plot(xq,Cap_CBF_tot_interp,'LineWidth',2);
plot(xq,Cap_LCMV_tot_interp,'LineWidth',2);
plot(xq,Cap_HEU_tot_interp,'LineWidth',2);


% EOF
