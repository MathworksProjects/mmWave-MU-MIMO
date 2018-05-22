% Setup the environment
clear; clc; close all;
addpath('utilities','-end');  % Add utilities folder at the end of search path
% Define several experiments here and override variable values accordingly
experimentList = 5;
if any(experimentList(:)==1);    experiment1();   end
if any(experimentList(:)==2);    experiment2();   end
if any(experimentList(:)==3);    experiment3();   end
if any(experimentList(:)==4);    experiment4();   end
if any(experimentList(:)==5)
    nUsers = 2;
%     nAntennasList = [4 5 6 7 8 9 10].^2;
    nAntennasList = [4 5 6].^2;
    nIter = 2;
    plotFLAG = true;
    [Cap,SINR_BB,SINR_PB,DirOK,DirNOK_gntd,DirNOK_pcvd] = ...
                          experiment5(nIter,nUsers,nAntennasList,plotFLAG);
    experiment5_plot(nUsers,nAntennasList,Cap,SINR_BB,SINR_PB,DirOK,DirNOK_gntd,DirNOK_pcvd);
end
if any(experimentList(:)==51)
    load('temp/exp5-results','nUsers','nAntennasList','Cap','SINR_BB','SINR_PB','DirOK','DirNOK_gntd','DirNOK_pcvd');
    experiment5_plot(nUsers,nAntennasList,Cap,SINR_BB,SINR_PB,DirOK,DirNOK_gntd,DirNOK_pcvd);
end

function experiment1(varargin)
    % EXPERIMENT 1 - Capacity offered
    % In this experiment we evaluate the capacity that the heuristics are able
    % to offer to the devices. Heuristics assigns antennas as a function of the
    % priority. The traffic is overloaded. The users location and channel
    % varies across simulations.
    %
    %------------- BEGIN CODE EXPERIMENT 1 --------------
    %
    fprintf('Running experiment 1...\n');
    % Load basic configuration - static and/or default
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override parameters
    problem.iat = 60;
    problem.deadline = 50;
    problem.payload = 1500*8*5e3;
    % Configure the simulation environment
    [problem,~,flows] = f_configuration(conf,problem);
    baseFlows = flows;  % For printing purposes at the end of execution
    % Main function
    [flows,CapTot,TXbitsTot,THTot,lastSlotSimm,lastSelFlow] = main(conf,problem,flows);
    % Report of single execution
    [ratioOK,ratioNOK] = f_generateReport(flows,DEBUG);
    % Plotting of single execution
    % main_plotting(problem,TXbitsTot,THTot,baseFlows,lastSelFlow);
end

% EXPERIMENT 2 - Chances of achieving the demanded throughput
% To-Do. It uses the whole simulator and takes the real traffic as input.

% EXPERIMENT 3 - Performance comparisson against null-forcing technique
% To-Do. The comparisson needs to be agains a more updated technique such
% as the JSDM. Waiting for reply from Kaushik to write to them and get the
% code.

% EXPERIMENT 4 - Convergency analysis
% To-Do. Santi takes care of it.

function [Cap,SINR_BB,SINR_PB,DirOK,DirNOK_gntd,DirNOK_pcvd] = experiment5(nIter,nUsers,nAntennasList,plotFLAG)
    % EXPERIMENT 5 -- 
    % 
    % Aim: Evaluate the average received power (Prx) at the intended users
    % and Interference (Int) generated at other users. The secuential
    % allocation policy leads to an unfair allocation policy, which leads
    % to different Prx and Int values. This experiment analyzes how the
    % size of the antenna array and the place in the priority list impact
    % on the performance of the system
    % 
	% Assumptions (Fixed):
    %   1. Number of antennas: Same across users and prop. to Array size.
    %   2. Number of users: nUsers.
    %   3. User location: From config file.
    %   4. Sub-array geometry: 'None'.
    %   5. Antenna Array geometry: Fixed to URA.
    %   6. Algorithm: GA
    % Variable:
    %   1. Antenna Array size variable: nAntennasList
    %   2. Population size: Prop. to nAntennas in Array
    % 
    % Syntax:  [CapTot,SINRTot,DirOKTot,DirOKAv,DirNOKTot,DirNOKAv] =
    % experiment5(nIter,nUsers,nAntennasList,plotFLAG)
    % 
    % Inputs:
    %    nIter - Number of iterations to extract average values
    %    nUsers - Number of users considered
    %    nAntennaList - Number of antenas
    %    plotFLAG - True for plotting directivity and antenna allocation
    %
    % Outputs: (all have dimensions [nUsers x nAntennasList])
    %    Cap - Capacity in b/Hz/s 
    %    SINR_BB - BaseBand SINR in dB
    %    SINR_PB - PassBand SINR in dB
    %    DirOK - Directivity to the intended transmitter in dB
    %    DirNOK_gntd - Directivity generated to other nodes in dB
    %    DirNOK_pcvd - Directivity perceived due to interfeering nodes in
    %                  dB (nUsers x nAntennasList)
    %
    %------------- BEGIN CODE EXPERIMENT 5 --------------
    %
    fprintf('Running experiment 5...\n');
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override (problem) parameters
    problem.nUsers = nUsers;  % Number of users in the simulation
    problem.MinObjFIsSNR = true;  % (arbitrary)
    problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
    problem.arrayRestriction = 'None';  % Possibilities: "None", "Localized", "Interleaved", "DiagInterleaved"
    % Override (conf) parameters
    conf.verbosity = 1;
    conf.algorithm = 'GA';  % Heuristic algorithm
    conf.NumPhaseShifterBits = 60;  % Number of 
    conf.FunctionTolerance_Data = 1e-10;  % Heuristics stops when not improving solution by this much
    conf.multiPath = false;  % LoS channel (for now)
	% Configure basic parameters
    candSet = (1:1:problem.nUsers);  % Set of users to be considered
	% Create output variables
    CapTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    SINRTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    DirOKTot = -Inf(problem.nUsers,length(nAntennasList),nIter);
    DirNOKTot = -Inf(problem.nUsers,problem.nUsers,length(nAntennasList),nIter);
    % Linearize combinations and asign Population size (To be replaced with
    % convergency analysis values)
%     totComb = log10(problem.nUsers.*factorial(ceil(nAntennasList/problem.nUsers)));
%     maxPop = 70;  % Maximum population size
%     minPop = 40;  % Minimum population size
%     slope = (maxPop - minPop) / (totComb(end)-totComb(1));
%     ordIdx = minPop - slope*totComb(1);
%     PopSizeList = ceil(slope*totComb + ordIdx);
    PopSizeList = 40*ones(length(nAntennasList),1);
    % Main execution
    for idxAnt = 1:length(nAntennasList)
        conf.PopulationSize_Data = PopSizeList(idxAnt);
        conf.Maxgenerations_Data = PopSizeList(idxAnt)*10;
        conf.EliteCount_Data = ceil(conf.PopulationSize_Data/2);
        conf.MaxStallgenerations_Data = conf.Maxgenerations_Data/10;  % Force it to cover all the generations
        for idxIter = 1:nIter
            fprintf('Iteration %d with PopSize %d\n',idxIter,PopSizeList(idxAnt));
            % Configure the simulation environment. Need to place users in new
            % locations and create new channels to have statistically
            % meaningful results
            [problem,~,~] = f_configuration(conf,problem);
            % Select number of antennas
            problem.N_Antennas = nAntennasList(idxAnt);
            % Adjust parameters
            problem.NxPatch = floor(sqrt(problem.N_Antennas));
            problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
            problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
            % Call heuristics
            fprintf('\t** %d Antennas and %d Users...\n',problem.N_Antennas,problem.nUsers);
            [~,W,~,estObj] = f_heuristics(problem,conf,candSet);
            % Heuristics - Post Processing
            if conf.MinObjFIsSNR;     CapTot(:,idxAnt,idxIter)  = log2(estObj+1);  % in bps/Hz
                                      SINRTot(:,idxAnt,idxIter) = 10*log10(estObj);  % in dB
            else;                     CapTot(:,idxAnt,idxIter)  = estObj;  % in bps/Hz
                                      SINRTot(:,idxAnt,idxIter) = 10*log10(2.^(estTH/problem.Bw) - 1);  % in dB
            end
            % Reconstruct array
            % Create handle per user
            problem1 = o_create_subarray_partition(problem);
            problem1.NzPatch = problem1.NxPatch;
            problem1.dz = problem1.dx;
            problem1.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                                    [problem1.freq-(problem1.Bw/2) problem1.freq+(problem1.Bw/2)],...
                                    'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque sí
            handle_ConformalArray = phased.URA([problem1.NyPatch,problem1.NzPatch],...
                                    'Lattice','Rectangular','Element',problem1.handle_Ant,...
                                    'ElementSpacing',[problem1.dy,problem1.dz]);
            problem1.possible_locations = handle_ConformalArray.getElementPosition;
            for id = 1:1:problem1.nUsers
                problem1.ant_elem = sum(W(id,:)~=0);
                relevant_positions = (W(id,:)~=0);
                Taper_user = W(id,relevant_positions);
                handle_Conf_Array = phased.ConformalArray('Element',problem1.handle_Ant,...
                                      'ElementPosition',...
                                      [zeros(1,problem1.ant_elem);...
                                      problem1.possible_locations(2,relevant_positions);...
                                      problem1.possible_locations(3,relevant_positions)],...
                                      'Taper',Taper_user);
                % Extract Rx Power (in dB)
                DirOKTot(id,idxAnt,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id),'Azimuth',problem.phiUsers(id),'Type','powerdb');
                fprintf('* Directivity IDmax: %.2f (dB)\n',DirOKTot(id,idxAnt,idxIter));
                % Extract interference generated to others (in dB)
                for id1 = 1:1:problem1.nUsers
                    if id1~=id
                        DirNOKTot(id,id1,idxAnt,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id1),'Azimuth',problem.phiUsers(id1),'Type','powerdb');
                        fprintf('  Directivity IDmin(%d): %.2f (dB)\n',id1,DirNOKTot(id,id1,idxAnt,idxIter));
                    end
                end
                problem1.IDUserAssigned = id;
                if plotFLAG
                    % Plot beam pattern obtained with assignation and BF configuration
                    o_plotAssignment_mod(problem1, handle_Conf_Array);
                end
            end
            if plotFLAG
                % Plot assignation
                px = problem1.possible_locations(3,:);  % Antenna allocation on x-axis
                py = problem1.possible_locations(2,:);  % Antenna allocation on y-axis
                pz = problem1.possible_locations(1,:);  % Antenna allocation on z-axispatch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
                patch = o_getPatch(problem1.NxPatch,problem1.NyPatch,px,py);
                arrays = o_getArrays(problem1.nUsers,W,px,py,pz);
                o_plot_feasible_comb(problem1,conf,patch,arrays);
            end
            save('temp/exp5-results_so_far','DirOKTot','DirNOKTot','nUsers','nAntennasList');
        end
    end
    % Convert back to Watts (from dB)
    DirOKTot_lin = 10.^(DirOKTot./10);
    DirNOKTot_lin = 10.^(DirNOKTot./10);
    % Compute average Directivities
    DirOK_lin = zeros(nUsers,length(nAntennasList));  % Directivity generated by intended user
    DirNOK_gntd_lin = zeros(nUsers,length(nAntennasList));  % Generated interference by intended user
    DirNOK_pcvd_lin = zeros(nUsers,length(nAntennasList));  % Perceived interference by intended user
    for antIdx = 1:length(nAntennasList)
        DirOK_lin(:,antIdx) = mean(DirOKTot_lin(:,antIdx,:),3);
        DirNOK_gntd_lin(:,antIdx) = sum(mean(DirNOKTot_lin(:,:,antIdx,:),4),1); % Generated interference 
        DirNOK_pcvd_lin(:,antIdx) = sum(mean(DirNOKTot_lin(:,:,antIdx,:),4),2); % Perceived interference
    end
    DirOK = 10*log10(DirOK_lin);  % Directivity generated to intended user
    DirNOK_gntd = 10*log10(DirNOK_gntd_lin);  % Directivity being generated by intended user
    DirNOK_pcvd = 10*log10(DirNOK_pcvd_lin);  % Directivity inflicted to intended user
    % Compute SINR and Capacities
    chLoss = 10*log10( ((4*pi*problem.dUsers(1:nUsers)) ./ problem.lambda).^2 ).';  % Losses
    chLoss = repmat(chLoss,1,length(nAntennasList));
    Ptx = repmat(problem.Ptx,1,length(nAntennasList));  % Initial transmit power
    SINR_PB = Ptx + DirOK - DirNOK_pcvd - chLoss - problem.Noise;  % Compute SINR Pass-Band (PB)
    SINR_BB = mean(SINRTot,3);  % Compute SINR Base-Band (BB)
    Cap = mean(CapTot,3);  % Compute Average Capacities in the system
    save('temp/exp5-results','Cap','SINR_BB','SINR_PB','DirOK','DirNOK_gntd','DirNOK_pcvd','DirOKTot','DirNOKTot','nUsers','nAntennasList');
end

function experiment5_plot(nUsers,nAntennasList,Cap,SINR_BB,SINR_PB,DirOK,DirNOK_gntd,DirNOK_pcvd)
    % EXPERIMENT 5 - Plotting results
    % Get figure number
    h = findobj('type','figure');
    figNum = length(h) + 1;
    % Plot Directivities
    figure(figNum);  figNum = figNum + 1;
    leg = cell(nUsers,1);
    for id = 1:nUsers
        subplot(1,3,1); hold on;
        plot(nAntennasList,DirOK(id,:),'LineWidth',2,'Marker','s');
        subplot(1,3,2); hold on;
        plot(nAntennasList,DirNOK_gntd(id,:),'LineWidth',2,'Marker','s');
        subplot(1,3,3); hold on;
        plot(nAntennasList,DirNOK_pcvd(id,:),'LineWidth',2,'Marker','s');
        leg{id} = cell2mat(strcat('user',{' '},num2str(id)));
    end
    subplot(1,3,1);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    title('Directivity to intended user','FontSize',12);
    legend(leg,'FontSize',12);
    subplot(1,3,2);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    title('Interference generated to other users','FontSize',12);
    legend(leg,'FontSize',12);
    subplot(1,3,3);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    title('Interference generated to intended user','FontSize',12);
    legend(leg,'FontSize',12);
    % Plot perceived SINRs
    figure(figNum);  figNum = figNum + 1;
    for id = 1:nUsers
        subplot(1,2,1); hold on;
        plot(nAntennasList,SINR_BB(id,:),'LineWidth',2,'Marker','s');
        subplot(1,2,2); hold on;
        plot(nAntennasList,SINR_PB(id,:),'LineWidth',2,'Marker','s');
    end
    subplot(1,2,1);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    title('Base-Band (BB) SINR','FontSize',12);
    legend(leg,'FontSize',12);
    subplot(1,2,2);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    title('Pass-Band (PB) SINR','FontSize',12);
    legend(leg,'FontSize',12);
    % Plot perceived Capacities
    figure(figNum);  figNum = figNum + 1;                              %#ok
    for id = 1:nUsers
        hold on;
        plot(nAntennasList,Cap(id,:),'LineWidth',2,'Marker','s');
    end
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Capacity in bits/Hz/s','FontSize',12);
    title('Capacity achieved in the system','FontSize',12);
    legend(leg,'FontSize',12);
end