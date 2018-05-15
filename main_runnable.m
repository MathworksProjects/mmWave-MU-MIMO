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
    %     nAntennasList = 2.^4;
    nAntennasList = 2.^[3 4];
    nIter = 3;
    plotFLAG = true;
    [CapTot,SINRTot,PrxTot,IntTot] = experiment5(nIter,nAntennasList,plotFLAG);
end
if any(experimentList(:)==51);    experiment51();   end

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

function [CapTot,SINRTot,PrxTot,IntTot] = experiment5(nIter,nAntennasList,plotFLAG)
    % EXPERIMENT 5 - Capacity achieved per #antennas and priority
    % In this experiment, we evaluate the impact of the priority given by
    % the scheduler to users on the number of antennas allocated per users
    % and, in turn, in the capacity achieved in the system.
    %
    %------------- BEGIN CODE EXPERIMENT 1 --------------
    %
    fprintf('Running experiment 5...\n');
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override (problem) parameters
    problem.nUsers = 2;
    problem.MinObjFIsSNR = true;  % (arbitrary)
    problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
    problem.arrayRestriction = 'None';  % Possibilities: "None", "Localized", "Interleaved", "DiagInterleaved"
    % Override (conf) parameters
    conf.algorithm = 'GA';
    conf.PopulationSize_Data = 40;
    conf.Maxgenerations_Data = 20;  % Increase the number of generations for GA
    conf.EliteCount_Data = ceil(conf.PopulationSize_Data/2);
    conf.MaxStallgenerations_Data = conf.Maxgenerations_Data;  % Force it to cover all the generations
    conf.multiPath = false;  % LoS channel (for now)
	% Configure basic parameters
    candSet = (1:1:problem.nUsers);  % Set of users to be considered
	% Create output variables
    CapTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    SINRTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    PrxTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    IntTot = zeros(problem.nUsers,problem.nUsers,length(nAntennasList),nIter);
    % Main execution
    for idxAnt = 1:length(nAntennasList)
        for idxIter = 1:nIter
            fprintf('Iteration %d\n',idxIter);
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
                PrxTot(id,idxAnt,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id),'Azimuth',problem.phiUsers(id),'Type','powerdb');
                fprintf('* Directivity IDmax: %.2f (dB)\n',PrxTot(id,idxAnt,idxIter));
                % Extract interference generated to others (in dB)
                for id1 = 1:1:problem1.nUsers
                    if id1~=id
                        IntTot(id,id1,idxAnt,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id1),'Azimuth',problem.phiUsers(id1),'Type','powerdb');
                        fprintf('  Directivity IDmin(%d): %.2f (dB)\n',id1,IntTot(id,id1,idxAnt,idxIter));
                    end
                end
                problem1.IDUserAssigned = id;
                if plotFLAG
                    % Plot beam pattern obtained with assignation and BF configuration
                    o_plotAssignment_mod(problem1, handle_Conf_Array);
                    % Plot assignation
                    px = problem1.possible_locations(3,:);  % Antenna allocation on x-axis
                    py = problem1.possible_locations(2,:);  % Antenna allocation on y-axis
                    pz = problem1.possible_locations(1,:);  % Antenna allocation on z-axispatch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
                    patch = o_getPatch(problem1.NxPatch,problem1.NyPatch,px,py);
                    arrays = o_getArrays(problem1.nUsers,max(problem1.NmaxArray),W,px,py,pz);
                    o_plot_feasible_comb(problem1,conf,patch,arrays);
                end
            end
            save('temp/exp5-results_so_far','PrxTot','IntTot');
        end
    end
end

function experiment51(varargin)
    % EXPERIMENT 5 - Capacity achieved per #antennas and priority
    % In this experiment, we evaluate the impact of the priority given by
    % the scheduler to users on the number of antennas allocated per users
    % and, in turn, in the capacity achieved in the system.
    %
    %------------- BEGIN CODE EXPERIMENT 1 --------------
    %
    fprintf('Running experiment 51...\n');
    load('temp/playful_results','problem','W');
    % Create handle per user
    problem = o_create_subarray_partition(problem);  %#ok
    problem.NzPatch = problem.NxPatch;
    problem.dz = problem.dx;
    problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                            [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
                            'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque sí
    handle_ConformalArray = phased.URA([problem.NyPatch,problem.NzPatch],...
                            'Lattice','Rectangular','Element',problem.handle_Ant,...
                            'ElementSpacing',[problem.dy,problem.dz]);
    problem.possible_locations = handle_ConformalArray.getElementPosition;

    % Plot beam pattern obtained with assignation and BF configuration
    for id = 1:1:problem.nUsers
        problem.ant_elem = sum(W(id,:)~=0);
        relevant_positions = (W(id,:)~=0);
        Taper_user = W(id,relevant_positions);
        handle_Conf_Array = phased.ConformalArray('Element',problem.handle_Ant,...
                              'ElementPosition',...
                              [zeros(1,problem.ant_elem);...
                              problem.possible_locations(2,relevant_positions);...
                              problem.possible_locations(3,relevant_positions)],...
                              'Taper',Taper_user);
        problem.IDUserAssigned = id;
        o_plotAssignment_mod(problem, handle_Conf_Array);
    end
    % Plot assignation
%     px = problem.possible_locations(3,:);  % Antenna allocation on x-axis
%     py = problem.possible_locations(2,:);  % Antenna allocation on y-axis
%     pz = problem.possible_locations(1,:);  % Antenna allocation on z-axispatch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
%     arrays = o_getArrays(problem.nUsers,max(problem.NmaxArray),W,px,py,pz);
%     o_plot_feasible_comb(problem,conf,patch,arrays);
end