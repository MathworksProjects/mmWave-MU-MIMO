% Setup the environment
clear; clc; close all;
addpath('utilities','-end');  % Add utilities folder at the end of search path
% Define several experiments here and override variable values accordingly
experimentList = 5;
if any(experimentList(:)==1);    experiment1();   end
if any(experimentList(:)==2);    experiment2();   end
if any(experimentList(:)==3);    experiment3();   end
if any(experimentList(:)==4);    experiment4();   end
if any(experimentList(:)==5);    experiment5();   end

function experiment1(varargin)
    %% EXPERIMENT 1 - Capacity offered
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

%% EXPERIMENT 2 - Chances of achieving the demanded throughput
% To-Do. It uses the whole simulator and takes the real traffic as input.
%% EXPERIMENT 3 - Performance comparisson against null-forcing technique
% To-Do. The comparisson needs to be agains a more updated technique such
% as the JSDM. Waiting for reply from Kaushik to write to them and get the
% code.
%% EXPERIMENT 4 - Convergency analysis
% To-Do. Santi takes care of it.
function experiment5(varargin)
    %% EXPERIMENT 5 - Capacity achieved per #antennas and priority
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
    % Override parameters
    problem.nUsers = 3;
    problem.MinObjFIsSNR = true;
%     problem.MinObjF = (problem.nUsers:-1:1);
    delta = 0.1;
    SNRdBTops = 30;  % in dB
    problem.MinObjF = delta.*(problem.nUsers:-1:1)*(10^(SNRdBTops/10));
    candSet = (1:1:problem.nUsers);
%     N_AntennasList = 2.^[3 4 5 6 7 8 9 10];
    N_AntennasList = 2.^8;
    % Configure the simulation environment
    [problem,~,~] = f_configuration(conf,problem);
    % Main execution
    for idx = 1:length(N_AntennasList)
        % Select number of antennas
        problem.N_Antennas = N_AntennasList(idx);
        % Adjust parameters
        problem.NxPatch = floor(sqrt(problem.N_Antennas));
        problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
        problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
        % Call heuristics
        fprintf('Solving problem with %d antennas and %d users...\n',problem.N_Antennas,problem.nUsers);
        [sol_found,W,~,~] = f_heuristics(problem,conf,candSet);
        save('temp/playful_results');
        %%
        load('temp/playful_results');
        fprintf('Solution found: %d\n',sol_found);
        % Create handle per user
        problem = o_create_subarray_partition(problem);
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
        px = problem.possible_locations(3,:);  % Antenna allocation on x-axis
        py = problem.possible_locations(2,:);  % Antenna allocation on y-axis
        pz = problem.possible_locations(1,:);  % Antenna allocation on z-axispatch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
        arrays = o_getArrays(problem.nUsers,max(problem.NmaxArray),W,px,py,pz);
        o_plot_feasible_comb(problem,conf,patch,arrays);
    end
end