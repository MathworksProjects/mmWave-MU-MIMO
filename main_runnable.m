%% Setup workspace
clear; clc; close all;
addpath('utilities','-end');  % Add utilities folder at the end of search path
%% Load basic configuration - static and/or default
problem = o_read_input_problem('data/metaproblem_test.dat');
conf = o_read_config('data/config_test.dat');

% Define several experiments here and override variable values accordingly
%% EXPERIMENT 1 - Capacity offered
% In this experiment we evaluate the capacity that the heuristics are able
% to offer to the devices. Heuristics assigns antennas as a function of the
% priority. The traffic is overloaded. The users location and channel
% varies across simulations.
%
%------------- BEGIN CODE EXPERIMENT 1 --------------
% 
% Override parameters
problem.iat = 60;
problem.deadline = 50;
problem.payload = 1500*8*5e3;
% Configure the simulation environment
[problem,traffic,flows] = f_configuration(conf,problem);
baseFlows = flows;  % For printing purposes at the end of execution
% Main function
[flows,CapTot,TXbitsTot,THTot,lastSlotSimm,lastSelFlow] = main(conf,problem,flows);
% Report of single execution
[ratioOK,ratioNOK] = f_generateReport(flows);
% Plotting of single execution
% main_plotting(problem,TXbitsTot,THTot,baseFlows,lastSelFlow);

%% EXPERIMENT 2 - Chances of achieving the demanded throughput
% To-Do
%% EXPERIMENT 3 - Performance comparisson against null-forcing technique
% To-Do
%% EXPERIMENT 4 - Convergency analysis
% To-Do
