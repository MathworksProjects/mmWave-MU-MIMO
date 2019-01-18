function [problem,traffic,flows] = f_configuration(conf,problem)
% f_configuration - Sets the stage for the DES (main function) by
% generating PHY-flows (application bits) per user, generating the 5G
% wireless channel, etcetera.
%
% Syntax:  [problem,traffic,flows] = f_configuration(conf,problem)
%
% Inputs:
%    conf - struct containint configuration in data/config_test.dat
%    problem - struct containint configuration in data/metaproblem_test.dat
%
% Outputs:
%    problem - Updated configuration
%    traffic - todo
%    flows - Array of structs of length equal to the number of
%            users. For each user, each flow (belonging to each packet) is
%            characterized by the amount of bits that needs to be delivered. The
%            amount of bits are distributed uniformly across the slots until
%            reaching the deadline. Thus, the variable flows contains four features:
%            - slots:     For each flow, the slots across it.
%            - TH:        The average throughput demanded for that flow over the slots
%                         indicated in the slots field.
%            - remaining: We set the remaining field of the flow to be the value of
%                         the Payload.
%            - deadlines: The deadline of the actual packet in slot ID. 
%                         This is used in the future to set priorities.
%            - failed:    Mark whether the flow has failed to be served before the
%                         deadline (0 or 1).
%            - success:   Mark whether the flow has been succesfully served before the
%                         deadline (0 or 1).
%            - maxSlot:   Maximum deadline slot used as simulation time or Tsym in
%                         the system in case trafficType is 'dataSet'.
%
% Example:
%           addpath('data','-end');
%           problem = o_read_input_problem('metaproblem_test.dat');
%           conf = o_read_config('config_test.dat');
%           [problem,~,flows] = f_configuration(conf, problem);  % Struct with configuration parameters
%
% Other m-files required: o_generate_5Gchannels, o_generate_positions, f_configureTraffic
%                         f_genDetTraffic, f_arrivalToFlow
% Subfunctions: None
% MAT-files required: 'TABLE-SNR-loc4.mat'
%
% See also: main
%
%------------- BEGIN CODE --------------
% Geographic distribution of users
[problem.thetaUsers, problem.phiUsers, problem.dUsers] = ...
    o_generate_positions(conf, problem);
% Generate channels per user
if conf.Use5GChannel
    % 5G 3GPP ETSI TR 38.901 compliant channel model
    [problem.fullChannels, problem.thetaChannels, problem.phiChannels, problem.alphaChannels] = ...
                        o_generate_5Gchannels(conf,problem.nUsers,...
                        problem.thetaUsers, problem.phiUsers, ...
                        problem.NxPatch, problem.NyPatch, problem.freq,...
                        problem.maxnChannelPaths);
else
    % Basic channel model (for the moment, including channel gains (alphas) !!!
    [problem.thetaChannels, problem.phiChannels, problem.alphaChannels] = ...
                        o_generate_channels(conf,problem.nUsers,...
                        problem.maxnChannelPaths);
end
% Load the PER-MCS table in the problem struct
if strcmp(conf.DelayProfile,'CDL-D'); load('MCSPERTable-CDLD-SC.mat','mcsTable');
elseif strcmp(conf.DelayProfile,'CDL-C'); load('MCSPERTable-CDLC-SC.mat','mcsTable');
end
%     load('data/MCSPERTable.mat','mcsTable');  % Simplistic channel model
problem.MCSPER = mcsTable;
% Load the SNR table for faster execution
load('TABLE-SNR-loc4.mat','SINRLCMV','SINRCBF','SINRHEU','nAntennasList','nUsersList');
[~,idxUsr] = min(abs(problem.nUsers - nUsersList));
[~,idxAnt] = min(abs(problem.N_Antennas - nAntennasList));
problem.SINR_CBF  = db2pow(SINRCBF(idxUsr,idxAnt));  % in linear scale
problem.SINR_LCMV = db2pow(SINRLCMV(idxUsr,idxAnt));  % in linear scale
problem.SINR_HEU  = db2pow(SINRHEU(idxUsr,idxAnt));  % in linear scale
% Configure Traffic
problem = f_configureTraffic(problem);  % Struct with configuration parameters
% Generate upper layer traffic
[traffic] = f_genDetTraffic(problem.class,problem.trafficType,problem.loadTraffic,problem.PLOT_DEBUG);
% Generate PHY-flows: Convert traffic (arrivals) into individual Flow for
% each user. Flows may overlap in time as the inter-arrival time may be
% less than Tslot
[flows,maxSlot] = f_arrivalToFlow(problem.Tslot,traffic,problem.class);
% Adequate Simulation time to the last packet arrival
if maxSlot~=0; problem.Tsym = maxSlot; end



% EOF