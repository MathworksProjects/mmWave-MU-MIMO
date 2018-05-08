function [flows,CapTot,TXbitsTot,THTot,lastSlotSim,lastSelFlow,varargout] = main(varargin)
%MAIN - This is the main runnable in the Project. The code is built as a Discrete
%Event Simulator (DES) that emulates a 5G Base Station (BS) operating in
%the millimeter wave band. The code iterates over the following steps: 
%
%  1. Generate 5G traffic closely ressembling the real traffic based on
%     information extracted from Data Bases. 
%  2. Generate traffic flow (in Throughput) per time slot and per user. 
%  3. Allocate users per time slot based on their application demands
%     packet deadline to meet required latency) and the throughput 
%     necessary based on step 2.
%  4. Configure subarrays (antenna sets) at the BS and allocate them to the
%     users to transmit concurrently. We apply an evolutionary algorithm known
%     as Particle-Swarm-Optimization (PSO) to obtain the optimal solution.
%  5. Evaluate the PER using standard-compliant frame and mmWave channel.
%  6. Update traffic flow and iterate back to step 1.
%
% Syntax:  [flows,TXbitsTot,THTot,lastSlotSim,lastSelFlow,varargout] = main([])
%          [flows,TXbitsTot,THTot,lastSlotSim,lastSelFlow,varargout] = ...
%                         main(runnable,conf,problem,flows)
%
% Inputs:
%    conf [optional] - Description
%    problem [optional] - Description
%    flows [optional] - Description
%
% Outputs:
%    flows - Description
%    CapTot - Description
%    TXbitsTot - Description
%    THTot - Description
%    lastSlotSim - Description
%    lastSelFlow - Description
%    baseFlow [optional] - Description
%
% Example: 
%    [flows,TXbitsTot,THTot,lastSlotSim,lastSelFlow,varargout] = main([])
%    main_plotting(problem,TXbitsTot,THTot,baseFlows,lastSelFlow);
%
% Other m-files required: Requires most of the o_*, s_* and f_* functions
% Subfunctions: Calls most of the o_*, s_* and f_* functions
% MAT-files required: ToDo
%
% See also: main_runnable
%
%------------- BEGIN CODE --------------

% Check correct number of input arguments. If none, then set-up default
% parameters from configuration functions
if (nargin==3)
    fprintf('Main: Parameters extracted from runnable. Running main...\n');
    % Check inputted parameters have the correct format and assign values
    % in current function if every requirement is met
    if ~isstruct(varargin{1}) || any( structfun(@isempty, varargin{1}) )
        error('ERROR: 1st input parameter should be a non-empty struct ' + ...
            'contaitning the configuration parameters of heuristics\n');
    elseif ~isstruct(varargin{2}) || any( structfun(@isempty, varargin{1}) )
        error('ERROR: 2nd input parameter should be a non-empty struct ' + ...
            'contaitning the configuration parameters of heuristics\n');
    elseif ~isstruct(varargin{3})
        error('ERROR: 3rd input parameter should be a non-empty struct\n');
    elseif isstruct(varargin{3})
        for k = 1:length(varargin{3})
            if any( structfun(@isempty, varargin{3}(k)) )
                error('ERROR: 3rd input parameter should be a non-empty struct ' + ...
                    'contaitning the PHY-flow information\n');
            end
        end
    end
    % Rename input variables
    conf = varargin{1};
    problem = varargin{2};
    flows = varargin{3};
elseif (nargin==0)
    % No parameters were inputted to the main
    %% Clear workspace
    clear; clc; close all;
    addpath('utilities','-end');  % Add utilities folder at the end of search path
    %% Load configuration
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    %% Input parameters
    problem = f_configuration(problem);  % Struct with configuration parameters
    %% Geographic distribution of users
    [problem.thetaUsers, problem.phiUsers, problem.dUsers] = ...
        o_generate_positions(conf, problem.nUsers, problem.maxdUsers,problem.mindUsers);
    %% Generate channels per user
    if conf.Use5GChannel
        % 5G 3GPP ETSI TR 38.901 compliant channel model
        [~, problem.thetaChannels, problem.phiChannels, problem.alphaChannels] = ...
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
    %% Handle traffic
    [traffic,maxTime] = f_genDetTraffic(problem.class,problem.trafficType,problem.DEBUG);
    if maxTime~=0
        problem.Tsym = maxTime;  % Adequate Simulation time to the last packet arrival
    end
    % Convert traffic (arrivals) into individual Flow for each user. Flows may
    % overlap in time as the inter-arrival time may be less than Tslot
    [flows] = f_arrivalToFlow(problem.Tslot,traffic);
    % Copy initial flows for plotting purposes (to see progression of sim)
    varargout{1} = flows;  % baseFlows: for printing purposes at the end of execution
    fprintf('Main: Parameters ovewritten in Main. Running main...\n');
else
    error('ERROR: The number of input parameters should be either 0 or 3\n' + ...
    'Please, check the parameter formatting on the description of this function.\n');
end

%% Main simulator - MAIN SECTION
Tsym = problem.Tsym;
Tslot = problem.Tslot;
% Represent the time (in slot ID) throughout the execution. It is the even
% in our DES
t = 1;
TXbitsTot = [];
THTot = [];
CapTot = [];
lastSelFlow = zeros(problem.nUsers,1);  % For printing purposes at the end of execution
while(t<Tsym)
%     fprintf('** SLOT ID: %d\n',t);
    % Distribute Flow accross users. Either we aggregate or disaggregate
    % overlapping flows in the current slot. Select the current flow for
    % each user
    [flows,selFlow] = f_distFlow(t,flows,Tslot,problem.FLAGagg,problem.DEBUG);
    % Compute priorities. As of now, priorities are computed as the inverse
    % of the time_to_deadline (for simplicity)
    if any(selFlow)~=0
        [combIds,combTH] = f_candidateSet(t,flows,selFlow);
        % Iterate over t
        for k = 1:length(combIds)
            % Retrieve user set to be scheduled and remove 0s from the set
            candSet = combIds(k,:);
            candSet = candSet(candSet~=0);
            candTH = combTH(k,:);
%             candTH = candTH(candTH~=0);
            %% SECTION HEURISTICS METHOD
            % Heuristics - Preprocessing
            problem.MaxObjF = Inf(1,length(candSet));
            problem.MinObjF = candTH/problem.Bw;
            if conf.MinObjFIsSNR;     problem.MinObjF = 2.^problem.MinObjF - 1;
            end
            % Heuristics - Call
            if problem.heuristicsDummy && ~isempty(candSet)
                % Dummy heuristics
                [estObj] = f_heuristicsDummy(problem.MinObjF,conf.MinObjFIsSNR,problem.MCSPER.snrRange);
            elseif ~problem.heuristicsDummy && ~isempty(candSet)
                % Real Heuristics
                [~,~,~,estObj] = f_heuristics(problem,conf,candSet);
            end
            % Heuristics - Post Processing
            if conf.MinObjFIsSNR;     estTH   = log2(estObj+1)*problem.Bw;  % in bps
                                      Caps = zeros(1,problem.nUsers);
                                      Caps(1,candSet) = log2(estObj+1);
                                      CapTot  = [CapTot ; Caps];  %#ok<AGROW> % in bps/Hz
                                      SNRList = 10*log10(estObj);  % in dB
            else;                     estTH   = estObj*problem.Bw;  % in bps
                                      Caps = zeros(1,problem.nUsers);
                                      Caps(1,candSet) = estObj;
                                      CapTot  = [CapTot ; Caps];  %#ok<AGROW> % in bps/Hz
                                      SNRList = 10*log10(2.^(estTH/problem.Bw) - 1);  % in dB
            end
            %% SECTION POST HEURISTICS
            % Decide whether to take the tentative TH or give it a
            % another round (This is Policy PLk)
            threshold = 0.7;  % Represents the ratio between the demanded 
                              % and the tentative achievable TH
            if ~any(estTH./candTH)<threshold
                % Select MCS for estimated SNR
                [~,PER] = f_selectMCS(candSet,SNRList,problem.targetPER,problem.MCSPER,problem.DEBUG);
                % Compute bits that can be transmitted and map it with the 
                % bits remaining to be transmitted
                TXbits = zeros(1,problem.nUsers);  % Reality - bits
                THiter = zeros(1,problem.nUsers);  % Reality - throughput
                for id = candSet
                    estTXbits = estTH(id==candSet).*Tslot.*1e-3;
                    if estTXbits > flows(id).remaining(selFlow(id))  %#ok  % estTXbits is always scalar
                        % Bits transmitted in slot
                        TXbits(1,id) = flows(id).remaining(selFlow(id));
                        % Throughput achieved in slot
                        THiter(id) = TXbits(1,id)./(Tslot.*1e-3);
                    else
                        % Bits transmitted in slot
                        TXbits(1,id) = estTXbits;
                        % Throughput achieved in slot
                        THiter(id) = estTH(id==candSet);
                    end
                end
                % Append bits transmitted and throughput achieved in slot. 
                TXbitsTot = [TXbitsTot ; TXbits];               %#ok<AGROW>
                % Append throughput achieved in slot
                THTot = [THTot ; THiter];                       %#ok<AGROW>
                % Evaluate PER
                finalSet = f_PERtentative(candSet,PER);
                if ~isempty(finalSet); finalTH = THiter(finalSet);
                else;                  finalTH = [];
                end
                TXbitsTot(end,setdiff(candSet,finalSet)) = 0;
                THTot(end,setdiff(candSet,finalSet)) = 0;
                % Update remaining bits to be sent upon tx success
                flows = f_updateFlow(t,flows,selFlow,finalSet,finalTH,candSet,Tslot,problem.DEBUG);
                lastSelFlow(selFlow~=0) = selFlow(selFlow~=0);
                % Exit the for loop - we have served in this time slot,
                % there's no way back even though some pkts didn't make it
                break;
            else
                TXbitsTot = [TXbitsTot ; zeros(1,problem.nUsers)];   %#ok<AGROW>
                THTot = [THTot ; zeros(1,problem.nUsers)];           %#ok<AGROW>
            end
        end
    else
        TXbitsTot = [TXbitsTot ; zeros(1,problem.nUsers)];   %#ok<AGROW>
        THTot = [THTot ; zeros(1,problem.nUsers)];           %#ok<AGROW>
    end
    % Increment variable event in DES
    t = t + 1;
end

lastSlotSim = t - 1;

%EOF