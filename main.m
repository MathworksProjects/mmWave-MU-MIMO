%% ======================== IMPORTANT INFORMATION ====================== %%
% This is the main runnable in the Project. The code is built as a Discrete
% Event Simulator (DES) that emulates a 5G Base Station (BS) operating in
% the millimeter wave band. The code iterates over the following steps: 
%
% 1) Generate 5G traffic closely ressembling the real traffic based on
% information extracted from Data Bases. 
% 2) Generate traffic flow (in Throughput) per time slot and per user. 
% 3) Allocate users per time slot based on their application demands
% (packet deadline to meet required latency) and the throughput necesaary
% based on step 2.
% 4) Configure subarrays (antenna sets) at the BS and allocate them to the
% users to transmit concurrently. We apply an evolutionary algorithm known
% as Particle-Swarm-Optimization (PSO) to obtain the optimal solution.
% 5) Evaluate the PER using standard-compliant frame and mmWave channel.
% 6) Update traffic flow and iterate back to step 1.
%% ===================================================================== %%
%% Clear workspace
clear; clc; close all;
addpath('utilities','-end');  % Add utilities folder at the end of search path
%% Load configuration
problem = o_read_input_problem('data/metaproblem_test.dat');
conf = o_read_config('data/config_test.dat');
%% Input parameters
problem = f_configuration(problem);  % Struct with configuration parameters
Tslot = 10;  % Time slot in milliseconds
Tsym  = 3e2;  % Total simulation time in milliseconds
%% Geographic distribution of users
[problem.thetaUsers, problem.phiUsers, problem.dUsers] = ...
    o_generate_positions(conf, problem.nUsers, problem.maxdUsers,problem.mindUsers);
%% Generate channels per user
if conf.Use5GChannel
    % 5G 3GPP ETSI TR 38.901 compliant channel model
    [FullChannels, problem.thetaChannels, problem.phiChannels, problem.alphaChannels] = ...
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
    Tsym = maxTime;  % Adequate Simulation time to the last packet arrival
end
% Convert traffic (arrivals) into individual Flow for each user. Flows may
% overlap in time as the inter-arrival time may be less than Tslot
[flows] = f_arrivalToFlow(Tslot,traffic);
% Copy initial flows for plotting purposes (to see progression of sim)
baseFlows = flows;  % For printing purposes at the end of execution
lastSelFlow = zeros(problem.nUsers,1);  % For printing purposes at the end of execution
%% Main simulator - MAIN SECTION
% Represent the time (in slot ID) throughout the execution. It is the even
% in our DES
t = 1;
TXbitsTot = [];
THTot = [];
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
                [sol_found,W,array_Handle,estObj] = f_heuristics(problem,conf,candSet);
            end
            % Heuristics - Post Processing
            if conf.MinObjFIsSNR;     estTH   = log2(estObj+1)*problem.Bw;  % in bps/Hz
                                      SNRList = 10*log10(estObj);  % in dB
            else;                     estTH   = estObj*problem.Bw;  % in bps/Hz
                                      SNRList = 10*log10(2.^(estTH/problem.Bw) - 1);  % in dB
            end
            %% SECTION POST HEURISTICS
            % Decide whether to take the tentative TH or give it a
            % another round (This is Policy PLk)
            threshold = 0.7;  % Represents the ratio between the demanded 
                              % and the tentative achievable TH
            if ~any(estTH./candTH)<threshold
                % Select MCS for estimated SNR
                [MCS,estPER] = f_selectMCS(candSet,SNRList,problem.targetPER,problem.MCSPER,problem.DEBUG);
                % Compute bits that can be transmitted and map it with the 
                % bits remaining to be transmitted
                estTXbits = zeros(1,problem.nUsers);  % Possibility - bits
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
                finalSet = f_PERtentative(candSet,[],[],[]);
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

%% REPORT
[~,~] = f_generateReport(flows);

%% PLOTTING
for id = 1:problem.nUsers
    figure(1); subplot(problem.nUsers,2,2*id - 1); hold on;
    bar((1:1:Tsym-1),TXbitsTot(:,id),'LineWidth',2,'EdgeColor','none','FaceColor','red');
    bitsToTxTot = zeros(1,max(baseFlows(id).slots{lastSelFlow(id)}));
    for fl = 1:lastSelFlow(id)
        slots = baseFlows(id).slots{fl};
        reqTXbits = baseFlows(id).TH(fl).*Tslot.*1e-3.*ones(1,length(slots));
        bitsToTxTot(slots) = bitsToTxTot(slots) + reqTXbits;
    end
    lastSlot = max(slots);
    bar((1:1:lastSlot),bitsToTxTot,'LineWidth',3,'EdgeColor','none','FaceColor','blue','FaceAlpha',0.4);
%     title('Number of Bits transmitted','FontSize',12);
    xlabel('Slot index','FontSize',12);
    ylabel('Bits TX','FontSize',12);
    lg = legend('TX Bits over the channel','Baseline required Bits to be TX');
    set(lg,'FontSize',12,'Location','Northeast');
    grid minor;
    
    figure(1); subplot(problem.nUsers,2,2*id); hold on;
    bitsToTxTot = zeros(1,max(baseFlows(id).slots{lastSelFlow(id)}));
    for fl = 1:lastSelFlow(id)
        slots = baseFlows(id).slots{fl};
        reqTXbits = baseFlows(id).TH(fl).*Tslot.*1e-3.*ones(1,length(slots));
        bitsToTxTot(slots) = bitsToTxTot(slots) + reqTXbits;
        bar(slots ,bitsToTxTot(slots),'LineWidth',3,'EdgeColor','none','FaceAlpha',0.8);
    end
%     lastSlot = max(slots);
%     bar((1:1:lastSlot),reqTXbits,'LineWidth',3,'EdgeColor','none','FaceAlpha',0.1);
    xlabel('Slot index','FontSize',12);
    ylabel('Bits to be TX','FontSize',12);
    lg = legend('Agg. bits flow 1','Agg. bits flow 1+2','Agg. bits flow 1+2+3');
    set(lg,'FontSize',12,'Location','Northeast');
    grid minor;

    figure(2); subplot(problem.nUsers,1,id); hold on;
    bar((1:1:Tsym-1),THTot(:,id).*1e-6,'LineWidth',3,'EdgeColor','none','FaceColor','red');
    for fl = 1:lastSelFlow(id)
        slots = baseFlows(id).slots{fl};
        reqTH = baseFlows(id).TH(fl).*1e-6.*ones(1,length(slots));
        bar(slots,reqTH,'LineWidth',2,'EdgeColor','none','FaceColor','blue','FaceAlpha',0.4);
    end
%     title('Throughput required per slot','FontSize',12);
    xlabel('Slot index','FontSize',12);
    ylabel('Throughput (Mbps)','FontSize',12);
    lg = legend('Achieved','Initially Required');
    set(lg,'FontSize',12,'Location','Northeast');
    grid minor;
    a = get(gcf,'Position');
    set(gcf,'Position',[a(1) a(2) 560 168])
end

%% TODO LIST
% TODO: Come up with a closed for of the priority calculation, more complex
%       than just the inverse of the time_to_deadline.
% TODO: Implement the non_aggregate method inside f_distFlow()
% TODO: Combine Santi's, Zhengnan's and main code. Debugg and improve!