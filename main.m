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
%% Load configuration
problem = o_read_input_problem('data/metaproblem_test.dat');
conf = o_read_config('data/config_test.dat');
%% Input parameters
problem = f_configuration(problem);  % Struct with configuration parameters
Tslot = 10;  % Time slot in milliseconds
Tsym  = 3e2;  % Total simulation time in milliseconds
%% Georgraphic distribution of users
[problem.thetaUsers, problem.phiUsers, problem.dUsers] = ...
    o_generate_positions(conf, problem.nUsers, problem.maxdUsers,problem.mindUsers);
%% First channel model (for the moment, including channel gains (alphas) !!!
[problem.thetaChannels, problem.phiChannels, problem.alphaChannels] = ...
                    o_generate_channels(conf,problem.nUsers,...
                    problem.maxnChannelPaths);
%% Handle traffic
[traffic] = f_genDetTraffic(problem.class);
% Convert traffic (arrivals) into individual Flow for each user. Flows may
% overlap in time as the inter-arrival time may be less than Tslot
[flows] = f_arrivalToFlow(Tslot,traffic);
% Copy initial flows for plotting purposes (to see progression of sim)
baseFlows = flows;
%% Main simulator
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
    [flows,selFlow] = f_distFlow(t,flows,Tslot,problem.FLAGagg);
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
            candTH = candTH(candTH~=0);
            problem.MinThr = candTH/problem.Bw;
            problem.MaxThr = problem.MinThr*Inf;
            % Call Heuristic method
%             [estSet,estTH] = f_heuristicsDummy(candSet,candTH);
            [sol_found,W,Cap] = f_heuristics(problem,conf,candSet);
            estTH = Cap*problem.Bw;
            % Decide whether to take the tentative TH or give it a
            % another round (This is Policy PLk)
            threshold = 0.7;  % Represents the ratio between the demanded 
                              % and the tentative achievable TH
            if ~any(estTH./candTH)<threshold
                % Evaluate PER
                TXbits = estTH.*Tslot.*1e-3;
                TXbitsTot = [TXbitsTot ; TXbits];               %#ok<AGROW>
                THTot = [THTot ; estTH];                        %#ok<AGROW>
                finalSet = f_PERtentative(candSet,[],[],[]);
                if ~isempty(finalSet); finalTH = estTH(candSet(finalSet));
                else;                  finalTH = [];
                end
                TXbitsTot(end,setdiff(candSet,finalSet)) = 0;
                THTot(end,setdiff(candSet,finalSet)) = 0;
                % Update remaining bits to be sent upon tx success
                flows = f_updateFlow(t,flows,selFlow,finalSet,finalTH,candSet,Tslot,problem.DEBUG);
                lastSelFlow = selFlow;
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
    figure(1); subplot(problem.nUsers,1,id); hold on;
    bar((1:1:Tsym-1),TXbitsTot(:,id),'LineWidth',2,'EdgeColor','none','FaceColor','red');
    for fl = 1:lastSelFlow
        slots = baseFlows(id).slots{fl};
        reqTXbits = baseFlows(id).TH(fl).*Tslot.*1e-3.*ones(1,length(slots));
        bar(slots,reqTXbits,'LineWidth',3,'EdgeColor','none','FaceColor','blue','FaceAlpha',0.4);
    end
%     title('Number of Bits transmitted','FontSize',12);
    xlabel('Slot index','FontSize',12);
    ylabel('Bits transmitted','FontSize',12);
    lg = legend('TX Bits over the channel','Baseline required Bits to be TX');
    set(lg,'FontSize',12,'Location','Northeast');
    grid minor;

    figure(2); subplot(problem.nUsers,1,id); hold on;
    bar((1:1:Tsym-1),THTot(:,id).*1e-6,'LineWidth',3,'EdgeColor','none','FaceColor','red');
    for fl = 1:lastSelFlow
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
% TODO: Implement a more realistic traffic model (Apply Stratis' advice)
% TODO: Come up with a closed for of the priority calculation, more complex
%       than just the inverse of the time_to_deadline.
% TODO: Implement the non_aggregate method inside f_distFlow()
% TODO: Combine Santi's, Zhengnan's and main code. Debugg and improve!