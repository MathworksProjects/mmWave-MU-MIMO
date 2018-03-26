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
%% Input parameters
conf  = f_configuration;  % Struct with configuration parameters
Tslot = 10;  % Time slot in milliseconds
Tsym  = 1e2;  % Total simulation time in milliseconds
%% Handle traffic
[traffic] = f_genDetTraffic(conf.class);
% Convert traffic (arrivals) into individual Flow for each user. Flows may
% overlap in time as the inter-arrival time may be less than Tslot
[flows] = f_arrivalToFlow(Tslot,traffic);
%% Main simulator
% Represent the time (in slot ID) throughout the execution. It is the even
% in our DES
t = 0;
while(t<Tsym)
%     fprintf('** SLOT ID: %d\n',t);
    % Distribute Flow accross users. Either we aggregate or disaggregate
    % overlapping flows in the current slot. Select the current flow for
    % each user
    [flows,selFlow] = f_distFlow(t,flows,Tslot,conf.FLAGagg);
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
            % Call Heuristic method
            [estSet,estTH] = f_heuristicsDummy(candSet,candTH);
            % Decide whether to take the tentative TH or give it a
            % another round (This is Policy PLk)
            threshold = 0.7;  % Represents the ratio between the demanded 
                              % and the tentative achievable TH
            if ~any(estTH./candTH)<threshold
                % Evaluate PER
                finalSet = f_PERtentative(candSet);
                if ~isempty(finalSet); finalTH = estTH(candSet(finalSet));
                else;                  finalTH = [];
                end
                % Update remaining bits to be sent upon tx success
                flows = f_updateFlow(t,flows,selFlow,finalSet,finalTH,candSet,Tslot,conf.DEBUG);
                % Exit the for loop - we have served in this time slot,
                % there's no way back even though some pkts didn't make it
                break;
            end
        end
    end
    % Increment variable event in DES
    t = t + 1;
end

%% REPORT
[~,~] = f_generateReport(flows);

%% TODO LIST
% TODO: Implement a more realistic traffic model (Apply Stratis' advice)
% TODO: Come up with a closed for of the priority calculation, more complex
%       than just the inverse of the time_to_deadline.
% TODO: Implement the non_aggregate method inside f_distFlow()
% TODO: Combine Santi's, Zhengnan's and main code. Debugg and improve!