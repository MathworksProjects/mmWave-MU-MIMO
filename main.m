%% Clear workspace
clear; clc; close all;
%% Input parameters
conf  = f_configuration;  % Struct with configuration parameters
Tslot = 10;  % Time slot in milliseconds
Tsym  = 1e3;  % Total simulation time in milliseconds
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
    fprintf('--- Sim time %d\n',t);
    % Distribute Flow accross users. Either we aggregate or disaggregate
    % overlapping flows in the current slot. Select the current flow for
    % each user
    [flows,selFlow] = f_distFlow(t,flows,conf.FLAGagg);
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
                if ~isempty(finalSet); finalTH = estTH(finalSet==candSet);
                else;                  finalTH = [];
                end
                % Update remaining bits to be sent upon tx success
                flows = f_updateFlow(flows,selFlow,finalSet,finalTH,candSet,Tslot);
                % Exit the for loop - we have served in this time slot
                break;
            end
        end
    end
    % Increment variable event in DES
    t = t + 1;
end

%% REPORT
% (TODO) Debugg f_distFlow Line 55 - Some issue with traffic aggregation
%% TODO LIST
% TODO: Debug 
% TODO: Syncronize it with repository on github
% TODO: Come up with a closed for of the priority calculation, more complex
%       than just the inverse of the time_to_deadline.
% TODO: Implement the non_aggregate method inside f_distFlow()
% TODO: Implement a more realistic traffic model (Talk to Stratis)
