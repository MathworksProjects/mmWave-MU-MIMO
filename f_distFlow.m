%% function [flows,selectedFlows] = f_distFlow(t,flows,aggregate)
% This function is intended to be called at the very beginning of the
% iterative loop in the simulation. It controls the flows of the users. If
% more than one flow is present in the current time slot 't', then the
% flows are either aggregated or disaggregated. The aggregation is
% controlled by the flag 'aggregate'. The flows are updated and outputed
% accordingly. The function also returns the flow_ID for every user. If the
% ID is 0, no flow is selected (indices start in 1).
%   INPUTS:
%   - t <----------> The current slot.
%   - flows: <-----> Information about the flow. For instance, the slots
%                    over which it is intended to transmit, the remaining 
%                    amount of bits, the expected throughput at a certain 
%                    time slot, etcetera.
%   - aggregate: <-> Boolean Flag (True or False) that controls the flow
%                    aggregation policy.
%   OUTPUTS:
%   - flows: <-----> Flows are updated accordingly in case there is an overlap in
%                    time.
%   - selFlows: <--> Contains the flow IDs for each user. If a user do not
%                    have a flow available at instant 't', the value is 0. 
%                    The ID is used to map the requirements from variable 
%                    'flow' further in the system.
function [flows,selectedFlows] = f_distFlow(t,flows,Tslot,aggregate)
    Nusers = length(flows);  % Total number of users
    selectedFlows = zeros(Nusers,1);
    for id = 1:Nusers
        numPkt = length(flows(id).slots);  % Maximum number of packets to iterate over
        for pkt = 1:numPkt
            % Initialize control variables
            idx2 = false;
            % idx0 controls weather a flow is present in the slot
            idx0 = find(ismember(flows(id).slots{pkt},t)~=0,1);
            % Only evaluate flow aggregation if we are not the last flow to
            % avoid array indexation out of available range (pkt+1)
            if pkt<numPkt
                % idx1 controls overlapping between flows
                idx1 = find(ismember(flows(id).slots{pkt},flows(id).slots{pkt+1})~=0,1);
                % idx2 controls if the overlapping happens within current slot
                idx2 = isequal(flows(id).slots{pkt}(idx1),t);
            end
            if idx2 && aggregate
                % We DO aggregate traffic. This option increases the efficiency of the
                % traffic distribution. We follow a uniform traffic distribution
                % policy across the deadlines. Traffic is evenly distributed
                % between two periods:
                % Delta = 2nd deadline - current
                Delta = flows(id).deadlines(pkt+1) - t + 1;
                % Delta1 = 1nd deadline - current
                Delta1 = flows(id).deadlines(pkt) - t + 1;
                % Delta2 = 2nd deadline - 1nd deadline
                Delta2 = flows(id).deadlines(pkt+1) - flows(id).deadlines(pkt);
                % Define proportional variables
                X = flows(id).remaining(pkt);
                Y = flows(id).remaining(pkt+1);
                alpha = Delta1/Delta;
                beta = Delta2/Delta;
                % Redefine remaining bits to be transmitted in each flow
                flows(id).remaining(pkt)   = alpha * (X+Y);
                flows(id).remaining(pkt+1) =  beta * (X+Y);
                % Redefine the TH
                flows(id).TH(pkt) = flows(id).remaining(pkt) / (Delta1*Tslot*1e-3);
                flows(id).TH(pkt+1) = flows(id).remaining(pkt+1) / (Delta2*Tslot*1e-3);
                % Cut the number of slots for the second overlapping flow
                flows(id).slots{pkt+1} = (max(flows(id).slots{pkt}) + 1 ...
                                          : 1 : max(flows(id).slots{pkt+1}));
                selectedFlows(id) = pkt;
            elseif idx2 && ~aggregate
                % (TODO)
                % We DO NOT aggregate traffic, meaning that we disaggregate any
                % overlapping flow in time (slot). In other words, we wait until we
                % finish the transmission of one flow before start transmitting the
                % next one. While we decrease the load on the antenna selection
                % side (less throughput) on the first interval (1nd deadline -
                % current), we increase substantially the demanded throughput in
                % the second one (2nd deadline - 1nd deadline).
            elseif ~isempty(idx0)
                % There is only one flow within this time slot. Its ID is
                % determined by the flow (pkt) number
                selectedFlows(id) = pkt;
                break;
            end
        end
    end
end
