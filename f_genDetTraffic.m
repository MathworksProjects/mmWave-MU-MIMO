%% function [traffic] = f_genDetTraffic(traffic)
% This function generates traffic distributed deterministically in time.
% That is, the inter-arrival times between packets is constant accordingly
% to its deadline.
% INPUTS:
% - traffic: it is a vector of structs, each element containing the following
%            fields:
%            - arrivalTimes: Packet inter-arrival time in milliseconds
%            - deadline: Packet deadlines in millisecond
%            - numPkts: Number of packets to generate for that class
%            - Payload: Number of bits in each Network packet
% OUTPUTS:
% - traffic: Same as input but with two new fields:
%            - arrivalTimes: A vector of length numPkts containing the 
%                            packet arrival times.
%            - deadlines: A vector of length numPkts containing the 
%                         deadling that needs to be met for each packet in 
%                         time.
function [traffic] = f_genDetTraffic(traffic)
    Nclasses = length(traffic);  % Number of classes traffic (users)

    for id = 1:Nclasses
        traffic(id).arrivalTimes = zeros(traffic(id).numPkts,1);
        traffic(id).deadlines = zeros(traffic(id).numPkts,1);
        for pkt = 1:traffic(id).numPkts
            traffic(id).arrivalTimes(pkt,1) = pkt*traffic(id).iat;
            traffic(id).deadlines(pkt,1) = traffic(id).arrivalTimes(pkt) + traffic(id).deadline;
        end
    end
end
