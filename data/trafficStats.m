% Clear workspace
clear; close all; clc;
% Load Data from UPC-Dataset
load('TrafficDataSetUPC.mat');
% Attach important iat and payload information
nTraffType = length(traffic);
for id = 1:nTraffType
    myTraffic = traffic{id};
    myTraffic_iat = [];  % distribution of inter-arrival times
    myTraffic_payloads = [];  % distribution of payloads in Bytes
    for idFlow = 1:myTraffic.numFlows
        myTraffic_iat = [myTraffic_iat ; myTraffic.times{idFlow}(2:end) - myTraffic.times{idFlow}(1:end-1)];  %#ok<AGROW>
        myTraffic_payloads = [myTraffic_payloads ; myTraffic.payload{idFlow}];  %#ok<AGROW>
    end
    traffic{id}.timesTot = myTraffic_iat;  %#ok
    traffic{id}.payloadTot = myTraffic_payloads;  %#ok
end
% Save with different name
save('TrafficDataSetUPC2.mat','traffic','appColorList','appNameList');