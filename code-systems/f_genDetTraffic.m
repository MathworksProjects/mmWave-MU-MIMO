function [traffic] = f_genDetTraffic(trafficClass,trafficType,loadTraffic,PLOT_DEBUG)
% f_genDetTraffic - generates traffic distributed deterministically in
% time. That is, the inter-arrival times between packets is constant
% accordingly to its deadline (for synthetic) and parse randomly select
% flows from the real DataSet (UPC). Also, distribute application across
% users depending upon the application presence, which can be either
% manually controlled (configuration script) or is based on the ratio
% numFlowApp/numFlowTotal from the DataSet.
%
% Syntax:  problem = f_configureTraffic(problem)
%
% Inputs:
%    traffic - vector of structs, each element has the following fields:
%            1. arrivalTimes: Packet inter-arrival time in milliseconds
%            2. deadline: Packet deadlines in millisecond
%            3. numPkts: Number of packets to generate for that class
%            4. Payload: Number of bits in each Network packet
%    trafficType -  Defines the origin of the data ('synthetic' or
%                   'dataSet')
%    loadTraffic -  Boolean that defines whether we generate the traffic
%                   from scatch based on DataSet or we load pre-stored
%                   values with same distribution
%    DEBUG - True for debugging purposes. False otherwise.
%
% Outputs:
%    traffic - Same as input but with two new fields:
%              1. arrivalTimes: A vector of length numPkts containing the 
%                               packet arrival times.
%              2. deadlines: A vector of length numPkts containing the 
%                            deadling that needs to be met for each packet 
%                            in time.
%
% Example: 
%       addpath('data','-end');
%       problem = o_read_input_problem('metaproblem_test.dat');
%       conf = o_read_config('config_test.dat');
%       problem = f_configureTraffic(problem);  % Struct with configuration parameters
%       [traffic,maxTime] = f_genDetTraffic(problem.class,problem.trafficType,false);
%
% Other m-files required: f_configureTraffic
% Subfunctions: none
% MAT-files required: data/metaproblem_test.dat,  data/config_test.dat
%
% See also: f_configure,  f_genDetTraffic

%------------- BEGIN CODE --------------

Nclasses = length(trafficClass);  % Number of classes traffic (users)
traffic = struct('name',string(),...
                 'payload',{},...
                 'arrivalTimes',[],...
                 'deadlines',[],...
                 'appColor',[],...
                 'numPkts',[]);
if strcmp(trafficType,'synthetic')
    % Traffic is generated in a synthetic and deterministic way
    for id = 1:Nclasses
        traffic(id).name = trafficClass(id).name;
        traffic(id).numPkts = trafficClass(id).numPkts;
        traffic(id).payload = zeros(trafficClass(id).numPkts,1);
        traffic(id).arrivalTimes = zeros(trafficClass(id).numPkts,1);
        traffic(id).deadlines = zeros(trafficClass(id).numPkts,1);
        for pkt = 1:trafficClass(id).numPkts
            traffic(id).payload(pkt,1) = trafficClass(id).payload;
            traffic(id).arrivalTimes(pkt,1) = pkt*trafficClass(id).iat;
            traffic(id).deadlines(pkt,1) = traffic(id).arrivalTimes(pkt) + trafficClass(id).deadline;
        end
    end
elseif strcmp(trafficType,'dataSet')
    % Traffic is generated following the empirical distribution in dataSet
    if PLOT_DEBUG
        ha = tight_subplot(Nclasses,1,[.05 .03],[.05 .01],[.1 .2]);
    end
    if loadTraffic
        % Load pre-generated traffic
        load('TrafficPreStored.mat','myIAT_Pre','myPayload_Pre');
    end
    for id = 1:Nclasses
        traffic(id).name = trafficClass(id).name;
        myTraffic = trafficClass(id).trafficFlows;
        trim = find(myTraffic.timesTot<var(myTraffic.timesTot)/1e4);
        if loadTraffic
            myIAT = myIAT_Pre(:,id);  %#ok
            myPayload = myPayload_Pre(:,id);  %#ok
        else
            % Create IAT and Payloads based on real distribution
            myIAT = emprand(myTraffic.timesTot(trim),trafficClass(id).numPkts,1);
            myPayload = emprand(myTraffic.payloadTot(trim)*8,trafficClass(id).numPkts,1);
        end
        traffic(id).arrivalTimes = ( mean(myIAT)*rand(1,1) + cumsum(myIAT) ).*1e+3;  % in ms
        traffic(id).payload = myPayload*8;  % Conversion to bits (DataSet in Bytes)
        traffic(id).deadlines = traffic(id).arrivalTimes + trafficClass(id).deadline;  % in ms
        traffic(id).appColor = trafficClass(id).appColor;
        traffic(id).numPkts = trafficClass(id).numPkts;
        if PLOT_DEBUG
            axes(ha(id));  %#ok
            stem(traffic(id).arrivalTimes,traffic(id).payload,'LineWidth',2,'Marker','none','Color',traffic(id).appColor);
            xlab = xlabel('time (ms)','FontSize',12);
            leg = legend(traffic(id).name);
            set(leg,'FontSize',12);
            % Manipulate y-axis
            yyaxis left
            ylabel('Bits','FontSize',12,'FontWeight','bold','Color','k');
            yyaxis right
            ylab = ylabel({'Packet Arrivals';'in Bits for';sprintf('User %d', id)},'FontSize',12,'FontWeight','bold','Color','k');
            set(ylab,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','left');
            ylim([0 300]); set(gca,'ytick',[]);
            p = get(xlab,'position'); p(2) = 0.9*p(2); set(xlab,'position',p);
            p1 = p.*[2.15 .5 1];  % Horizontal x Vertical x Depth
            set(xlab,'Position',p1);
            grid minor;
        end
    end
else
    error('Traffic type incorrect - Ending execution\n');
end


% EOF