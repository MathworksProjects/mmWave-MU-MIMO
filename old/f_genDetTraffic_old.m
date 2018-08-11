function [traffic,maxTime] = f_genDetTraffic(trafficClass,trafficType,DEBUG)
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
%    DEBUG - True for debugging purposes. False otherwise.
%
% Outputs:
%    traffic - Same as input but with two new fields:
%              1. arrivalTimes: A vector of length numPkts containing the 
%                               packet arrival times.
%              2. deadlines: A vector of length numPkts containing the 
%                            deadling that needs to be met for each packet 
%                            in time.
%    maxTime -  Maximum deadline time (used as simulation time or Tsym in
%               the system in case trafficType is 'dataSet').
%
% Example: 
%       problem = o_read_input_problem('data/metaproblem_test.dat');
%       conf = o_read_config('data/config_test.dat');
%       problem = f_configureTraffic(problem);  % Struct with configuration parameters
%       [traffic,maxTime] = f_genDetTraffic(problem.class,problem.trafficType,false);
%
% Other m-files required: f_configureTraffic
% Subfunctions: none
% MAT-files required: data/metaproblem_test.dat,  data/config_test.dat
%
% See also: f_configure,  f_genDetTraffic
% 
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
    maxTime = 0;
elseif strcmp(trafficType,'dataSet')
    % Traffic is generated in a synthetic and deterministic way
    if DEBUG
        ha = tight_subplot(Nclasses,1,[.05 .03],[.05 .01],[.1 .2]);
    end
    for id = 1:Nclasses
        % v1 traffic generation
%         traffic(id).payload = [];
%         traffic(id).arrivalTimes = [];
%         timeReference = 0;  % Initial time
%         maxSpacing = 100;  % Spacing between two different flows
%         flowSet = length(trafficClass(id).trafficFlows.times);
%         numPktTot = 0;
%         while(numPktTot<trafficClass(id).numPkts)
%             if DEBUG
%                 if strcmp(trafficClass(id).name,'Youtube')
%                     flowID = 1803;  % Testing for Youtube
%                 elseif strcmp(trafficClass(id).name,'Justin TV')
%                     flowID = 1000;  % Testing for Justin TV
%                 elseif strcmp(trafficClass(id).name,'Facebook')
%                     flowID = 1116;  % Testing for Facebook
%                 elseif strcmp(trafficClass(id).name,'Web Browsing')
%                     flowID = 1500;  % Testing for Web Browsing
%                 elseif strcmp(trafficClass(id).name,'Twitter')
%                     flowID = 470;   % Testing for Twitter
%                 end
%             else
%                 flowID = randi([1 flowSet]);
%             end
%             traffic(id).payload      = [traffic(id).payload ; ...
%                                         trafficClass(id).trafficFlows.payload{flowID}];
%             traffic(id).arrivalTimes = [traffic(id).arrivalTimes ; ...
%                                         timeReference + trafficClass(id).trafficFlows.times{flowID}./1e-3];
%             timeReference = traffic(id).arrivalTimes(end) + randi([1 maxSpacing]);
%             numPktTot = numPktTot + length(trafficClass(id).trafficFlows.times{flowID});
%         end
%         traffic(id).deadlines = traffic(id).arrivalTimes + trafficClass(id).deadline;
        % v2 traffic generation
        traffic(id).name = trafficClass(id).name;
        myTraffic = trafficClass(id).trafficFlows;
        myIAT = emprand(myTraffic.timesTot,trafficClass(id).numPkts,1);
        traffic(id).arrivalTimes = mean(myIAT)*rand(1,1) + cumsum(myIAT);
        traffic(id).payload = emprand(myTraffic.payloadTot,trafficClass(id).numPkts,1);
        traffic(id).deadlines = traffic(id).arrivalTimes + trafficClass(id).deadline;
        traffic(id).appColor = trafficClass(id).appColor;
        if DEBUG
            axes(ha(id));
            stem(traffic(id).arrivalTimes,traffic(id).payload,'LineWidth',2,'Marker','none','Color',traffic(id).appColor);
            xlab = xlabel('time (s)','FontSize',12);
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
    maxTime = 0;
    for id = 1:Nclasses
       maxTime = max(maxTime,max(traffic(id).deadlines));  % in ms
    end
else
    fprintf('Traffic type incorrect - Ending execution\n');
    return;
end


% EOF