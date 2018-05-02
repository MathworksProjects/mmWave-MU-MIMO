%% function [traffic] = f_genDetTraffic(traffic)
% This function generates traffic distributed deterministically in time.
% That is, the inter-arrival times between packets is constant accordingly
% to its deadline (for synthetic) and parse randomly select flows from the
% real DataSet (UPC). Also, distribute application across users depending
% upon the application presence, which can be either manually controlled
% (configuration script) or is based on the ratio numFlowApp/numFlowTotal
% from the DataSet.
% INPUTS:
% - traffic: it is a vector of structs, each element containing the following
%            fields:
%            - arrivalTimes: Packet inter-arrival time in milliseconds
%            - deadline: Packet deadlines in millisecond
%            - numPkts: Number of packets to generate for that class
%            - Payload: Number of bits in each Network packet
% - trafficType: Can either be 'synthetic' or 'dataSet'
% - DEBUG:       True for debugging purposes. False otherwise.
% OUTPUTS:
% - traffic: Same as input but with two new fields:
%            - arrivalTimes: A vector of length numPkts containing the 
%                            packet arrival times.
%            - deadlines: A vector of length numPkts containing the 
%                         deadling that needs to be met for each packet in 
%                         time.
% - maxTime: Maximum deadline time (used as simulation time or Tsym in the
% system in case trafficType is 'dataSet'.
function [traffic,maxTime] = f_genDetTraffic(traffic,trafficType,DEBUG)
    Nclasses = length(traffic);  % Number of classes traffic (users)
    if strcmp(trafficType,'synthetic')
        % Traffic is generated in a synthetic and deterministic way
        for id = 1:Nclasses
            payload = traffic(id).payload;
            traffic(id).payload = zeros(traffic(id).numPkts,1);
            traffic(id).arrivalTimes = zeros(traffic(id).numPkts,1);
            traffic(id).deadlines = zeros(traffic(id).numPkts,1);
            for pkt = 1:traffic(id).numPkts
                traffic(id).payload(pkt,1) = payload;
                traffic(id).arrivalTimes(pkt,1) = pkt*traffic(id).iat;
                traffic(id).deadlines(pkt,1) = traffic(id).arrivalTimes(pkt) + traffic(id).deadline;
            end
        end
        maxTime = 0;
    elseif strcmp(trafficType,'dataSet')
        % Traffic is generated in a synthetic and deterministic way
        if DEBUG
            ha = tight_subplot(Nclasses,1,[.05 .03],[.05 .01],[.1 .2]);
        end
        for id = 1:Nclasses
            traffic(id).payload = [];
            traffic(id).arrivalTimes = [];
            timeReference = 0;  % Initial time
            maxSpacing = 100;  % Spacing between two different flows
            flowSet = length(traffic(id).trafficFlows.times);
            numPktTot = 0;
            while(numPktTot<traffic(id).numPkts)
                if DEBUG
                    if strcmp(traffic(id).name,'Youtube')
                        flowID = 1803;  % Testing for Youtube
                    elseif strcmp(traffic(id).name,'Justin TV')
                        flowID = 1000;  % Testing for Justin TV
                    elseif strcmp(traffic(id).name,'Facebook')
                        flowID = 1116;  % Testing for Facebook
                    elseif strcmp(traffic(id).name,'Web Browsing')
                        flowID = 1500;  % Testing for Web Browsing
                    elseif strcmp(traffic(id).name,'Twitter')
                        flowID = 470;   % Testing for Twitter
                    end
                else
                    flowID = randi([1 flowSet]);
                end
                traffic(id).payload      = [traffic(id).payload ; ...
                                            traffic(id).trafficFlows.payload{flowID}];
                traffic(id).arrivalTimes = [traffic(id).arrivalTimes ; ...
                                            timeReference + traffic(id).trafficFlows.times{flowID}./1e-3];
                timeReference = traffic(id).arrivalTimes(end) + randi([1 maxSpacing]);
                numPktTot = numPktTot + length(traffic(id).trafficFlows.times{flowID});
            end
            traffic(id).deadlines = traffic(id).arrivalTimes + traffic(id).deadline;
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
           maxTime = max(maxTime,max(traffic(id).deadlines));
        end
    else
        fprintf('Traffic type incorrect - Ending execution\n');
        return;
    end
end