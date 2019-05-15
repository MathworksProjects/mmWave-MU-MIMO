function [flows,maxSlot] = f_arrivalToFlow(Tslot,traffic,trafficClass)
% F_ARRIVALTOFLOW - % This function converts the packets into individual
% flows. Packets arrive at a discret time with an exact amount of bits. A
% flow consist of a desirable Throughput accross time (slots).
% 
% The function returns 'flows', which is an array of structs of length
% equal to the number of users. For each user, each flow (belonging to each
% packet) is characterized by the amount of bits that needs to be
% delivered. The amount of bits are distributed uniformly across the slots
% until reaching the deadline. Thus, the variable flows contains four
% features:
% - slots:     For each flow, the slots across it.
% - TH:        The average throughput demanded for that flow over the slots
%              indicated in the slots field.
% - remaining: We set the remaining field of the flow to be the value of
%              the Payload.
% - deadlines: The deadline of the actual packet in slot ID. 
%              This is used in the future to set priorities.
% - failed:    Mark whether the flow has failed to be served before the
%              deadline (0 or 1).
% - success:   Mark whether the flow has been succesfully served before the
%              deadline (0 or 1).
% - maxSlot:   Maximum deadline slot used as simulation time or Tsym in
%              the system in case trafficType is 'dataSet'.
%
% Syntax:  [flows,maxSlot] = f_arrivalToFlow(Tslot,traffic,trafficClass)
%
% Inputs:
%    Tslot - struct containint configuration in data/metaproblem_test.dat
%    traffic - struct containint configuration in data/config_test.dat
%    trafficClass - Vector containing the users IDs
%
% Outputs:
%    flows - flag indicating whether the solution was found
%    maxSlot - beamforming weights
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: none
%
% See also: 
%
% Example:  flows(1).slots(4) = [11 12 13 14];  % Flow 4 from user 1 occupies 
%                                               % slots 11, 12, 13 and 14 
%                                               % with an average traffic of flows(1).TH.

%------------- BEGIN CODE --------------
Nusers = length(traffic);  % Total number of users. Consider 1 
                           % application per user.

% Prealocate memory for structure users. Initialize the two fields to
% empty since the length for slots is variable for each user (number of
% packets). i.e. users(U).slots{5} is a list containing the slots where
% the packet (flow) 5 from user U is present.
flows = struct('slots',{},...
               'TH',[],...
               'remaining',[],...
               'remainingPerSlot',[],...
               'deadlines',[],...
               'failed',[],...
               'success',[],...
               'numPkts',[]);

% Convert packets into demanded throughput following a uniform
% distribution accros the interval <tArrival,tDeadline>
maxSlot = 0;
for id = 1:Nusers
    arrivals = traffic(id).arrivalTimes;
    arrivals(arrivals==0) = eps;  % Infinitesimal time correction
    deadlines = traffic(id).deadlines;
    Npkt = traffic(id).numPkts;
    % Pre-allocate memory for a faster execution
    valor.slots = cell(Npkt,1);
    valor.TH = zeros(Npkt,1);
    valor.remaining = zeros(Npkt,1);
    valor.deadlines = zeros(Npkt,1);
    valor.numPkts = Npkt;
    for pkt = 1:Npkt
        % The first slot is the right next available slot after the
        % moment of arrival (ceil). The last slot is the slots that
        % comprises the deadline time (floor)
        valor.slots{pkt} = (ceil(arrivals(pkt)/Tslot) : 1 : ...
                            floor(deadlines(pkt)/Tslot));
        nSlots = length(valor.slots{pkt});
        % The requested Throughput is computed in bits per second (bps)
        valor.TH(pkt) = traffic(id).payload(pkt) / (nSlots*Tslot*1e-3);  % Conversion from bytes to bitss
        valor.remaining(pkt) = traffic(id).payload(pkt);  % Conversion from bytes to bits
        valor.remainingPerSlot(pkt) = valor.remaining(pkt) / nSlots;
%             valor.deadlines(pkt) = ceil(traffic(id).deadlines(pkt)/Tslot);
        valor.deadlines(pkt) = floor(traffic(id).deadlines(pkt)/Tslot);
        valor.failed(pkt) = 1;
        valor.success(pkt) = 0;
    end
    flows(id) = valor;
    maxSlot = max(maxSlot,max(flows(id).deadlines)+5);
end



% EOF