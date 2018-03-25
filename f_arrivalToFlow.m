%% function [flows] = f_arrivalToFlow(Tslot,traffic)
% This function converts the packets into individual flows. Packets arrive
% at a discret time with an exact amount of bits. A flow consist of a
% desirable Throughput accross time (slots).
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
%
% Example: flows(1).slots(4) = [11 12 13 14]. This means that flow 4 from
% user 1 occupies slots 11, 12, 13 and 14 with an average traffic of
% flows(1).TH.
function [flows] = f_arrivalToFlow(Tslot,traffic)
    Nusers = length(traffic);  % Total number of users. Consider 1 
                               % application per user.

    % Prealocate memory for structure users. Initialize the two fields to
    % empty since the length for slots is variable for each user (number of
    % packets). i.e. users(U).slots{5} is a list containing the slots where
    % the packet (flow) 5 from user U is present.
    flows = struct('slots',{},'TH',[],'remaining',[],'deadlines',[],'failed',[],'success',[]);
    
    % Convert packets into demanded throughput following a uniform
    % distribution accros the interval <tArrival,tDeadline>
    for id = 1:Nusers
        arrivals = traffic(id).arrivalTimes;
        deadlines = traffic(id).deadlines;
        Npkt = traffic(id).numPkts;
        % Pre-allocate memory for a faster execution
        valor.slots = cell(Npkt,1);
        valor.TH = zeros(Npkt,1);
        valor.remaining = zeros(Npkt,1);
        valor.deadlines = zeros(Npkt,1);
        for pkt = 1:Npkt
            % The first slot is the right next available slot after the
            % moment of arrival (ceil). The last slot is the slots that
            % comprises the deadline time (floor)
            valor.slots{pkt} = (ceil(arrivals(pkt)/Tslot) : 1 : ...
                                floor(deadlines(pkt)/Tslot));
            nSlots = length(valor.slots{pkt});
            % The requested Throughput is computed in bits per second (bps)
            valor.TH(pkt) = traffic(id).Payload / (nSlots*Tslot*1e-3);
            valor.remaining(pkt) = traffic(id).Payload;
            valor.deadlines(pkt) = ceil(traffic(id).deadlines(pkt)/Tslot);
            valor.failed(pkt) = 0;
            valor.success(pkt) = 0;
        end
        flows(id) = valor;
    end
end