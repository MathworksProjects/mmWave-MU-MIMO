function flows = f_updateFlow(t,flows,selFlow,finalSet,THiter,candSet,Tslot,DEBUG)
% F_UPDATEFLOW - Updates the remaining bits to be sent for each flow and
% returns the updated flow information. The transmitted bits are computed
% as the expected throughput.
%
% In addition, the function redistributes the requires throughput for those
% packets that didn't make it to the receiver. For those flows that have
% failed and it was the last chance to make it (las slot available), give
% up and increase the counter of failed flows
%
% Syntax:  flows = f_updateFlow(t,flows,selFlow,finalSet,THiter,candSet,Tslot,DEBUG)
%
% Inputs:
%    t - Scalar value that represents the slot ID in the DES
%    flows - Array of structs of length equal to the number of
%            users. For each user, each flow (belonging to each packet) is
%            characterized by the amount of bits that needs to be delivered. The
%            amount of bits are distributed uniformly across the slots until
%            reaching the deadline. Thus, the variable flows contains four features:
%            - slots:     For each flow, the slots across it.
%            - TH:        The average throughput demanded for that flow over the slots
%                         indicated in the slots field.
%            - remaining: We set the remaining field of the flow to be the value of
%                         the Payload.
%            - deadlines: The deadline of the actual packet in slot ID. 
%                         This is used in the future to set priorities.
%            - failed:    Mark whether the flow has failed to be served before the
%                         deadline (0 or 1).
%            - success:   Mark whether the flow has been succesfully served before the
%                         deadline (0 or 1).
%            - maxSlot:   Maximum deadline slot used as simulation time or Tsym in
%                         the system in case trafficType is 'dataSet'.
%    selFlow - For each user, we select a particular flow that maps with
%              the current time slot in the simulator. Only one flow per user can be
%              selected, per slot. Aggregation of upper layer flows has already
%              happened at this point.
%    finalSet - 
%    THiter - Required Throughput per user
%    candSet - Users ID being considered in the current slot
%    Tslot - Time slot duration in ms
%    DEBUG - Debug flag (True or False)
%
% Outputs:
%    flows -Updated flow structure
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: main
%
%------------- BEGIN CODE --------------

% Packets with PER=0. Update flows of packets that made it
selFinalFlow = selFlow(finalSet);  % Iterate over the meaningful set
for id = finalSet
    succFlow = selFinalFlow(id==finalSet);
    TXbits = THiter(id)*(Tslot*1e-3);
    flows(id).remaining(succFlow) = flows(id).remaining(succFlow) - TXbits;
    if DEBUG; fprintf('\tTS=%d\tID=%d\tFLOW=%d\tTX=%.1f(bits)\tremain=%.0f(bits)\n',t,id,succFlow,TXbits,flows(id).remaining(succFlow)); end
    % Check if the flow has no time left to be transmitted
    deadlineSlot = flows(id).deadlines(succFlow);
    nSlotsRem = deadlineSlot - t;
    % Determine if we have met the deadline for the indicated packet
    if (nSlotsRem==0) && (flows(id).remaining(succFlow)>1)
        % We succesfully transmitted the packet but there are still
        % bits to be tx and this was the last chance we've got.
        % Increment failed flow count
        if DEBUG; fprintf('NOK\tTS=%d\tID=%d\tFLOW=%d\n',t,id,succFlow); end
    elseif (flows(id).remaining(succFlow)<1)
        % We increment the number of flows that we served succesfully.
        % (It can happen that we serve them way before the deadline)
        flows(id).success(succFlow) = 1;
        flows(id).failed(succFlow)  = 0;
        if DEBUG; fprintf('OK\tTS=%d\tID=%d\tFLOW=%d\n',t,id,succFlow); end
    end
    % Update slots in the flow. If we have served all the bits before
    % the deadline, need to remove the slots pertaining to that flow.
    if nSlotsRem>0 && flows(id).remaining(succFlow)<=0
        updSlots = (t+1:deadlineSlot);
        [~,idxUpdate] = intersect(flows(id).slots{succFlow},updSlots);
        flows(id).slots{succFlow}(idxUpdate) = [];
    end
end
% Packets with PER=1. Redistribute flows for packets that didn't made it
failSet = setdiff(candSet,finalSet);
selFailFlow = selFlow(failSet);  % Iterate over the meaningful set
for id = failSet
    failFlow = selFailFlow(id==failSet);
    rem = flows(id).remaining(failFlow);
    % Check if the flow has no time left to be transmitted
    deadlineSlot = flows(id).deadlines(failFlow);
    nSlotsRem = deadlineSlot - t;
    if nSlotsRem==0
        % The flow was not served succesfully
        if DEBUG; fprintf('NOK\tTS=%d\tID=%d\tFLOW=%d\n',t,id,failFlow); end
    else
        % Redistribute flow across slots
        if DEBUG; flows(id).TH(failFlow) = rem/(nSlotsRem*Tslot*1e-3); end
    end
end



% EOF