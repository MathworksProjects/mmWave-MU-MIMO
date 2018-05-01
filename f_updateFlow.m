%% flows = f_updateFlow(flows,selFlow,finalSet,Tslot)
% This function updates the remaining bits to be sent for each flow and
% returns the updated flow information. The transmitted bits are computed
% as the expected throughput.
%
% In addition, the function redistributes the requires throughput for those
% packets that didn't make it to the receiver. For those flows that have
% failed and it was the last chance to make it (las slot available), give
% up and increase the counter of failed flows
function flows = f_updateFlow(t,flows,selFlow,finalSet,finalTH,candSet,Tslot,DEBUG)
    % Packets with PER=0. Update flows of packets that made it
    selFinalFlow = selFlow(finalSet);  % Iterate over the meaningful set
    for id = finalSet
        succFlow = selFinalFlow(id==finalSet);
        TXbits = finalTH(id==finalSet)*(Tslot*1e-3);
        flows(id).remaining(succFlow) = flows(id).remaining(succFlow) - TXbits;
        if DEBUG; fprintf('\tTS=%d\tID=%d\tFLOW=%d\tTX=%.1f(bits)\tremain=%.0f(bits)\n',t,id,succFlow,TXbits,flows(id).remaining(succFlow)); end
        % Check if the flow has no time left to be transmitted
        deadlineSlot = max(flows(id).slots{succFlow});
        nSlotsRem = deadlineSlot - t;
        % Determine if we have met the deadline for the indicated packet
        if (nSlotsRem==0) && (flows(id).remaining(succFlow)>0)
            % We succesfully transmitted the packet but there are still
            % bits to be tx and this was the last chance we've got.
            % Increment failed flow count
            flows(id).failed(succFlow) =  1;
            if DEBUG; fprintf('NOK\tTS=%d\tID=%d\tFLOW=%d\n',t,id,succFlow); end
        elseif flows(id).remaining(succFlow)<=0
            % We increment the number of flows that we served succesfully.
            % (It can happen that we serve them way before the deadline)
            flows(id).success(succFlow) = 1;
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
        % DeadlineSlot is the last slot available for transmission
        deadlineSlot = max(flows(id).slots{failFlow});
        nSlotsRem = deadlineSlot - t;
        % Check if the flow has no time left to be transmitted
        if nSlotsRem==0
            % The flow was not served succesfully
            flows(id).failed(failFlow) = 1;
            if DEBUG; fprintf('NOK\tTS=%d\tID=%d\tFLOW=%d\n',t,id,failFlow); end
        else
            % Redistribute flow across slots
            if DEBUG; flows(id).TH(failFlow) = rem/(nSlotsRem*Tslot*1e-3); end
        end
    end
end