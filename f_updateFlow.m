%% flows = f_updateFlow(flows,selFlow,finalSet,Tslot)
% This function updates the remaining bits to be sent for each flow and
% returns the updated flow information. The transmitted bits are computed
% as the expected throughput.
%
% In addition, the function redistributes the requires throughput for those
% packets that didn't make it to the receiver. For those flows that have
% failed and it was the last chance to make it (las slot available), give
% up and increase the counter of failed flows
function flows = f_updateFlow(flows,selFlow,finalSet,finalTH,candSet,Tslot)
    % Update the flows of those packets that make it through
    selFinalFlow = selFlow(finalSet);  % Iterate over the meaningful set
    for id = finalSet
        TXbits = finalTH(id==finalSet)*(Tslot*1e-3);
        flows(id).remaining = flows(id).remaining(selFinalFlow) - TXbits;
        fprintf('User %d TX %.1f bits - %d bits remain\n',id,TXbits,flows(id).remaining(selFinalFlow(id==finalSet)));
    end
    % Redistribute the flows for those packets that didn't make it through
    failSet = setdiff(candSet,finalSet);
    selFailFlow = selFlow(finalSet);  % Iterate over the meaningful set
    for id = failSet
        failFlow = selFailFlow(id==failSet);
        rem = flows(id).remaining(failFlow);
        deadline = flows(id).deadlines(failFlow);
        nSlotsRem = deadline - Tslot;
        % Check if the flow has no time left to be transmitted
        if nSlotsRem==0
            % Increment failed flow count
            flows(id).failed(selFailFlow(id==failSet)) = ...
                flows(id).failed(selFailFlow(id==failSet)) + 1;
        else
            % Redistribute flow across slots
            flows(id).TH(selFailFlow(id==failSet)) = rem/(nSlotsRem*Tslot*1e-3);
        end
    end
end