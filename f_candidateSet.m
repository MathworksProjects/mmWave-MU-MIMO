%% TODO description
function [combIds,combTH] = f_candidateSet(t,flows,selFlow)
    % Initialize variables
    Nusers = length(flows);  % Number of users
    priority = zeros(1,Nusers);  % Priority vector
    userID = (1:1:Nusers);  % User ID vector
    % priority assignation
    for id = 1:Nusers
        if selFlow(id)~=0
            deadline = flows(id).deadlines(selFlow(id));
            Delta = deadline - t;
            % Priorities are inversively proportional to the remaining time
            % (in slots) to deadline (This is Policy PLt)
            priority(id) = 1/Delta;
        end
    end
    % Sort Priorities and return candidate set
    idxPrior = priority~=0;      % Index of the users with non 0 priority
    prior = priority(idxPrior);  % Priorities
    prior(isinf(prior)) = 1.1;   % Inf is substituted by value slightly bigger than 1.
    priorID = userID(idxPrior);  % User IDs with non-0 priorities. 
                                 % This is the candidate set to be
                                 % scheduled in the slot
    % Create combinatorial matrix
    N = 2^length(prior);
    binMatrix = de2bi((0:1:N-1),length(prior));
    % Compute possible combinations
    combIds = binMatrix .* repmat(priorID,N,1);
    % Compute aggregate priority amongst users
    sumVec = binMatrix*prior.';
    % Replace NaN (priority=Inf when the Deadline is right now)
    sumVec(isnan(sumVec)) = 0;
    % Sort the aggregated priority in a descending way
    [~,idxSort] = sort(sumVec,'descend');
    % Apply Greedy policy to create a sorted candidate set matrix
    combIds = combIds(idxSort,:);
    % Create Matrix with the demanded Throughput per user per combination
    combTH = zeros(size(combIds));
    for idx = 1:size(combIds,1)
        idList = combIds(idx,:);
%         idList = idList(idList~=0);
        for id = idList
            if id~=0
                combTH(idx,id==idList) = flows(id).TH(selFlow(id));
            end
        end
    end
end