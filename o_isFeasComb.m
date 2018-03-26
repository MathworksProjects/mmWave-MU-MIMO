function [ bool ] = o_isFeasComb(patch,comb,combMat,NmaxArray,c)
    if isempty(c)
        bool = true;
        return
    end
    N = length(c);
    cb = zeros(N,1);
    maxNarray = zeros(N,1);
    for ind = 1:1:length(c)
        temp = comb{c(ind),NmaxArray(c(ind))};
        maxNarray(ind) = size(temp,3); % Number of possible combinations for this user (cb is shifted to start from 0)
    end
    nCombGrTZero = maxNarray(maxNarray>0);
    nComb = prod(nCombGrTZero);
    bool = false;
    
    for i = 1:nComb
        totsum = patch;
        % Contains the index of every combination for every node. If N=2
        % and there are 4 combinations, cb(1) is a vector defining the
        % combination indices for UE1 and UE2.
        for ind = 1:1:N
            % totSum contains the number of times the antennas have been
            % used amongst the N UE's
            if maxNarray(ind) > 0
                temp = combMat{c(ind),NmaxArray(c(ind))};  % Get combinations for user n
                mat = temp(:,:,cb(ind)+1);  % Get combination cb(n) for user n
                totsum(3,:) = totsum(3,:) + mat(3,:);  % Sum antenna usage
            end
        end
        % Find feasible set - If all of them have only been used once, then
        % the combination is feasible
        if ~any(totsum(3,:) > 1)
            bool = true;
            return
        end
        cb = o_sumOneToCombination(cb, maxNarray-1);
    end
end