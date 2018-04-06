function [problem] = o_compute_antennas_per_user(problem,usersToBeAssigned)
% Let's first compute the Nmax for every user
        problem.NmaxArray = zeros(1,problem.nUsers);
        minSub = 0;
        [~,orderedIndices] = sort(problem.MinThr,'descend');
        for u = orderedIndices
            if problem.MinThr(u) > 0 && minSub < problem.N_Subarrays
                problem.NmaxArray(u) = 1; % At least 1
                minSub = minSub + 1;
            end
        end
        sum_MinThr = sum(problem.MinThr(usersToBeAssigned));
        if sum_MinThr > 0
            thr_ant_converter = (problem.N_Subarrays-minSub)/sum_MinThr;
        else
            thr_ant_converter = (problem.N_Subarrays-minSub)/numel(usersToBeAssigned);
        end
        for u = usersToBeAssigned
            % Now NmaxArray is the Nmax value chosen for each user
            problem.NmaxArray(u) = max(0,problem.NmaxArray(u) + ...
                floor(problem.MinThr(u)*thr_ant_converter));
        end
        % The previous division might have left some antennas / subarrays
        % unassigned, due to rounding errors. Let's assign the remaining 
        % antennas / subarrays to the user most in need
        sumAnt = 0;
        for u = usersToBeAssigned; sumAnt = sumAnt + problem.NmaxArray(u);end
        if sumAnt < problem.N_Subarrays % Se podría equilibrar el reparto...
            [~,i] = max(problem.MinThr(usersToBeAssigned));
            userToAssignExtraAnt = usersToBeAssigned(i);
            problem.NmaxArray(userToAssignExtraAnt) = ...
                problem.NmaxArray(userToAssignExtraAnt) + ...
                (problem.N_Subarrays-sumAnt);
        end
end

