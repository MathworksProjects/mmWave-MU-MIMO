function [sol_found,W,averageCap,totTime,usersAssigned] = o_best_or_remove_solver(problem,conf)
    %  Antenna allocation solver, Best or Remove method
    tic;

    %% Iteration over partial solutions to maximize utilisation and users assigned

    % Array containing the users to be assigned (it will be reduced at every
    % iteration, if no solution found. In principle, we should try assigning
    % all
    usersToBeAssigned = 1:problem.nUsers;

    % Boolean flag indicating if we have already found a feasible solution
    sol_found = false;
    while ~sol_found && ~isempty(usersToBeAssigned)
        [sol_found,W,Cap] = f_heuristics(problem,conf,usersToBeAssigned);

        % For the next iteration (we will only run it if no sol_found) we
        % remove the most consuming user
        if ~sol_found
            [~, maxIndex] = max(problem.MinThr(usersToBeAssigned));
            usersToBeAssigned(maxIndex) = [];
            if isempty(usersToBeAssigned)
                fprintf('No feasible solution found to serve any user!\n');
            end
        end
    end
    averageCap = mean(Cap);
    totTime = toc;
    usersAssigned = usersToBeAssigned;
end

