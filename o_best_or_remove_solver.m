function [sol_found,W,averageCap,totTime] = o_best_or_remove_solver(problem,conf)
    %  Antenna allocation solver, Best or Remove method
    tic;

    %% Iteration over partial solutions to maximize utilisation and users assigned

    % Array containing the users to be assigned (it will be reduced at every
    % iteration, if no solution found. In principle, we should try assigning
    % all
    usersToBeAssigned = 1:problem.nUsers;
    
    % change of orientation required for PhasedArrayTBx
    problem.NzPatch = problem.NxPatch; 
    % change of orientation required for PhasedArray TBx
    problem.dz = problem.dx;

    % Boolean flag indicating if we have already found a feasible solution
    sol_found = false;
    while ~sol_found || isempty(usersToBeAssigned)
        [sol_found,W,Cap] = f_heuristics(problem,conf,usersToBeAssigned);
        
        if conf.verbosity >= 1
            % px, py and px are independent from Nmax or user ID:
            pz = problem.possible_locations(1,:);
            py = problem.possible_locations(2,:);
            px = problem.possible_locations(3,:);

            % Create Patch: 3 rows (x, y, z coordinates), 
            % problem.NxPatch*problem.NyPatch columns (antennas in the patch)
            patch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);

            % Group array assignments: 3 rows (x, y, z coordinates), problem.Nmax 
            % columns (antennas selected from the patch), and N layers (number of users
            % i.e. number of different assignations from the same patch and Nmax)
            arrays = o_getArrays(problem.nUsers,max(problem.NmaxArray),W,px,py,pz);
            o_plot_feasible_comb(problem,conf,patch,arrays);
        end

        % For the next iteration (we will only run it if no sol_found) we
        % remove the most consuming user
        [~, maxIndex] = max(problem.MinThr(usersToBeAssigned));
        usersToBeAssigned(maxIndex) = [];
        if isempty(usersToBeAssigned)
            fprintf('No feasible solution found to serve any user!\n');
        end
    end
    averageCap = mean(Cap);
    totTime = toc;
end

