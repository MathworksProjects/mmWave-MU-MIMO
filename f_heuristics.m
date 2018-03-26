function [sol_found,W,Cap] = f_heuristics(problem,conf,usersToBeAssigned)
    % We will paralelize the solution computations: we need (if not already
    % created) a parallelization processes pool
    p = gcp('nocreate');
    
    %% Create subarray partition
    problem = o_create_subarray_partition(problem);

    %% Create the antenna handler and the data structure with all possible pos.
    problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
        [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
        'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque sí
    handle_URA = phased.URA([problem.NyPatch,problem.NzPatch],...
        'Lattice','Rectangular','Element',problem.handle_Ant,...
    'ElementSpacing',[problem.dy,problem.dz]);

    problem.possible_locations = handle_URA.getElementPosition;

    initial_partition = problem.Partition;
    initial_N_Subarrays = problem.N_Subarrays;

    % Boolean flag indicating if we have already found a feasible solution
    sol_found = false;
    if conf.verbosity >= 1
        fprintf('=========================================\n');
        fprintf('Solving problem with the following users:\n');
        disp(usersToBeAssigned);
    end
    problem.Partition = initial_partition;
    problem.N_Subarrays = initial_N_Subarrays;
    problem = o_compute_antennas_per_user(problem,usersToBeAssigned);
    RoomForImprovement = true;
    while RoomForImprovement
        % Restore the initial values
            % Complex antenna weights matrix
        W = zeros(problem.nUsers, problem.NxPatch*problem.NyPatch);
            % Received power matrix
        PRx = ones(1,problem.nUsers)*-Inf;
            % Interferences matrix
        I = ones(problem.nUsers, problem.nUsers)*-Inf;
            % structs cell containing every parameter in the solution files
        sol = cell(1,problem.nUsers);
        problem.Partition = initial_partition;
        problem.N_Subarrays = initial_N_Subarrays;
        % We will accumulate in the assignments_status var the
        % antennas / subarrays assigned as soon as we assign them
        assignments_status = zeros(1,problem.N_Subarrays);
        if conf.verbosity >= 1
            display(problem.NmaxArray);
        end
        for u = 1:problem.nUsers
            problem.IDUserAssigned = u;      
            if conf.verbosity >= 1
                fprintf('------------------------------------\n');
                fprintf('Solving sub-problem for the user %d:\n',...
                    problem.IDUserAssigned);            
            end
            % After the first user has already chosen some antennas /
            % subarrays, we remove them from the Partition in order to not to
            % select them again
            already_assigned_elem = nonzeros(assignments_status)';
            if u > 1
                problem.Partition = initial_partition;
                problem.Partition(:,already_assigned_elem) = [];
                problem.N_Subarrays = initial_N_Subarrays - ...
                    numel(already_assigned_elem);
            end
            W_temp = zeros(1,problem.NxPatch*problem.NyPatch);
            PRx_temp = -Inf;
            I_temp = ones(1,problem.nUsers)*-Inf;
            % if the user is not to be assigned, no need to do anything else
            if ismember(u,usersToBeAssigned)
                % Actually solve the user's assignment (observe that Nmax index
                % is 1)
                if conf.randomSolution
                    [sol_temp,W_temp,PRx_temp,I_temp] = ...
                        o_solveSingleNmaxUserInstance(conf,problem,...
                        problem.NmaxArray(problem.IDUserAssigned),...
                        'random');
                else
                    [sol_temp,W_temp,PRx_temp,I_temp] = ...
                        o_solveSingleNmaxUserInstance(conf,problem,...
                        problem.NmaxArray(problem.IDUserAssigned));
                end
                % Update the already assigned antennas
                assignments_status = padarray([already_assigned_elem,...
                    sol_temp(1:problem.NmaxArray(u))], [0 initial_N_Subarrays - ...
                    problem.NmaxArray(u) - nnz(assignments_status)],'post');
            end
            sol{u} = sol_temp;
            W(u,:) = W_temp;
            PRx(u) = PRx_temp;
            I(u,:) = I_temp;
            if conf.verbosity >= 1
                fprintf('Solved!\n')
            end
        end

        [aveCap, Cap] = o_compute_averageCap_maxminthr(PRx,I,problem.Noise,...
            problem.MaxThr,problem.MinThr,usersToBeAssigned);
        if aveCap ~= -Inf
            sol_found = true;
            if conf.verbosity >= 1
                fprintf('The average capacity achieved is %f\n', aveCap);
                display(Cap);
            end
        else
            if conf.verbosity >= 1
                fprintf('The solution does not satisfy the capacity constraints.\n');
            end
        end

        if conf.RefineSolution && (aveCap == -Inf) % No sol found...
            % Let's check how good is the antenna distribution
            capPerAnt = Cap./problem.NmaxArray;
            DeltaCap = Cap-problem.MinThr;
            % The DeltaCap of the users that are not to be assigned should be 0
            DeltaCap(setdiff(1:problem.nUsers,usersToBeAssigned)) = 0;
            AvailableAnt = floor(DeltaCap./capPerAnt);
            usersNotInNeed = find(AvailableAnt > 0);
            SpareAntennas = AvailableAnt(usersNotInNeed);
            totAvailAnt = sum(SpareAntennas);
            if conf.verbosity >= 1
                display(AvailableAnt);
            end
            if totAvailAnt == 0
                % If there are no antennas available, there's nothing to do
                RoomForImprovement = false;
                if conf.verbosity >= 1
                    fprintf("Nothing to do to refine the solution\n");
                end
            else
                usersInNeed = find(AvailableAnt < 0);
                NeededAntennas = AvailableAnt(usersInNeed);
                [~,indexes] = sort(NeededAntennas);
                u = 1;
                if conf.verbosity >= 1
                    display(problem.NmaxArray);
                end
                totalAntMoved = 0;
                while totAvailAnt > 0 && u <= length(usersInNeed)
                   antennasToMove = min(totAvailAnt,...
                       abs(NeededAntennas(indexes(u))));
                   problem.NmaxArray(usersInNeed(indexes(u))) = ...
                       problem.NmaxArray(usersInNeed(indexes(u))) + ...
                       antennasToMove;
                   totAvailAnt = totAvailAnt - antennasToMove;
                   totalAntMoved = totalAntMoved + antennasToMove;
                   u = u + 1;
                end
                [~,indexes] = sort(SpareAntennas,'descend');
                u = 1;
                while totalAntMoved > 0 && u <= length(usersNotInNeed)
                   antennasToMove = min(totalAntMoved,...
                       SpareAntennas(indexes(u)));
                   problem.NmaxArray(usersNotInNeed(indexes(u))) = ...
                       problem.NmaxArray(usersNotInNeed(indexes(u))) - ...
                       antennasToMove;
                   totalAntMoved = totalAntMoved - antennasToMove;
                   u = u + 1;
                end
                if conf.verbosity >= 1
                    display(problem.NmaxArray);
                end
            end
        else
            RoomForImprovement = false;
        end
    end
end

