function [sol_found,W,handle_ConformalArray,estObj] = f_heuristics(problem,conf,usersToBeAssigned)
% F_HEURISTICS - Computes the antenna weights and an estimation of the 
%                capacity for the MU-MIMO antenna allocation problem in mmWave
%
% Syntax:  [sol_found,W,handle_ConformalArray,estObj] = f_heuristics(problem,conf,usersToBeAssigned)
%
% Inputs:
%    conf - struct containint configuration in data/config_test.dat
%    problem - struct containint configuration in data/metaproblem_test.dat
%    usersToBeAssigned: [1 x N] vector containing the N users that are to
%                       be allocated antennas
%
% Outputs:
%    sol_found - boolean variable stating whether a solution meeting the re-
%                quirements of the users has been found
%    W - [N x A] matrix containing the weights applied per user to
%        each antenna. That is, the element Wij corresponds to the
%        weight that the user i will apply to antenna j. Conversely,
%        row i represents the weights applied for user i, and the
%        column j will contain only ONE element different from 0.
%        An example for N=3 and A=4:
%
%        W = [0.0+0.0j, 0.0+0.0j, 0.1+0.4j, 0.0+0.0j;
%             0.0+0.0j, 0.4+0.2j, 0.0+0.0j, 0.8+0.3j;
%             0.4+0.7j, 0.0+0.0j, 0.0+0.0j, 0.0+0.0j]
%
%    handle_ConformalArray - phased.ConformalArray object handle containing
%                            the array used in the system
%    Cap - [1 x N] vector containing the capacity (bits/s/Hz)
%          estimation per user for the assignment obtained as a result
%          of the optimization
%          If conf.MinObjFIsSNR = true, this is no longer a capacity but
%          a SINR in dBs
%
% Other m-files required: o_solveSingleNmaxUserInstance,
%                         o_antennas_to_subarrays, f_BF_results,
%                         o_getPatch, o_getArrays, o_plot_feasible_comb
% Subfunctions: None
% MAT-files required: None
%
% See also: main
%
%------------- BEGIN CODE --------------

% Check if we inherit Weights from conventional beamformer
if isfield(problem,'initialW') && ~strcmp(problem.arrayRestriction, 'None')
    error('When inheriting Beamforming weights from LCMV or conventional beamforming, ' + ...
          'we need to set the subarray geometry to None. Please, try again.\n');
end

% Correct number of users to just focus on the ones scheduled in this slot
nUsers = length(usersToBeAssigned);

% We will paralelize the solution computations: we need (if not already
% created) a parallelization processes pool
gcp;

%% Create subarray partition
problem = o_create_subarray_partition(problem);

problem.NzPatch = problem.NxPatch;
problem.dz = problem.dx;

%% Create the antenna handler and the data structure with all possible pos.
problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
    [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
    'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque s�
handle_ConformalArray = phased.URA([problem.NyPatch,problem.NzPatch],...
    'Lattice','Rectangular','Element',problem.handle_Ant,...
    'ElementSpacing',[problem.dy,problem.dz]);

problem.possible_locations = handle_ConformalArray.getElementPosition;

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
    W = zeros(nUsers, problem.NxPatch*problem.NyPatch);
    % Received power matrix
    PRx = ones(1,nUsers)*-Inf;
    % Interferences matrix
    I = ones(nUsers, nUsers)*-Inf;
    % structs cell containing every parameter in the solution files
    sol = cell(1,nUsers);
    problem.Partition = initial_partition;
    problem.N_Subarrays = initial_N_Subarrays;
    % We will accumulate in the assignments_status var the
    % antennas / subarrays assigned as soon as we assign them
    assignments_status = zeros(1,problem.N_Subarrays);
    if conf.verbosity >= 1
        display(problem.NmaxArray);
    end
    [~,orderedIndices] = sort(problem.MinObjF,'descend');
    for u = orderedIndices
        problem.IDUserAssigned = usersToBeAssigned(u);
        if conf.verbosity >= 1
            fprintf('------------------------------------\n');
            fprintf('Solving sub-problem for the user %d:\n',...
                problem.IDUserAssigned);
        end
        % After the first user has already chosen some antennas /
        % subarrays, we remove them from the Partition in order to not to
        % select them again
        already_assigned_elem = nonzeros(assignments_status)';
        problem.Partition = initial_partition;
        problem.Partition(:,already_assigned_elem) = [];
        problem.N_Subarrays = initial_N_Subarrays - numel(already_assigned_elem);
        % Call heuristics for individual user
        if conf.randomSolution
            [sol_temp,W_temp,PRx_temp,I_temp,~] = ...
                o_solveSingleNmaxUserInstance(conf,problem,...
                problem.NmaxArray(problem.IDUserAssigned),...
                'random');
        else
            [sol_temp,W_temp,PRx_temp,I_temp,~] = ...
                o_solveSingleNmaxUserInstance(conf,problem,...
                problem.NmaxArray(problem.IDUserAssigned));
        end
        % Update the already assigned antennas
        n_selected = length(sol_temp)/3;
        subarrays_selected = o_antennas_to_subarrays(sol_temp(1:n_selected),problem.Partition);
        assignments_status = padarray([already_assigned_elem,...
            sol_temp(1:n_selected)], [0 initial_N_Subarrays - ...
            length(subarrays_selected) - ...
            nnz(assignments_status)],'post');
        sol{u} = sol_temp;
        W(u,:) = W_temp;
        PRx(u) = PRx_temp;
        I(u,:) = I_temp;
        if conf.verbosity >= 1
            fprintf('Solved!\n');
        end
    end

    % Normalize weigths
    for id = 1:nUsers
        idx = find(W(id,:)~=0.0);
        W(id,idx) = (1/sqrt(W(id,idx)*W(id,idx)')) * W(id,idx);
    end

    % Compute Capacities and SINR
    oldVerbosity = conf.verbosity;  conf.verbosity = 0;
    [~,~,Cap,SINR] = f_BF_results(W,handle_ConformalArray,usersToBeAssigned,problem,conf,false);
    conf.verbosity = oldVerbosity;

    % Parse Results properly
    if conf.MinObjFIsSNR;    estObj = db2pow(SINR).';  % Linear Scale
    else;                    estObj = db2pow(Cap).';  % Linear Scale
    end
    if any((estObj-problem.MinObjF)<0) || any((problem.MaxObjF-estObj)<0)
        aveEstObj = -Inf;
        if conf.verbosity >= 1
            fprintf('The solution does not satisfy the capacity constraints.\n');
        end
    else
        aveEstObj = mean(estObj);
        sol_found = true;
        if conf.verbosity >= 1
            fprintf('The average capacity achieved is %f\n', pow2db(mean(db2pow(Cap))));
            fprintf('The achieved capacity per user is:\n');
            display(Cap.');
        end
    end

    if conf.RefineSolution && (aveEstObj == -Inf) % No sol found...
        % Let's check how good is the antenna distribution
        capPerAnt = estObj./problem.NmaxArray;
        DeltaCap = estObj-problem.MinObjF;
        % The DeltaCap of the users that are not to be assigned should be 0
        DeltaCap(setdiff(1:nUsers,usersToBeAssigned)) = 0;
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

if conf.verbosity > 1
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
    arrays = o_getArrays(nUsers,W,px,py,pz);
    o_plot_feasible_comb(problem,conf,patch,arrays);
end

if problem.DEBUG
    nAntennasTot = 0;
    for id =1:nUsers
        nAntennas = length(W(id,W(id,:)~=0));
        fprintf("\t\t\tID=%d\t#Ant=%d\n",id,nAntennas);
        nAntennasTot = nAntennasTot + nAntennas;
    end
end


% EOF