function [sol_found,W,averageCap,totTime,usersAssigned] = o_luck_and_choice_solver(problem,conf)
    % We will paralelize the solution computations: we need (if not already
    % created) a parallelization processes pool
    gcp;
    
    %  Antenna allocation solver, Best or Remove method
    tic;
    % change of orientation required for PhasedArrayTBx
    problem.NzPatch = problem.NxPatch; 
    % change of orientation required for PhasedArray TBx
    problem.dz = problem.dx;
    
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

    %% Algorithm 1: solutions for every <User, Nmax> pair
    % Complex antenna weights matrix
    W = zeros(problem.nUsers, length(problem.NmaxArray), ...
        problem.NxPatch*problem.NyPatch);
        % Received power matrix
    PRx = zeros(problem.nUsers, length(problem.NmaxArray));
        % Interferences matrix
    I = zeros(problem.nUsers, length(problem.NmaxArray), problem.nUsers);
        % structs cell containing every parameter in the solution files
    sol = cell(problem.nUsers, length(problem.NmaxArray));
    % For the parallelization to be possible, we need to copy the conf struct
    % for each user, specifying the value of IDUserAssigned
    cfg = cell(1,problem.nUsers);
    prblm = cell(1,problem.nUsers);
    for u = 1:problem.nUsers
        cfg{u} = conf;
        prblm{u} = problem;
        prblm{u}.IDUserAssigned = u;
    end
    parfor u = 1:problem.nUsers
        sol_temp = cell(1,length(prblm{u}.NmaxArray));
        W_temp = zeros(length(prblm{u}.NmaxArray),...
            prblm{u}.NxPatch*prblm{u}.NyPatch);
        PRx_temp = zeros(1,length(prblm{u}.NmaxArray));
        I_temp = zeros(length(prblm{u}.NmaxArray),prblm{u}.nUsers);
        for i = 1:length(prblm{u}.NmaxArray)
            if cfg{u}.randomSolution
                [sol_temp{i},W_temp(i,:),PRx_temp(i),I_temp(i,:)] = ...
                    o_solveSingleNmaxUserInstance(cfg{u},prblm{u},...
                    prblm{u}.NmaxArray(i),'random');
            else
                [sol_temp{i},W_temp(i,:),PRx_temp(i),I_temp(i,:)] = ...
                    o_solveSingleNmaxUserInstance(cfg{u},prblm{u},...
                    prblm{u}.NmaxArray(i));
            end
        end
        sol{u,:} = sol_temp;
        W(u,:,:) = W_temp;
        PRx(u,:) = PRx_temp;
        I(u,:,:) = I_temp;
    end

    % px, py and px are independent from Nmax or user ID:
    pz = problem.possible_locations(1,:);
    py = problem.possible_locations(2,:);
    px = problem.possible_locations(3,:);

    % Create Patch: 3 rows (x, y, z coordinates), problem.NxPatch*problem.NyPatch
    % columns (antennas in the patch)
    patch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);

    % Group array assignments: 3 rows (x, y, z coordinates), problem.Nmax 
    % columns (antennas selected from the patch), and N layers (number of users
    % i.e. number of different assignations from the same patch and Nmax)
    arraysCell = cell(length(problem.NmaxArray),1);
    for i = 1:length(problem.NmaxArray)
        arraysCell{i} = o_getArrays(problem.nUsers,...
            squeeze(W(:,i,:)),px,py,pz);
    end

    %% Get all combinations of displacements for every pair <user,Nmax>
    % Total number of combinations per <UE,problem.Nmax>
    nComb = zeros(problem.nUsers,length(problem.NmaxArray)); 
    % Combinations per <UE,problem.Nmax>
    comb = cell(problem.nUsers,length(problem.NmaxArray));
    % Binary combinations per <UE,problem.Nmax>
    combMat = cell(problem.nUsers,length(problem.NmaxArray));
    for n = 1:length(problem.NmaxArray)
        for u = 1:problem.nUsers
            % Compute the actual possible combinations
            array = arraysCell{n};
            % The following function fills in the respective positions of the
            % nComb, comb and combMat arrays. conf is needed for some debug and
            % plot config. parameters
            [nComb(u,n),comb{u,n},combMat{u,n},conf] = o_getComb(patch,problem.dx,...
                problem.dy,array(:,:,u),conf);
        end
    end
    %% Alg. 2: We search for the best combination complying with the requirements
    [aveCap, assignation] = o_find_best_combination(problem,conf,PRx,I,patch,comb,combMat,arraysCell);
    sol_found = false;
    if aveCap ~= -Inf
        sol_found = true;
        if conf.verbosity >= 1
            fprintf('The average capacity achieved is %f\n', aveCap);
        end
    else
        fprintf('The solution does not satisfy the capacity constraints.\n');
    end
    W = [];
    averageCap = aveCap;
    totTime = toc;
    usersAssigned = sum(find(assignation~=0));
end

