%% Antenna allocation solver for massive arrays in mm-Wave and THz-band
    % This script is the main program of the MATLAB solver for the antenna
    % allocation problem in massive arrays for mm-Wave and THz-band.
    % Version 2, not backward-compatible.
    % Other solver configuration parameters, such as the policy used or the 
    % DEBUG and plot flags, are also part of the conf structure, and are 
    % explained further below.

%   Nothing is left behind...
clear; close all; clc;
    
% Let's read the configuration parameters for the solver. Edit this script
% to change the parameters
read_config; % Every config parameter is contained in the conf struct

% Let's process the input parameters and read the input file (problem). Edit 
% this script to change the parameters
read_input_problem; % Everything is stored in the problem struct

%% Now it's time to execute the actual solver
%  Antenna allocation solver, Best or Remove method

tic;
problem.Bw = 0.2e9; % to be further investigated...
% change of orientation required for PhasedArrayTBx
problem.NzPatch = problem.NxPatch; 
% change of orientation required for PhasedArray TBx
problem.dz = problem.dx;
problem.N_Antennas = problem.NyPatch*problem.NzPatch;
problem.N_Subarrays = problem.NxSubarrays*problem.NySubarrays;
% problem.Partition contains the antenna arrangement into the different
% subarrays, i.e. for each subarrays, you will find the set of antenna
% indexes associated to it.
problem.Partition = cell(1,problem.N_Subarrays);

% Now we fill in the problem.Partition cell, depending on the restriction
% imposed on it:
% - 'Localized': the subarrays are rectangular contiguous divisions of the
%    array
% - 'DiagInterleaved': the subarrays are interleaved in a diagonal fashion
% - 'Interleaved': the subarrays are interleaved forming rectangular arrays
%    with larger distance between elements.
if strcmp(problem.arrayRestriction,'Localized')
    for nx=0:(problem.NxSubarrays-1)
        for ny=0:(problem.NySubarrays-1)
            x0 = floor(nx/problem.NxSubarrays*problem.NxPatch);
            xend = floor((nx+1)/problem.NxSubarrays*problem.NxPatch);
            for x=x0:(xend-1)
                y0 = floor(ny/problem.NySubarrays*problem.NyPatch);
                yend = floor((ny+1)/problem.NySubarrays*problem.NyPatch);
                for y=y0:(yend-1)
                    problem.Partition{nx*problem.NySubarrays+ny+1} = ...
                        [problem.Partition{nx*problem.NySubarrays+ny+1},...
                        (x*problem.NyPatch+y+1)];
                end
            end
        end
    end
elseif strcmp(problem.arrayRestriction,'DiagInterleaved')
    if problem.NxSubarrays ~= problem.NySubarrays
        if  (problem.NxSubarrays <= problem.NySubarrays && ...
                problem.NxPatch ~= problem.NxSubarrays) || ...
                (problem.NySubarrays <= problem.NxSubarrays && ...
                problem.NyPatch ~= problem.NySubarrays)
            fprintf(['ERROR: problem.NxSubarrays (%d) and',...
                ' problem.NySubarrays',...
                ' (%d) do not match with a correct diagonally interleaved',...
                ' configuration\n'], ...
            problem.NxSubarrays, problem.NySubarrays);
            return
        end
    end
    problem.N_Subarrays = max(problem.NxSubarrays,problem.NySubarrays);
    for x=0:(problem.NxPatch-1)
        for y=0:(problem.NyPatch-1)
            problem.Partition{mod(x+y,problem.NySubarrays)+1} = ...
                [problem.Partition{mod(x+y,problem.NySubarrays)+1},...
                x*problem.NyPatch+y+1];
        end
    end
elseif strcmp(problem.arrayRestriction,'Interleaved')
    for x=0:(problem.NxPatch-1)
        for y=0:(problem.NyPatch-1)
            problem.Partition{mod(y,problem.NySubarrays)+1+...
                mod(x,problem.NxSubarrays)*problem.NySubarrays} = ...
                [problem.Partition{mod(y,problem.NySubarrays)+1+...
                mod(x,problem.NxSubarrays)*problem.NySubarrays},...
                x*problem.NyPatch+y+1];
        end
    end    
elseif strcmp(problem.arrayRestriction,'None')
    problem.Partition = cell(1,problem.N_Antennas);
    for a=1:problem.N_Antennas
        problem.Partition{a} = a;
    end
    problem.N_Subarrays = problem.N_Antennas;
end
%% Randomization and conversion of the users' positions in the space
%%%% Only in the case where no determinist position is given in the input
%%%% file defining the problem. The conversion to degrees is performed if
%%%% angles were given in radians
pd = makedist('Normal');
pd.sigma = 45;
t = truncate(pd,-45,45);
if isempty(problem.thetaUsers)
    fprintf('===================================================\n');
    fprintf('The users were not assigned positions in the space:\n');
    fprintf('Assigning random values...\n');
    pause(2);
    problem.thetaUsers = random(t,1,problem.nUsers);
    fprintf('New elevations assigned:\n');
    display(problem.thetaUsers);
elseif problem.anglesInRadians
    problem.thetaUsers = problem.thetaUsers/(2*pi)*360;
end
if isempty(problem.phiUsers)
    problem.phiUsers = random(t,1,problem.nUsers);
    fprintf('New azimuths assigned:\n');
    display(problem.phiUsers);
elseif problem.anglesInRadians
    problem.phiUsers = problem.phiUsers/(2*pi)*360;
end
% Now same thing with channel paths AoA. t could be changed if needed
% pd = makedist('Normal');
% pd.sigma = 45;
% t = truncate(pd,-45,45);

% Boolean to check whether we have already selected the channels to remove
% for each user
indexToRemoveSelected = false;
if isempty(problem.thetaChannels)
    fprintf('===================================\n');
    fprintf('The channels were not assigned AoA:\n');
    fprintf('Assigning random values...\n');
    pause(2);
    problem.thetaChannels = random(t,problem.nUsers,problem.maxnChannelPaths);
    userWithMaxNChannels = randi(problem.nUsers);
    indexToRemove = ones(1,problem.nUsers)*(problem.maxnChannelPaths+1);
    for i=setdiff(1:problem.nUsers,userWithMaxNChannels)
        % Now we select some channels (from indexToRemove to end) and remove 
        % them (not every user should have maxnChannelPaths
        indexToRemove(i) = max(1,randi(problem.maxnChannelPaths)+1);
        problem.thetaChannels(i,indexToRemove(i):end) = -Inf;
        problem.phiChannels(i,indexToRemove(i):end) = -Inf;
        problem.alphaChannels(i,indexToRemove(i):end) = -Inf;
    end
    indexToRemoveSelected = true;
    fprintf('New elevations assigned:\n');
    display(problem.thetaChannels);
elseif problem.anglesInRadians
    problem.thetaChannels = problem.thetaChannels/(2*pi)*360;
end
if isempty(problem.phiChannels)
    problem.phiChannels = random(t,problem.nUsers,problem.maxnChannelPaths);
    if indexToRemoveSelected
        for i=1:problem.nUsers
            problem.phiChannels(i,indexToRemove(i):end) = -Inf;
        end
    else
        userWithMaxNChannels = randi(problem.nUsers);
        indexToRemove = ones(1,problem.nUsers)*(problem.maxnChannelPaths+1);
        for i=setdiff(1:problem.nUsers,userWithMaxNChannels)
            % Now we select some channels (from indexToRemove to end) and remove 
            % them (not every user should have maxnChannelPaths
            indexToRemove(i) = max(1,randi(problem.maxnChannelPaths)+1);
            problem.thetaChannels(i,indexToRemove(i):end) = -Inf;
            problem.phiChannels(i,indexToRemove(i):end) = -Inf;
            problem.alphaChannels(i,indexToRemove(i):end) = -Inf;
        end
        indexToRemoveSelected = true;
    end    
    fprintf('New azimuths assigned:\n');
    display(problem.phiChannels);
elseif problem.anglesInRadians
    problem.phiChannels = problem.phiChannels/(2*pi)*360;
end
% Finally, alpha gains of each path
pd = makedist('Normal');
pd.mu = 0.5;
pd.sigma = 0.25;
t = truncate(pd,0,1);
if isempty(problem.alphaChannels)
    problem.alphaChannels = random(t,problem.nUsers,problem.maxnChannelPaths);
    if indexToRemoveSelected
        for i=1:problem.nUsers
            problem.alphaChannels(i,indexToRemove(i):end) = -Inf;
        end
    else
        userWithMaxNChannels = randi(problem.nUsers);
        indexToRemove = ones(1,problem.nUsers)*(problem.maxnChannelPaths+1);
        for i=setdiff(1:problem.nUsers,userWithMaxNChannels)
            % Now we select some channels (from indexToRemove to end) and remove 
            % them (not every user should have maxnChannelPaths
            indexToRemove(i) = max(1,randi(problem.maxnChannelPaths)+1);
            problem.thetaChannels(i,indexToRemove(i):end) = -Inf;
            problem.phiChannels(i,indexToRemove(i):end) = -Inf;
            problem.alphaChannels(i,indexToRemove(i):end) = -Inf;
        end
        indexToRemoveSelected = true;
    end    
    fprintf('New gains assigned:\n');
    display(problem.alphaChannels);
end

%% Create the antenna handler and the data structure with all possible pos.
handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
    [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
    'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque sí
handle_URA = phased.URA([problem.NyPatch,problem.NzPatch],...
    'Lattice','Rectangular','Element',handle_Ant,...
'ElementSpacing',[problem.dy,problem.dz]);

problem.possible_locations = handle_URA.getElementPosition;

%% Iteration over partial solutions to maximize utilisation and users assigned

% Array containing the users to be assigned (it will be reduced at every
% iteration, if no solution found. In principle, we should try assigning
% all
usersToBeAssigned = 1:problem.nUsers;

initial_partition = problem.Partition;
initial_N_Subarrays = problem.N_Subarrays;

% Boolean flag indicating if we have already found a feasible solution
sol_found = false;
while ~sol_found
    fprintf('=========================================\n');
    fprintf('Solving problem with the following users:\n');
    disp(usersToBeAssigned);
    problem.Partition = initial_partition;
    problem.N_Subarrays = initial_N_Subarrays;
    % Let's first compute the Nmax for every user
    problem.NmaxArray = zeros(1,problem.nUsers);
    minSub = 0;
    for u = usersToBeAssigned
        if problem.MinThr(u) > 0
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
        problem.NmaxArray(u) = problem.NmaxArray(u) + ...
            floor(problem.MinThr(u)*thr_ant_converter);
    end
    % The previous division might have left some antennas / subarrays
    % unassigned, due to rounding errors. Let's assign the remaining 
    % antennas / subarrays to the user most in need
    sumAnt = 0;
    for u = usersToBeAssigned; sumAnt = sumAnt + problem.NmaxArray(u);end
    if sumAnt < problem.N_Subarrays
        [~,i] = max(problem.MinThr(usersToBeAssigned));
        problem.NmaxArray(i) = problem.NmaxArray(i) + (problem.N_Subarrays-sumAnt);
    end
    RoomForImprovement = true;
    while conf.RefineSolution && RoomForImprovement
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
        display(problem.NmaxArray);
        for u = 1:problem.nUsers
            problem.IDUserAssigned = u;
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
                        solveSingleNmaxUserInstance(conf,problem,handle_Ant,...
                        'random');
                else
                    [sol_temp,W_temp,PRx_temp,I_temp] = ...
                        solveSingleNmaxUserInstance(conf,problem,handle_Ant);
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
        end

        [aveCap, Cap] = compute_averageCap_maxminthr(PRx,I,problem.Noise,...
            problem.MaxThr,problem.MinThr,usersToBeAssigned);
        if aveCap ~= -Inf
            sol_found = true;
        end
        fprintf('The average capacity achieved is %f\n', aveCap);

        capPerAnt = Cap./problem.NmaxArray;
        display(problem.NmaxArray);
        display(capPerAnt);
        DeltaCap = Cap-problem.MinThr;
        DeltaCap(setdiff(1:problem.nUsers,usersToBeAssigned)) = 0;
        display(DeltaCap);
        AvailableAnt = floor(DeltaCap./capPerAnt);
        display(AvailableAnt);
        usersWithAvailAnt = find(DeltaCap > 0);
        totAvailAnt = sum(AvailableAnt(usersWithAvailAnt));
        if totAvailAnt == 0
            RoomForImprovement = false;
        else
            for u=usersWithAvailAnt
                problem.NmaxArray(u) = problem.NmaxArray(u)-AvailableAnt(u);
            end
            usersInNeed = find(DeltaCap < 0);
            totNeededCap = sum(abs(DeltaCap(usersInNeed)));
            for u=usersInNeed
               problem.NmaxArray(u) = problem.NmaxArray(u) + ...
                   totAvailAnt*abs(DeltaCap(u))/totNeededCap;
            end
            display(problem.NmaxArray);
        end
    end
    % px, py and px are independent from Nmax or user ID:
    pz = problem.possible_locations(1,:);
    py = problem.possible_locations(2,:);
    px = problem.possible_locations(3,:);

    % Create Patch: 3 rows (x, y, z coordinates), 
    % problem.NxPatch*problem.NyPatch columns (antennas in the patch)
    patch = getPatch(problem.NxPatch,problem.NyPatch,px,py);

    % Group array assignments: 3 rows (x, y, z coordinates), problem.Nmax 
    % columns (antennas selected from the patch), and N layers (number of users
    % i.e. number of different assignations from the same patch and Nmax)
    arrays = getArrays(problem.nUsers,max(problem.NmaxArray),W,px,py,pz);
    
    plot_feasible_comb;
    
    % For the next iteration (we will only run it if no sol_found) we
    % remove the most consuming user
    [~, maxIndex] = max(problem.MinThr(usersToBeAssigned));
    usersToBeAssigned(maxIndex) = [];
    if isempty(usersToBeAssigned)
        fprintf('No feasible solution found to serve any user!\n');
        sol_found = true; % We want to break cleanly
    end
end
toc;