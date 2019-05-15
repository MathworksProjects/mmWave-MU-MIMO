function [W_LCMV,W_CBF,handle_ConformalArray,estObj_LCMV,estObj_CBF,candSet] = f_conventionalBF(problem,conf,candSet,refine)
%
% Syntax:  [W_LCMV,W_CBF,handle_ConformalArray,estObj_LCMV,estObj_CBF] = ...
%                                    f_conventionalBF(problem,conf,candSet)
%
% Inputs:
%    problem - struct containint configuration in data/metaproblem_test.dat
%    conf - Struct containing configuration in data/config_test.dat
%    candSet - Vector containing the users ID being considered in the 
%              current slot
%    refine - RBoolean, Refine to find solution, select users with highest
%             priority
%
% Outputs:
%    W_LCMV - Matrix [nUser x nAntennas] containg the weights for LCMV
%    W_CBF - Matrix [nUser x nAntennas] containg the weights for CBF
%    handle_ConformalArray - Initial conformal array
%    estObj_LCMV - Onjective function value for LCMV. Could be either SINR
%                  or Capacity (controlled by conf.MinObjFIsSNR)
%    estObj_CBF - Onjective function value for CBF. Could be either SINR
%                 or Capacity (controlled by conf.MinObjFIsSNR)
%    candSet - We refine the users allocated given the limitations of the
%    current beamforming algorithms. For instance, a good solution is found
%    when the number of allocated antennas is larger than the number of
%    constraints (angles to maximize and minimize)
%
% Example: 
%   problem = o_read_input_problem('data/metaproblem_test.dat');
%   conf = o_read_config('data/config_test.dat');
%   conf.verbosity = 1;  % To visualize metrics on command line
%   nUsers = 5;  % Fix number of users manually for example
%   [problem,~,~] = f_configuration(conf,problem);
%   problem.N_Antennas = nAntennas;  % Select number of antennas
%   problem.NxPatch = floor(sqrt(problem.N_Antennas));  % Adjust
%   problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);  % Adjust
%   problem.N_Antennas = problem.NxPatch.*problem.NyPatch;  % Adjust
%   candSet = [2 4 5];  % Users 1 and 3 are left out
%   [W_LCMV,W_CBF,handle_ConformalArray,estObj_LCMV,estObj_CBF] = ...
%                                    f_conventionalBF(problem,conf,candSet);
%
% Other m-files required: f_configuration
% Subfunctions: none
% MAT-files required: none
% DAT-files required: data/metaproblem_test.dat and data/config_test.dat
%
% See also: main.m,  f_BF_results,  f_heuristics.m,  main_runnable.m

%------------- BEGIN CODE --------------

% Initial number of users
nUsers = length(candSet);

% Restrict sub-arrays to Localized for LCMV
problem.arrayRestriction = 'Localized';

% Compute number of sub-arrays to assign per user. We ensure each user
% receives one array and locate them horizontaly
if ceil(sqrt(problem.N_Antennas)) > nUsers
    problem.NySubarrays = nUsers;
    problem.NxSubarrays = 1;
else
    if mod(nUsers,2)~=0;  t = factor(nUsers + 1);  % odd
    else;                 t = factor(nUsers);  % even
    end
    t = [t(1) prod(t(2:end))];
    problem.NxSubarrays = t(1);
    problem.NySubarrays = t(2);
end
problem.N_Subarrays = problem.NxSubarrays*problem.NySubarrays;
problem = o_compute_antennas_per_user(problem,candSet); 

% Create initial assignation to users
problem = o_create_subarray_partition(problem);

% Distribute partitions amongst users
totSubArrays_1 = floor((problem.NxSubarrays * problem.NySubarrays) / nUsers);
remainder = rem((problem.NxSubarrays * problem.NySubarrays),nUsers);
totSubArrays = totSubArrays_1 .* ones(1,nUsers);
totSubArrays(1:remainder) = totSubArrays(1:remainder) + 1;

% Antennas assigned to each user (initial, prior refinement in case needed)
mySubArray = (1:1:(problem.NxSubarrays * problem.NySubarrays));
relevant_positions = cell(nUsers,1);
for valID = 1:nUsers
    partAssignation = mySubArray(1:totSubArrays(valID));
    temp = [];
    for ass = partAssignation
        temp = [temp problem.Partition{ass}];  %#ok<AGROW>
    end
    relevant_positions{valID} = temp;
    for k = partAssignation
        mySubArray(mySubArray==k) = [];  % delete assigned antennas
    end
end

% Initialize refinement parameters for user selection in case of shortage
% of antenna resources. Some users will be left out.
Delta1 = 2e-1;  % 20% tolerance = 5x max difference on number of antennas
Delta2 = 1;     % Number of devices to cut out for next iteration
Delta3 = 1.3;   % Tolerance of number of antennas
% Iterate until finding a candSet that allows for good cancellation amongst
% users. The less number of users scheduled to transmit, the more antennas
% for each, the more antennas per constraint, the least interference
% generated, the higher performance.
while refine
    % Compute number of sub-arrays to assign per user. We ensure each user
    % receives one array and locate them horizontaly
    if ceil(sqrt(problem.N_Antennas)) > nUsers
        problem.NySubarrays = nUsers;
        problem.NxSubarrays = 1;
    else
        if mod(nUsers,2)~=0;  t = factor(nUsers + 1);  % odd
        else;                 t = factor(nUsers);  % even
        end
        t = [t(1) prod(t(2:end))];
        problem.NxSubarrays = t(1);
        problem.NySubarrays = t(2);
    end
    problem.N_Subarrays = problem.NxSubarrays*problem.NySubarrays;
    problem = o_compute_antennas_per_user(problem,candSet); 
    
    % Create subarray partition
    problem = o_create_subarray_partition(problem);
    
    % Distribute partitions amongst users
    totSubArrays_1 = floor((problem.NxSubarrays * problem.NySubarrays) / nUsers);
    remainder = rem((problem.NxSubarrays * problem.NySubarrays),nUsers);
    totSubArrays = totSubArrays_1 .* ones(1,nUsers);
    totSubArrays(1:remainder) = totSubArrays(1:remainder) + 1;
    
    % Antennas assigned to each user (fixed)
    mySubArray = (1:1:(problem.NxSubarrays * problem.NySubarrays));
    relevant_positions = cell(nUsers,1);
    for valID = 1:nUsers
        partAssignation = mySubArray(1:totSubArrays(valID));
        temp = [];
        for ass = partAssignation
            temp = [temp problem.Partition{ass}];  %#ok<AGROW>
        end
        relevant_positions{valID} = temp;
        for k = partAssignation
            mySubArray(mySubArray==k) = [];  % delete assigned antennas
        end
    end
    
    % Check if the configuration is valid
    cellsz = cellfun(@length,relevant_positions,'uni',false);
    n_ant_users = cell2mat(cellsz);
    totConstraints = nUsers;
    if any(n_ant_users./totConstraints < Delta3)
        % Need to refine the solution - even traffic demands
        candIdxObjF = problem.MinObjF ./ max(problem.MinObjF);
        idx = find(candIdxObjF>=Delta1);
        candSet = candSet(idx);  % adapt
        problem.MinObjF = problem.MinObjF(idx);  % adapt
        % Need to refine the solution - Cut down number of users scheduled
        [~,idxList] = sort(problem.MinObjF);
        if length(idxList) >= 3
            idxList(end-Delta2+1:end) = [];  % take out the last N users 
        end
        candSet = candSet(idxList);  % update
        problem.MinObjF = problem.MinObjF(idxList);  % update
        % Re-arrange indexes
        [candSet,idxSort] = sort(candSet);
        problem.MinObjF = problem.MinObjF(idxSort);
    else
        % We've found a good solution
        refine = false;  % Stop the loop
    end
    
    % Update number of users
    nUsers = length(candSet);
end

% Recreate the conformal array
problem.NzPatch = problem.NxPatch;
problem.dz = problem.dx;
problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                        [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
                        'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque si
handle_ConformalArray = phased.URA([problem.NyPatch,problem.NzPatch],...
                          'Lattice','Rectangular','Element',problem.handle_Ant,...
                          'ElementSpacing',[problem.dy,problem.dz]);
problem.possible_locations = handle_ConformalArray.getElementPosition;

% Compute weights (beamforming)
W_LCMV = zeros(nUsers, problem.NxPatch*problem.NyPatch);
W_CBF = zeros(nUsers, problem.NxPatch*problem.NyPatch);
for id = 1:nUsers
    Nant_user = length(relevant_positions{id});
    % Convert antenna ID's into physical locations
    elementPos = [problem.possible_locations(1,relevant_positions{id});...
                  problem.possible_locations(2,relevant_positions{id});...
                  problem.possible_locations(3,relevant_positions{id})];
    elementPosNorm = elementPos./problem.lambda;

    % Determine which other users to null (based on limitation of LCMV)
    nullSet = (1:1:nUsers);
    maxUsersNull = Nant_user - 1;  % case Single-path
    if maxUsersNull >= (nUsers - 1)
        % We can nullify all of them
        nullSet(id) = []; % Take out current user
    else
        % We can only nullify a restricted set
        myListPrior = problem.MinObjF;  % Priorize by traffic
        myListPrior(id) = [];  % Take out current user
        [~,priorList] = sort(myListPrior);
        nullList = priorList(1:maxUsersNull);
        nullSet(id) = [];  % Take out current user
        nullSet = nullSet(nullList);
    end

    % Always locate the intented user on location 1 of the vector
    totSet = [id nullSet];
    PhiTheta = ([-problem.phiUsers(candSet(totSet)) ; -problem.thetaUsers(candSet(totSet))]);

    % Apply LCMV Beamformer for selected user
    sv = steervec(elementPosNorm,PhiTheta);
    Sn = eye(Nant_user);
    resp = zeros(length(totSet),1) + eps;
    resp(1) = 1;  % Maximum towards intended direction
    w_lcmv = lcmvweights(sv,resp,Sn);  % LCMV Beamformer method

    % Apply Convencional Beamformer for selected user
    w_cbf = cbfweights(elementPosNorm,PhiTheta(:,1));  % conventional beamformer

    % Normalize weights
    w_lcmv = (1/sqrt(w_lcmv'*w_lcmv)) * w_lcmv;
    w_cbf = (1/sqrt(w_cbf'*w_cbf)) * w_cbf;

    % Store results in global W
    W_LCMV(id,relevant_positions{id}) = w_lcmv.';
    W_CBF(id,relevant_positions{id}) = w_cbf.';
end

% Extract Directivities
conf.verbosity = 0;
[~,~,Cap_LCMV,SINRLCMV_PB]  = f_BF_results(W_LCMV,handle_ConformalArray,candSet,problem,conf,false);
[~,~,Cap_CBF,SINRCBF_PB]  = f_BF_results(W_LCMV,handle_ConformalArray,candSet,problem,conf,false);

if conf.MinObjFIsSNR
    estObj_LCMV = db2pow(SINRLCMV_PB).';
    estObj_CBF = db2pow(SINRCBF_PB).';
else
    estObj_LCMV = db2pow(Cap_LCMV).';
    estObj_CBF = db2pow(Cap_CBF).';
end

    
% EOF

