function [problem] = o_compute_antennas_per_user(problem,usersToBeAssigned)
%O_COMPUTE_ANTENNAS_PER_USER - The function computes the number of antennas
%to be allocated to each user based on either inheritted beamforming
%methods (CBF or LCMV) or the traffic demands (proposed mechanism)
%
% Syntax:  [problem] = o_compute_antennas_per_user(problem,usersToBeAssigned)
%
% Inputs:
%    problem - struct containing configuration in data/metaproblem_test.dat 
%    usersToBeAssigned - Vector containing the users IDs
%
% Outputs:
%    problem - New variable is created and appended into the struct:
%    NmaxArray, containg the number of antennas to be allocated
%
% Example: 
%    usersToBeAssigned = [1 3 4 6];
%    [problem] = o_compute_antennas_per_user(problem,usersToBeAssigned)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: f_heuristics

%------------- BEGIN CODE --------------

% Variable to store the total number of antennas to assign to each user
problem.NmaxArray = zeros(1,problem.nUsers);

if isfield(problem,'initialW')
    % Heuristics first assign antennas based onnthe configuration mandated
    % by conventional beamforming methods (either CBF or LCMV)
    for u = usersToBeAssigned
        % Need to index
        idx = (u==usersToBeAssigned);
        problem.NmaxArray(idx) = length(find(problem.initialW(idx,:)~=0));
    end
else
    % Heuristics do not inherit from any conventional beamforming
    % mechanis. Thus, we can allocate antennas more effectivelly in order
    % to meet the application requirements
    minSub = 0;
    [~,orderedIndices] = sort(problem.MinObjF,'descend');
    for u = orderedIndices
        if problem.MinObjF(u) > 0 && minSub < problem.N_Subarrays
            problem.NmaxArray(u) = 1; % At least 1
            minSub = minSub + 1;
        end
    end
    sum_MinObjF = sum(problem.MinObjF);
    if sum_MinObjF > 0
        thr_ant_converter = (problem.N_Subarrays-minSub)/sum_MinObjF;
    else
        thr_ant_converter = (problem.N_Subarrays-minSub)/numel(usersToBeAssigned);
    end
    for u = usersToBeAssigned
        % Need to index
        idx = u==usersToBeAssigned;
        % Now NmaxArray is the Nmax value chosen for each user
        problem.NmaxArray(idx) = max(0,problem.NmaxArray(idx) + ...
            floor(problem.MinObjF(idx)*thr_ant_converter));
    end
    % The previous division might have left some antennas / subarrays
    % unassigned, due to rounding errors. Let's assign the remaining 
    % antennas / subarrays to the user most in need
    sumAnt = 0;
    for u = usersToBeAssigned; sumAnt = sumAnt + problem.NmaxArray(u); end
    if sumAnt < problem.N_Subarrays % Se podría equilibrar el reparto...
        [~,i] = max(problem.MinObjF);
        userToAssignExtraAnt = usersToBeAssigned(i);
        problem.NmaxArray(userToAssignExtraAnt) = ...
            problem.NmaxArray(userToAssignExtraAnt) + ...
            (problem.N_Subarrays-sumAnt);
    end
end

    
% EOF

