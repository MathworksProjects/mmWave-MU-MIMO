function [W_LCMV,W_CBF,handle_ConformalArray,estObj_LCMV,estObj_CBF] = f_conventionalBF(problem,conf,candSet)
%
% Syntax:  [W_LCMV,W_CBF,handle_ConformalArray,estObj_LCMV,estObj_CBF] = ...
%                                    f_conventionalBF(problem,conf,candSet)
%
% Inputs:
%    problem - struct containint configuration in data/metaproblem_test.dat
%    conf - Struct containing configuration in data/config_test.dat
%    candSet - Vector containing the users ID being considered in the 
%              current slot
%
% Outputs:
%    W_LCMV - Matrix [nUser x nAntennas] containg the weights for LCMV
%    W_CBF - Matrix [nUser x nAntennas] containg the weights for CBF
%    handle_ConformalArray - Initial conformal array
%    estObj_LCMV - Onjective function value for LCMV. Could be either SINR
%                  or Capacity (controlled by conf.MinObjFIsSNR)
%    estObj_CBF - Onjective function value for CBF. Could be either SINR
%                 or Capacity (controlled by conf.MinObjFIsSNR)
%
% Example: 
%   problem = o_read_input_problem('data/metaproblem_test.dat');
%   conf = o_read_config('data/config_test.dat');
%   conf.verbosity = 1;  % To visualize metrics on command line
%   problem.nUsers = 5;  % Fix number of users manually for example
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


% Restrict sub-arrays to Localized for LCMV
problem.arrayRestriction = 'Localized';
% Compute number of sub-arrays to assign per user. We ensure each user
% receives one array and locate them horizontaly
if ceil(sqrt(problem.N_Antennas)) > problem.nUsers
    problem.NySubarrays = problem.nUsers;
    problem.NxSubarrays = 1;
else
    if mod(problem.nUsers,2)~=0;  t = factor(problem.nUsers + 1);  % odd
    else;                         t = factor(problem.nUsers);  % even
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
totSubArrays_1 = floor((problem.NxSubarrays * problem.NySubarrays) / problem.nUsers);
remainder = rem((problem.NxSubarrays * problem.NySubarrays),problem.nUsers);
totSubArrays = totSubArrays_1 .* ones(1,problem.nUsers);
totSubArrays(1:remainder) = totSubArrays(1:remainder) + 1;
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

% Antennas assigned to each user (fixed)
mySubArray = (1:1:(problem.NxSubarrays * problem.NySubarrays));
relevant_positions = cell(problem.nUsers,1);
for valID = 1:problem.nUsers
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

% Compute weights (beamforming)
PhiTheta = ([-problem.phiUsers ; -problem.thetaUsers]);
W_LCMV = zeros(problem.nUsers, problem.NxPatch*problem.NyPatch);
W_CBF = zeros(problem.nUsers, problem.NxPatch*problem.NyPatch);
for id = 1:problem.nUsers
    Nant_user = length(relevant_positions{id});
    % Convert antenna ID's into physical locations
    elementPos = [problem.possible_locations(1,relevant_positions{id});...
                  problem.possible_locations(2,relevant_positions{id});...
                  problem.possible_locations(3,relevant_positions{id})];
    elementPosNorm = elementPos./problem.lambda;

    % Apply LCMV Beamformer for selected user
    sv = steervec(elementPosNorm,PhiTheta);
    Sn = eye(Nant_user);
    resp = zeros(problem.nUsers,1) + eps;
    resp(id) = 1;  % Maximum restricted to limit (33dB)
    w_lcmv = lcmvweights(sv,resp,Sn);  % LCMV Beamformer method

    % Apply Convencional Beamformer for selected user
    w_cbf = cbfweights(elementPosNorm,PhiTheta(:,id));  % conventional beamformer

    % Normalize weights
    w_lcmv = (1/sqrt(w_lcmv'*w_lcmv)) * w_lcmv;
    w_cbf = (1/sqrt(w_cbf'*w_cbf)) * w_cbf;

    % Store results in global W
    W_LCMV(id,relevant_positions{id}) = w_lcmv.';
    W_CBF(id,relevant_positions{id}) = w_cbf.';
end

% Extract Directivities
conf.verbosity = 0;
[~,~,Cap_LCMV,SINRLCMV_PB]  = f_BF_results(W_LCMV,handle_ConformalArray,problem,conf,false);
[~,~,Cap_CBF,SINRCBF_PB]  = f_BF_results(W_LCMV,handle_ConformalArray,problem,conf,false);

if conf.MinObjFIsSNR
    estObj_LCMV = db2pow(SINRLCMV_PB).';
    estObj_CBF = db2pow(SINRCBF_PB).';
else
    estObj_LCMV = db2pow(Cap_LCMV).';
    estObj_CBF = db2pow(Cap_CBF).';
end

    
    
% EOF