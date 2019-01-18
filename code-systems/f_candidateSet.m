function [combIds,combTH] = f_candidateSet(t,flows,selFlow)
% F_CANDIDATESET - Selects the user ID that require service in the current
% time slot, given the PHY-flows
%
% Syntax:  [combIds,combTH] = f_candidateSet(t,flows,selFlow)
%
% Inputs:
%    t - Scalar value that represents the slot ID in the DES
%    flows - Array of structs of length equal to the number of
%            users. For each user, each flow (belonging to each packet) is
%            characterized by the amount of bits that needs to be delivered. The
%            amount of bits are distributed uniformly across the slots until
%            reaching the deadline. Thus, the variable flows contains four features:
%            - slots:     For each flow, the slots across it.
%            - TH:        The average throughput demanded for that flow over the slots
%                         indicated in the slots field.
%            - remaining: We set the remaining field of the flow to be the value of
%                         the Payload.
%            - deadlines: The deadline of the actual packet in slot ID. 
%                         This is used in the future to set priorities.
%            - failed:    Mark whether the flow has failed to be served before the
%                         deadline (0 or 1).
%            - success:   Mark whether the flow has been succesfully served before the
%                         deadline (0 or 1).
%            - maxSlot:   Maximum deadline slot used as simulation time or Tsym in
%                         the system in case trafficType is 'dataSet'.
%    selFlow - For each user, we select a particular flow that maps with
%              the current time slot in the simulator. Only one flow per user can be
%              selected, per slot. Aggregation of upper layer flows has already
%              happened at this point.
%
% Outputs:
%    combIds - List of user IDs that require service in the current slot
%    combTH - Demanded Throughput per user in the current slot
%
% Example:
%           addpath('data','-end');
%           problem = o_read_input_problem('metaproblem_test.dat);
%           conf = o_read_config('config_test.dat');
%           [problem,~,flows] = f_configuration(conf, problem);  % Struct with configuration parameters
%           t = 1;  % first slot
%           selFlow = zeros(problem.nUsers,1);  % Inizialization
%           [combIds,combTH] = f_candidateSet(t,flows,selFlow);
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: main
%
%------------- BEGIN CODE --------------

% Initialize variables
Nusers = length(flows);  % Number of users
priority = zeros(1,Nusers);  % Priority vector
userID = (1:1:Nusers);  % User ID vector
% priority assignation
for id = 1:Nusers
    if selFlow(id)~=0
        deadline = flows(id).deadlines(selFlow(id));
        Delta = deadline - t;
        % Priorities are inversively proportional to the remaining time
        % (in slots) to deadline (This is Policy PLt)
        priority(id) = 1/Delta;
    end
end
% Sort Priorities and return candidate set
idxPrior = priority~=0;      % Index of the users with non 0 priority
prior = priority(idxPrior);  % Priorities
prior(isinf(prior)) = 1.1;   % Inf is substituted by value slightly bigger than 1.
priorID = userID(idxPrior);  % User IDs with non-0 priorities. 
                             % This is the candidate set to be
                             % scheduled in the slot
% Create combinatorial matrix
N = 2^length(prior);
binMatrix = de2bi((0:1:N-1),length(prior));
% Compute possible combinations
combIds = binMatrix .* repmat(priorID,N,1);
% Compute aggregate priority amongst users
sumVec = binMatrix*prior.';
% Replace NaN (priority=Inf when the Deadline is right now)
sumVec(isnan(sumVec)) = 0;
% Sort the aggregated priority in a descending way
[~,idxSort] = sort(sumVec,'descend');
% Apply Greedy policy to create a sorted candidate set matrix
combIds = combIds(idxSort,:);
% Create Matrix with the demanded Throughput per user per combination
combTH = zeros(size(combIds));
for idx = 1:size(combIds,1)
    idList = combIds(idx,:);
%         idList = idList(idList~=0);
    for id = idList
        if id~=0
            combTH(idx,id==idList) = flows(id).TH(selFlow(id));
        end
    end
end
% Trim output - remove last combination (empty)
combIds(end,:) = [];
combTH(end,:) = [];



% EOF