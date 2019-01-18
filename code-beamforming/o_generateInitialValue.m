function InitialValue = o_generateInitialValue(NAntToBeAssigned,problem,conf)
% O_GENERATEINITIALVALUE - Generate initial values for Beamforming weights
% in Heuristics. The InitialValue contains information about the antennas
% being allocated to the user, the amplitudes and phases.
%
% Syntax:  InitialValue = o_generateInitialValue(NAntToBeAssigned,problem,conf)
%
% Inputs:
%    NAntToBeAssigned - Number of antennas to assign to selected user
%    problem - struct containing configuration in data/metaproblem_test.dat
%    conf - Struct containing configuration in data/config_test.dat
%
% Outputs:
%    InitialValue - Array of size NAntToBeAssigned x 3 + 1, where the first
%    set of NAntToBeAssigned contain the location of the antennas in the
%    array (indices) and the second and third contain the phases of the
%    beamformer for each antenna respectively. 
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
%   NAntToBeAssigned = problem.N_Antennas / problem.nUsers;
%   InitialValue = o_generateInitialValue(NAntToBeAssigned,problem,conf)
%
% Other m-files required: f_configuration
% Subfunctions: none
% MAT-files required: none
%
% See also: f_heuristics,  o_solveSingleNmaxUserInstance, o_subarrayToAntennaElements

%------------- BEGIN CODE --------------

if strcmp(conf.genStructure, 'nchoosek') %CHECK
    InitialValue = [randi([1 nchoosek(problem.N_Subarrays,...
                    NAntToBeAssigned)],1,1),...
                    ones(1,NAntToBeAssigned),...
                    zeros(1,NAntToBeAssigned)];
elseif strcmp(conf.genStructure, 'allAntennas') % CHECK
    positions = randperm(problem.N_Subarrays,NAntToBeAssigned);
    weights_amp = zeros(1,problem.N_Subarrays);
    weights_amp(positions) = 1;
    InitialValue = [weights_amp,...
                    zeros(1,problem.N_Subarrays)];
elseif isfield(problem,'initialW') && problem.IDUserAssigned==1
    % We retrieve the initial configuration from the conventional
    % Beamformers. Since the antenna allocation changes in heuristics, the
    % baseline only applies for the first user
    antAlloc     = find(problem.initialW(problem.IDUserAssigned,:)~=0);
    antAllocPart = o_antennas_to_subarrays(antAlloc,problem.Partition);
    ampAlloc = zeros(1,length(antAllocPart));
    phaseAlloc = zeros(1,length(antAllocPart));
    for partID = antAllocPart
        antennasID = problem.Partition{partID};
        subpos = problem.possible_locations(:,antennasID);
        wT_analog = exp(1i*angle(steervec(subpos/problem.lambda,...
                        [problem.phiUsers(problem.IDUserAssigned);...
                        problem.thetaUsers(problem.IDUserAssigned)],...
                        conf.NumPhaseShifterBits)));
        Taper = problem.initialW(problem.IDUserAssigned,antennasID).';  % conjugate transpose
        temp  = Taper./wT_analog;
        ampAlloc(partID==antAllocPart) = abs(temp(1));  % supposedly, all values in temp are the same
        phaseAlloc(partID==antAllocPart) = angle(temp(1));  % supposedly, all values in temp are the same
    end
    Fbb = 1;
    InitialValue = [antAllocPart, ...
                    ampAlloc, ...
                    phaseAlloc, ...
                    Fbb];
else
    % select antennas(or subarray partitions) randomly amongst the
    % available set
    antAlloc     = randperm(problem.N_Subarrays,NAntToBeAssigned);
    ampAlloc     = ones(1,NAntToBeAssigned);
    phaseAlloc   = zeros(1,NAntToBeAssigned);
    InitialValue = [antAlloc, ...
                    ampAlloc, ...
                    phaseAlloc, ...
                    rand];
end
    
% EOF