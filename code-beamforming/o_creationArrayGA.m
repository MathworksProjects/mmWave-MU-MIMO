function InitialPopulation  = o_creationArrayGA(GenomeLength, ~, options, ...
    problem, conf)
% O_CREATIONARRAYGA - Configures the initial genome for the genetic
% algorithm. Contains information about the antenna (subarray partition)
% selection as well as the amplitudes and phases
%
% Syntax:  InitialPopulation  = o_creationArrayGA(GenomeLength, ~, options, problem, conf)
%
% Inputs:
%    GenomeLength - The population size
%
% Outputs:
%    InitialPopulation - The initial population. An array of size
%    options.PopulationSize.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

%------------- BEGIN CODE --------------

InitialPopulation = zeros(options.PopulationSize,GenomeLength);
if strcmp(conf.genStructure,'allAntennas')
    for i = 1:options.PopulationSize
        % Amplitude Selection
        amplitudes = zeros(1,GenomeLength/2);
        pos = randperm(GenomeLength/2,problem.Nmax);
        proposedAmplitudes = rand(1,problem.Nmax)*...
            (problem.maxAmpW-problem.minAmpW) + problem.minAmpW; %always > 0
        if conf.NbitsAmplitude ~= 0
            fixedPointAmplitudes = floor(proposedAmplitudes/problem.fixedPointScale)*problem.fixedPointScale;
        else
            fixedPointAmplitudes = proposedAmplitudes;
        end
        amplitudes(pos) = fixedPointAmplitudes;
        % Phase Selection
        phases = zeros(1,GenomeLength/2);
        phases(pos) = rand(1,problem.Nmax)*...
            (problem.maxPhaseW-problem.minPhaseW) + problem.minPhaseW;
        % Final initial configuration
        InitialPopulation(i,:) = [amplitudes,phases];
        if isempty(find(InitialPopulation(i,:),1))
            disp('creation:uyuyuy');
        end
    end
elseif isfield(problem,'initialW') && problem.IDUserAssigned==1
    % We retrieve the initial configuration from the conventional
    % Beamformers. Since the antenna allocation changes in heuristics, the
    % baseline only applies for the first user
    for i = 1:options.PopulationSize
        % Antenna Selection
        antAlloc   = find(problem.initialW(problem.IDUserAssigned,:)~=0);
        partition  = o_antennas_to_subarrays(antAlloc,problem.Partition);
        % Amplitude Selection
        amplitudes = abs(problem.initialW(problem.IDUserAssigned,antAlloc));
        % Phase Selection
        phases     = angle(problem.initialW(problem.IDUserAssigned,antAlloc));    
        % Baseband Fbb
        Fbb = 1;  % Already accounted for in the other terms
        % Final initial configuration
        InitialPopulation(i,:) = [partition,amplitudes,phases,Fbb];
    end
elseif strcmp(conf.genStructure,'onlyAssigned')
    for i = 1:options.PopulationSize
        % Antenna Selection
        antennasSelected = randperm(problem.N_Subarrays,problem.Nmax);
        % Amplitude Selection
        proposedAmplitudes = rand(1,problem.Nmax)*...
            (problem.maxAmpW-problem.minAmpW) + problem.minAmpW; %always > 0
        if conf.NbitsAmplitude ~= 0
            fixedPointAmplitudes = floor(proposedAmplitudes/problem.fixedPointScale)*problem.fixedPointScale;
        else
            fixedPointAmplitudes = proposedAmplitudes;
        end
        amplitudes = fixedPointAmplitudes;
        % Phase Selection
        phases = rand(1,problem.Nmax)*...
            (problem.maxPhaseW-problem.minPhaseW) + problem.minPhaseW;
        % Baseband Fbb
        proposedFbb = rand; %always > 0
        if conf.NbitsAmplitude ~= 0
            fixedPointFbb = floor(proposedFbb/problem.fixedPointScale)*problem.fixedPointScale;
        else
            fixedPointFbb = proposedFbb;
        end
        Fbb = fixedPointFbb;
        % Final initial configuration
        InitialPopulation(i,:) = [antennasSelected,amplitudes,phases,Fbb];
    end
end

    
    
% EOF
