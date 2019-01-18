function InitialPopulation  = CBG_creationArrayGA(problem,conf)
% CBG_creationArrayGA - Configures the initial genome for the genetic
% algorithm. Contains information about the antenna (subarray partition)
% selection as well as the amplitudes and phases
%
% Syntax:  InitialPopulation  = o_creationArrayGA(GenomeLength, ~, options,
% problem, conf)
%
% Inputs:
%    problem - struct containint configuration in data/metaproblem_test.dat
%    conf - struct containint configuration in data/config_test.dat
%
% Outputs:
%    InitialPopulation - The initial population. An array of size
%    options.PopulationSize.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CBG_solveit

%------------- BEGIN CODE --------------
    antPerUser = problem.NmaxArray;
    availableAnt = (1:1:problem.N_Antennas);
    availableAnt1 = availableAnt;
    assignation = zeros(1,problem.N_Antennas);
    if isfield(problem,'initialW')
        % Initial gene configuration
        for id = 1:problem.nUsers
            % Select relevant indices
            antennaSelected = (problem.initialW(id,:)~=0);
            % Allocate those to designated user
            assignation(antennaSelected) = id;
        end
    else
        % Initial gene configuration
        for id = 1:problem.nUsers
            % Select antenna indices to be allocated to selected user
            antennaSelected = randsample(availableAnt1,antPerUser(id));
            % Update gene with those indices
            assignation(antennaSelected) = id;
            % Remove selected antennas from active set
            availableAnt1 = setdiff(availableAnt1,antennaSelected);
        end
    end
    myGene = [availableAnt , assignation];
    InitialPopulation = repmat(myGene,conf.PopulationSize_Data,1);
end

    
    
% EOF
