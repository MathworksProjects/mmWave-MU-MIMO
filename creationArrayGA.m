function InitialPopulation  = creationArrayGA(GenomeLength, ~, options, problem)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    InitialPopulation = zeros(options.PopulationSize,GenomeLength);
    for i = 1:options.PopulationSize
        weights = zeros(1,GenomeLength/2);
        pos = randperm(GenomeLength/2,problem.Nmax);
        weights(pos) = rand(1,problem.Nmax)*(problem.maxAmpW-problem.minAmpW) + problem.minAmpW; % always > 0
        phases = zeros(1,GenomeLength/2);
        phases(pos) = rand(1,problem.Nmax)*(problem.maxPhaseW-problem.minPhaseW) + problem.minPhaseW; % always > 0
        InitialPopulation(i,:) = [weights,phases];
        if isempty(find(InitialPopulation(i,:),1))
            disp('creation:uyuyuy');
        end
    end
end

