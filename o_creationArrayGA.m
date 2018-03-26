function InitialPopulation  = o_creationArrayGA(GenomeLength, ~, options, ...
    problem, conf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    InitialPopulation = zeros(options.PopulationSize,GenomeLength);
    if strcmp(conf.genStructure,'allAntennas')
        for i = 1:options.PopulationSize
            amplitudes = zeros(1,GenomeLength/2);
            pos = randperm(GenomeLength/2,problem.Nmax);
            amplitudes(pos) = rand(1,problem.Nmax)*...
                (problem.maxAmpW-problem.minAmpW) + problem.minAmpW; %always > 0
            phases = zeros(1,GenomeLength/2);
            phases(pos) = rand(1,problem.Nmax)*...
                (problem.maxPhaseW-problem.minPhaseW) + problem.minPhaseW;
            InitialPopulation(i,:) = [amplitudes,phases];
            if isempty(find(InitialPopulation(i,:),1))
                disp('creation:uyuyuy');
            end
        end
    elseif strcmp(conf.genStructure,'onlyAssigned')
        for i = 1:options.PopulationSize
            antennasSelected = randperm(problem.N_Subarrays,problem.Nmax);
            amplitudes = rand(1,problem.Nmax)*...
                (problem.maxAmpW-problem.minAmpW) + problem.minAmpW;
            phases = rand(1,problem.Nmax)*...
                (problem.maxPhaseW-problem.minPhaseW) + problem.minPhaseW;
            Fbb = rand;
            InitialPopulation(i,:) = [antennasSelected,amplitudes,phases,Fbb];
        end
    end
end

