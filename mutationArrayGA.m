function mutationChildren  = mutationArrayGA(parents, options, nvars, ...
~, ~, ~, thisPopulation, problem)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    nEliteKids = options.EliteCount;
    nXoverKids = round(options.CrossoverFraction * (options.PopulationSize - nEliteKids));
    nMutateKids = options.PopulationSize - nEliteKids - nXoverKids;
    mutationChildren = zeros(nMutateKids,nvars);
    for i=1:length(parents)
        parent = thisPopulation(parents(i),:); % Normally thisPopulation(parents(i),:)
        % Swap two weights
        p = ceil(nvars/2 * rand(1,2));
        child = parent;
        child(p(1)) = parent(p(2));
        child(p(2)) = parent(p(1));
        % Randomly modify min(one tenth,1) of the weights
        mutations = max(ceil(problem.Nmax/2/10),1);
        m = ceil(nvars/2 * rand(1,mutations));
        for j = 1:mutations
            child(m(j)) = min(problem.maxAmpW,max(problem.minAmpW,child(m(j)) + ...
                normrnd(0,0.7*(problem.maxAmpW-problem.minAmpW))));
            if(child(m(j))==0) || child(m(j))>1
                disp('mutation:uyuyuy');
            end
        end
        pos = find(child(1:nvars/2));
        if isempty(pos)
            disp('mutation:uyuyuy');
        end
        ant = randperm(numel(pos));
        ant = pos(ant);
        child(ant(problem.Nmax+1:end)) = 0; % Only Nmax antennas!!
        % Now the same with phases
        p = min(nvars,round(nvars/2 * rand(1,2)) + nvars/2 + 1);
        child(p(1)) = parent(p(2));
        child(p(2)) = parent(p(1));
        % Randomly modify min(one tenth,1) of the weights
        for j = 1:mutations
            child(m(j) + nvars/2) = min(problem.maxPhaseW, max(problem.minPhaseW,child(m(j) + ...
                nvars/2)+normrnd(0,(problem.maxPhaseW-problem.minPhaseW)/16)));
            if(child(m(j) + nvars/2)<-pi)
                disp('mutation:uyuyuy');
            end
        end
        child(ant(problem.Nmax+1:end)+(nvars/2)) = 0; % Only Nmax antennas!!
        mutationChildren(i,:) = child; % Normally mutationChildren(i,:)
    end
end

