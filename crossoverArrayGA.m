function xoverKids = crossoverArrayGA(parents,options,nvars, ...
    ~,~,thisPopulation,problem)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    nEliteKids = options.EliteCount;
    nXoverKids = round(options.CrossoverFraction * (options.PopulationSize - nEliteKids));
    xoverKids = zeros(nXoverKids,nvars);

    for i=1:nXoverKids
        par = randperm(length(parents),2);
        parent1 = thisPopulation(parents(par(1)),:);
        parent2 = thisPopulation(parents(par(2)),:);
        p1 = ceil((nvars/2) * rand);
        child = zeros(1,nvars);
        child(p1:nvars/2) = parent2(p1:nvars/2);
        child(1:p1-1) = parent1(1:p1-1);
        p2 = p1+nvars/2;
        child(p2:end) = parent2(p2:end);
        child((nvars/2+1):p2-1) = parent1((nvars/2+1):p2-1);
        % Regular a nmax!!!
        pos = find(child(1:nvars/2));
        if isempty(pos)
            disp('crossover:uyuyuy');
        end
        ant = randperm(numel(pos));
        ant = pos(ant);
        child(ant(problem.Nmax+1:end)) = 0; % Only Nmax antennas!!
        child(ant(problem.Nmax+1:end)+(nvars/2)) = 0; % Only Nmax antennas!!
        xoverKids(i,:) = child;
    end
end

