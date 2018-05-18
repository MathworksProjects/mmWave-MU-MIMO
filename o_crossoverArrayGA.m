function xoverKids = o_crossoverArrayGA(parents,options,nvars, ...
    ~,~,thisPopulation,problem,conf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    nEliteKids = options.EliteCount;
    nXoverKids = round(options.CrossoverFraction * (options.PopulationSize - nEliteKids));
    xoverKids = zeros(nXoverKids,nvars);

    if strcmp(conf.genStructure,'allAntennas')
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
    elseif strcmp(conf.genStructure,'onlyAssigned')
        for i=1:nXoverKids
            par = randperm(length(parents),2);
            parent1 = thisPopulation(parents(par(1)),:);
            parent2 = thisPopulation(parents(par(2)),:);
            p1 = ceil(((nvars-1)/3) * rand);
            child = zeros(1,nvars);
            child(p1:((nvars-1)/3)) = parent2(p1:((nvars-1)/3));
            child(1:p1-1) = parent1(1:p1-1);
            p2 = p1+((nvars-1)/3);
            child(p2:(2*((nvars-1)/3))) = parent2(p2:(2*((nvars-1)/3)));
            child((((nvars-1)/3)+1):p2-1) = parent1((((nvars-1)/3)+1):p2-1);
            p3 = p2+((nvars-1)/3);
            child(p3:(3*((nvars-1)/3))) = parent2(p3:(3*((nvars-1)/3)));
            child((2*((nvars-1)/3)+1):p3-1) = parent1((2*((nvars-1)/3)+1):p3-1);
            if rand < 0.5
                child(end) = parent1(end);
            else
                child(end) = parent2(end);
            end
            % Avoid repeated antennas!!
            uniqueAnt = unique(child(1:((nvars-1)/3)));
            nonSelectedAntennas = setdiff(1:problem.N_Subarrays, uniqueAnt);
            newAntennas = nonSelectedAntennas(randsample(length(nonSelectedAntennas),...
                            ((nvars-1)/3)-length(uniqueAnt)));
            uniqueAnt = uniqueAnt*0;
            if length(uniqueAnt) ~= ((nvars-1)/3)
                u = 1;
                a = 1;
                for j = 1:((nvars-1)/3)
                    if isempty(find(uniqueAnt==child(j), 1))
                        uniqueAnt(u) = child(j);
                        u = u + 1;
                    else
                        child(j) = newAntennas(a);
                        a = a + 1;
                    end
                end
            end
            xoverKids(i,:) = child;
            uniqAnt = unique(child(1:((nvars-1)/3)));
            if conf.verbosity > 1 && length(uniqAnt) ~= ((nvars-1)/3)
                fprintf('crossover:uyuyuy\n');
            end
        end
    end
    
end

