%%
function [x,fval,bestScores,exitflag,output,population,score] = ...
                      CBG_CA_Position_Objective_optim_ga(conf,problem,InitialPopulation)
    % Store best results in a global variable
    global bestS;
    % Start with the default options
    options = optimoptions('ga');
    % Modify options setting
    options = optimoptions(options,'PopulationSize',      conf.PopulationSize_Data);
    options = optimoptions(options,'EliteCount',          conf.EliteCount_Data);
    options = optimoptions(options,'CrossoverFraction',   conf.CrossoverFraction_Data);
    options = optimoptions(options,'MaxGenerations',      conf.Maxgenerations_Data);
    options = optimoptions(options,'FunctionTolerance',   conf.FunctionTolerance_Data);
    options = optimoptions(options,'MaxStallGenerations', conf.MaxStallgenerations_Data);
	options = optimoptions(options,'OutputFcn',           @gaoutfun);
    options = optimoptions(options,'Display',             'off');
    options = optimoptions(options,'UseParallel',         true);
    if conf.verbosity >= 1
        options = optimoptions(options,'PlotFcn',         {@gaplotbestf,@gaplotstopping});
    end
    options = optimoptions(options,'CrossoverFcn',@(p,o,n,z1,z2,t)p_crossoverArrayGA(p,o,n,z1,z2,t,problem));
    options = optimoptions(options,'MutationFcn',@(p,o,n,z1,z2,z3,t)p_mutationArrayGA(p,o,n,z1,z2,z3,t,problem));
    options = optimoptions(options,'CreationFcn',@(g,f,o)CBG_creationArrayGA(g,f,o,problem,conf));
    options.InitialPopulation = InitialPopulation;
    InitialScores = CBG_Position_Objective_optim_cost_singlepath(InitialPopulation(1,:), conf, problem);
    options.InitialScores = repmat(InitialScores,conf.PopulationSize_Data,1);
    nvar = 2*sum(problem.NmaxArray);
    % Call optimization function
    [x,fval,exitflag,output,population,score] = ...
                           ga(@(x)CBG_Position_Objective_optim_cost_singlepath(x,conf,problem),...
                           nvar,[],[],[],[],[],[],[],[],options);
    % Retrieve the best scores
    bestScores = bestS;
end

%%
function xoverKids = p_crossoverArrayGA(parents,options,nvars,~,~,thisPopulation,problem)
    % Single point crossover - one crossover point is selected, till this
    % point the permutation is copied from the first parent, then the
    % second parent is scanned and if the number is not yet in the
    % offspring it is added
    
    nEliteKids = options.EliteCount;
    nXoverKids = round(options.CrossoverFraction * (options.PopulationSize - nEliteKids));
    xoverKids = zeros(nXoverKids,nvars);
    
    nAntennas = sum(problem.NmaxArray);

    for i=1:nXoverKids
        % Select parent to perform crossovers (always pairs)
        par = randperm(length(parents),2);
        parent1 = thisPopulation(parents(par(1)),:);
        parent2 = thisPopulation(parents(par(2)),:);
        
        % Number of antennas to be under each parent's control
        nAntP1 = ceil(nAntennas/2);
        nAntP2 = floor(nAntennas/2);
        
        % Initialize child
        child = [(1:1:nAntennas) , zeros(1,nAntennas)];
        
        % Parent 1 has priority
        child(nAntennas+1:nAntennas+nAntP1) = parent1(nAntennas+1:nAntennas+nAntP1);
        
        % Parent 2 sees what antennas have not been yet allocated and
        % locates them in its subset
        usrsAllocP1 = parent1(nAntennas+1:nAntennas+nAntP1);
        totUsers = parent1(nAntennas+1:end);
        
        % Create list of possibilities for parent 2
        for k = usrsAllocP1
            indxs = find(totUsers==k);
            totUsers(indxs(1)) = [];
        end
        
        % Start assigning possibilities to the child
        for idx = 1:nAntP2
            overallIndex = nAntennas+nAntP1+idx;
            [~,index] = find(parent2(overallIndex)==totUsers);
            if ~isempty(index)
                child(overallIndex) = parent2(overallIndex);  % Assign
                totUsers(index(1)) = [];  % Remove it from possibles
            end
        end
        
        % Assign the empty fields in child (Values in parent2 that could
        % not be assigned because parent 1 already assigned them)
        [~,idx] = find(child==0);
        if ~isempty(idx)
            % Shuffle list
            newTotUsers = totUsers(randperm(length(totUsers)));
            child(idx) = newTotUsers;
        end
        
       % Final Assignation
        xoverKids(i,:) = child;
    end
end

%%
function mutationChildren  = p_mutationArrayGA(parents,options,nvars,~,~,~,thisPopulation,problem)
    % Order changing - two numbers are selected and exchanged
    
    nEliteKids = options.EliteCount;
    nXoverKids = round(options.CrossoverFraction * (options.PopulationSize - nEliteKids));
    nMutateKids = options.PopulationSize - nEliteKids - nXoverKids;
    mutationChildren = zeros(nMutateKids,nvars);
    
    nAntennas = sum(problem.NmaxArray);
    
    % Antenna swapping array
    swap = ceil(nAntennas * rand(nAntennas,2));
    
    for i=1:length(parents)
        % Initialize before making changes
        parent = thisPopulation(parents(i),:);
        child = parent;
        
        % Swap two assignations
        child(nAntennas + swap(i,1)) = parent(nAntennas + swap(i,2));
        child(nAntennas + swap(i,2)) = parent(nAntennas + swap(i,1));
        
        % Final Assignation
        mutationChildren(i,:) = child; % Normally mutationChildren(i,:)
    end
end