function [gene,W,PRx,I,bestScores] = o_solveSingleNmaxUserInstance(conf,problem,...
                            NAntToBeAssigned,mode)
%SOLVESINGLENMAXUSERINSTANCE Obtain the antenna allocation partial solution
%   for a single <User, Nmax> pair
%   Detailed explanation goes here
    if nargin < 4
        random_calc = false;
    elseif strcmp(mode,'random')
        random_calc = true;
    else
        random_calc = false;
    end
    
    if conf.verbosity >= 1
        fprintf('We are going to assign %d\n',NAntToBeAssigned);
    end
    
    if NAntToBeAssigned == 0
        gene = zeros(0);
        W = zeros(1,problem.NxPatch*problem.NyPatch);
        PRx = -Inf;
        I = ones(1,problem.nUsers)*-Inf;
        return
    end
    
    X0 = o_generateInitialValue(NAntToBeAssigned,problem,conf);

    [X0_elem,problem] = o_subarrayToAntennaElements(X0,conf,problem);

    problem.Nmax = NAntToBeAssigned;
    
    %% 
    if ~random_calc
        % Configure initial conformal array with preliminary weigths and
        % phases
        [handle_Conf_Array,~,~,~] = o_geneToAssignment(X0_elem,problem,conf);

        if conf.plotAssignmentInitialAndFinal
            o_plotAssignment_mod(problem, handle_Conf_Array);
        end
        %% 
        % Define lower and upper bounds
        if strcmp(conf.genStructure, 'nchoosek')
            lb = zeros(1,length(X0));
            lb(1) = 1;
            lb(2:(1+NAntToBeAssigned)) = 0.1;
            lb(((NAntToBeAssigned+1)+1):end) = 0;%-pi;
            ub = zeros(1,length(X0));
            ub(1) = nchoosek(problem.N_Subarrays,NAntToBeAssigned);
            ub(2:(1+NAntToBeAssigned)) = 1;
            ub(((NAntToBeAssigned+1)+1):end) = 0;%pi;
        elseif strcmp(conf.genStructure, 'allAntennas')
            problem.maxAmpW = 1;
            problem.fixedPointScale = (problem.maxAmpW-0)/(2^conf.NbitsAmplitude);
            lb = zeros(1,length(X0));
            lb(1,problem.N_Subarrays) = 0;
            lb(1+problem.N_Subarrays,end) = 0;
            problem.minAmpW = problem.fixedPointScale;
            problem.minPhaseW = 0;%-pi;
            ub = zeros(1,length(X0));
            ub(1,problem.N_Subarrays) = 1;
            ub(1+problem.N_Subarrays,end) = 0;
            problem.maxPhaseW = 0;%pi;
        else % 'onlyAssigned'
            problem.maxAmpW = 1; % Constant modulus
            problem.fixedPointScale = (problem.maxAmpW-0)/(2^conf.NbitsAmplitude);
            problem.minAmpW = problem.fixedPointScale;  % Constant modulus
            problem.minPhaseW = 0; % Phase is computed
            problem.maxPhaseW = 0; % Phase is computed
            lb = zeros(1,length(X0));
            lb(1:NAntToBeAssigned) = 1;
            lb((1+NAntToBeAssigned):(2*NAntToBeAssigned)) = problem.minAmpW;
            lb((2*NAntToBeAssigned+1):(3*NAntToBeAssigned)) = problem.minPhaseW;
            lb((3*NAntToBeAssigned+1):end) = problem.minAmpW; %Fbb
            ub = zeros(1,length(X0));
            ub(1:NAntToBeAssigned) = problem.N_Subarrays;
            ub((1+NAntToBeAssigned):(2*NAntToBeAssigned)) = problem.maxAmpW;
            ub((2*NAntToBeAssigned+1):(3*NAntToBeAssigned)) = problem.maxPhaseW;
            ub((3*NAntToBeAssigned+1):end) = 1; %Fbb
        end
        %% 
        % Start optimization
        PopulationSize_Data = conf.PopulationSize_Data;
        EliteCount_Data = conf.EliteCount_Data;
        CrossoverFraction_Data = conf.CrossoverFraction_Data;
        Maxgenerations_Data = conf.Maxgenerations_Data;
        MaxStallgenerations_Data = conf.MaxStallgenerations_Data;
        FunctionTolerance_Data = conf.FunctionTolerance_Data;
        if strcmp(conf.algorithm,'GA')
            [x,~,bestScores,~,~,~,~] = o_CA_Position_Objective_optim_ga(conf,problem,lb,ub,...
                PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,...
                Maxgenerations_Data,MaxStallgenerations_Data,FunctionTolerance_Data);
        elseif strcmp(conf.algorithm,'GA-rnd')
            [x,~,~,~,~,~,~] = o_CA_Position_Objective_optim_ga(conf,problem,lb,ub,...
                PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,...
                Maxgenerations_Data,MaxStallgenerations_Data,FunctionTolerance_Data);
            x = apply_random_modification_assignment(x,conf.GArndImpact,...
                conf.genStructure,problem.N_Subarrays,conf.GArndmodifyAmpl,...
                problem.maxAmpW);
            bestScores = [];
        elseif strcmp(conf.algorithm,'PS')
            [x,~,~,~] = o_CA_Position_Objective_optim_pattern(X0,conf,problem,lb,ub,1e-1,1e-1,250);
        elseif strcmp(conf.algorithm,'PSO')
            InitialPopulation = zeros(40,length(lb));
            for i=1:40
                InitialPopulation(i,:) = o_generateInitialValue(NAntToBeAssigned,...
                    problem,conf);
            end
            [x,~,~,~] = o_CA_Position_Objective_optim_pso(X0,conf,problem,lb,ub,InitialPopulation);
        elseif strcmp(conf.algorithm,'ES')
            [x,bestScores] = o_CA_Position_Objective_optim_exhaustive(conf,problem,lb,ub);
        else    
            x = X0;
        end
        %%
        [gene,problem] = o_subarrayToAntennaElements(x,conf,problem);
    else
        gene = X0_elem;
        bestScores = [];
    end
    %%
    [handle_Conf_Array,W,PRx,I] = o_geneToAssignment(gene,problem,conf);
    
    %%
    if conf.plotAssignmentInitialAndFinal
        o_plotAssignment(problem,handle_Conf_Array);
    end

    %%
    if conf.plotAssignmentInitialAndFinal
        % plot of initial values, lb,ub and optimized values.

        f = figure;
        plot(X0,'DisplayName','X0');hold on;plot(x,'DisplayName','x');plot(ub,'DisplayName','ub');plot(lb,'DisplayName','lb');hold off;
        disp('Program paused, press any key to continue')
        pause
        if isvalid(f)
            close(f)
        end
    end
end

