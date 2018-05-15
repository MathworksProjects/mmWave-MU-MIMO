function [gene,W,PRx,I] = o_solveSingleNmaxUserInstance(conf,problem,...
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
        % constructing a conformal array

        % Array geometry and pattern for initial values.

        % Extracting taper avalues from the input vector
        Taper_value = X0_elem(problem.ant_elem+1:problem.ant_elem*2) .* ...
            exp(1i.*X0_elem(problem.ant_elem*2+1:problem.ant_elem*3));

        % Creating a Conformal Array with cosine elements.
        % Conformal array will be limited to a single plane
        handle_Conf_Array = phased.ConformalArray('Element',problem.handle_Ant,...
            'ElementPosition',[zeros(1,problem.ant_elem);...
            problem.possible_locations(2,X0_elem(1:problem.ant_elem));...
            problem.possible_locations(3,X0_elem(1:problem.ant_elem))],...
            'Taper',Taper_value);
        if conf.plotAssignmentInitialAndFinal
            o_plotAssignment(problem,handle_Conf_Array);
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
            lb = zeros(1,length(X0));
            lb(1,problem.N_Subarrays) = 0;
            lb(1+problem.N_Subarrays,end) = 0;
            problem.minAmpW = 0.1;%0.01
            problem.minPhaseW = 0;%-pi;
            ub = zeros(1,length(X0));
            ub(1,problem.N_Subarrays) = 1;
            ub(1+problem.N_Subarrays,end) = 0;
            problem.maxAmpW = 1;
            problem.maxPhaseW = 0;%pi;
        else % 'onlyAssigned'
            problem.minAmpW = 1;  % Constant modulus
            problem.minPhaseW = 0; % Phase is computed
            problem.maxAmpW = 1; % Constant modulus
            problem.maxPhaseW = 0; % Phase is computed
            lb = zeros(1,length(X0));
            lb(1:NAntToBeAssigned) = 1;
            lb((1+NAntToBeAssigned):(2*NAntToBeAssigned)) = problem.minAmpW;
            lb((2*NAntToBeAssigned+1):(3*NAntToBeAssigned)) = problem.minPhaseW;
            lb((3*NAntToBeAssigned+1):end) = 0; %Fbb
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
            [x,~,~,~,~,~] = o_CA_Position_Objective_optim_ga(conf,problem,lb,ub,...
                PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,...
                Maxgenerations_Data,MaxStallgenerations_Data,FunctionTolerance_Data);
        elseif strcmp(conf.algorithm,'PS')
            [x,~,~,~] = o_CA_Position_Objective_optim_pattern(X0,conf,problem,lb,ub,1e-1,1e-1,250);
        elseif strcmp(conf.algorithm,'PSO')
            InitialPopulation = zeros(40,length(lb));
            for i=1:40
                InitialPopulation(i,:) = o_generateInitialValue(NAntToBeAssigned,...
                    problem,conf);
            end
            [x,~,~,~] = o_CA_Position_Objective_optim_pso(X0,conf,problem,lb,ub,InitialPopulation);
        else
            x = X0;
        end
        %%
        [gene,problem] = o_subarrayToAntennaElements(x,conf,problem);
    else
        gene = X0_elem;
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

