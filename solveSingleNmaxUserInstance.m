function [sol,W,PRx,I] = solveSingleNmaxUserInstance(conf,problem,handle_Ant,mode)
%SOLVESINGLENMAXUSERINSTANCE Obtain the antenna allocation partial solution
%   for a single <User, Nmax> pair
%   Detailed explanation goes here
    if nargin < 6
        random_calc = false;
    elseif strcmp(mode,'random')
        random_calc = true;
    else
        random_calc = false;
    end
    
    NAntToBeAssigned = problem.NmaxArray(problem.IDUserAssigned);
    
    if NAntToBeAssigned == 0
        sol = zeros(0);
        W = zeros(1,problem.NxPatch*problem.NyPatch);
        PRx = -Inf;
        I = ones(1,problem.nUsers)*-Inf;
        return
    end
    
    X0 = generateInitialValue(NAntToBeAssigned,problem,conf);

    [X0_elem,problem] = subarrayToAntennaElements(X0,conf,problem);

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
        handle_Conf_Array = phased.ConformalArray('Element',handle_Ant,...
            'ElementPosition',[zeros(1,problem.ant_elem);...
            problem.possible_locations(2,X0_elem(1:problem.ant_elem));...
            problem.possible_locations(3,X0_elem(1:problem.ant_elem))],...
            'Taper',Taper_value);
        if conf.plotAssignmentInitialAndFinal
            plotAssignment(problem,handle_Conf_Array);
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
        else
            lb = zeros(1,length(X0));
            lb(1:NAntToBeAssigned) = 1;
            lb((1+NAntToBeAssigned):(2*NAntToBeAssigned)) = 0.1;
            lb((2*NAntToBeAssigned+1):end) = 0;%-pi;
            ub = zeros(1,length(X0));
            ub(1:NAntToBeAssigned) = problem.N_Subarrays;
            ub((1+NAntToBeAssigned):(2*NAntToBeAssigned)) = 1;
            ub((2*NAntToBeAssigned+1):end) = 0;%pi;
        end
        %% 
        % Start optimization
        PopulationSize_Data = 20;
        EliteCount_Data = 2;
        CrossoverFraction_Data = 0.7;
        MaxGenerations_Data = 400;
        MaxStallGenerations_Data = 40;
        if strcmp(conf.algorithm,'GA')
            [x,~,~,~,~,~] = CA_Position_Objective_optim_ga(conf,problem,lb,ub,PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,MaxGenerations_Data,MaxStallGenerations_Data);
        elseif strcmp(conf.algorithm,'PS')
            [x,~,~,~] = CA_Position_Objective_optim_pattern(X0,conf,problem,lb,ub,1e-1,1e-1,250);
        elseif strcmp(conf.algorithm,'PSO')
            InitialPopulation = zeros(40,length(lb));
            for i=1:40
                InitialPopulation(i,:) = generateInitialValue(NAntToBeAssigned,...
                    problem,conf);
            end
            [x,~,~,~] = CA_Position_Objective_optim_pso(X0,conf,problem,lb,ub,InitialPopulation);
        else
            x = X0;
        end
        %%
        [x_elem,problem] = subarrayToAntennaElements(x,conf,problem);
    else
        x_elem = X0_elem;
    end
    % Extracting taper avalues from the input vector
    Taper_value = x_elem(problem.ant_elem+1:problem.ant_elem*2) .* ...
        exp(1i.*x_elem(problem.ant_elem*2+1:problem.ant_elem*3));

    % Creating a Conformal Array with cosine elements.
    % Conformal array will be limited to a single plane
    handle_Conf_Array = phased.ConformalArray('Element',handle_Ant,...
        'ElementPosition',[zeros(1,problem.ant_elem);...
        problem.possible_locations(2,x_elem(1:problem.ant_elem));...
        problem.possible_locations(3,x_elem(1:problem.ant_elem))],...
        'Taper',Taper_value);

    if conf.plotAssignmentInitialAndFinal
        plotAssignment(problem,handle_Conf_Array);
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
    %%    
    sol = x_elem;
    PotRx = zeros(1,problem.nUsers);
    for u1=1:problem.nUsers
        if conf.multiPath
            for i=1:problem.maxnChannelPaths
                if problem.phiChannels(u1,i) ~= -Inf
                    PotRx(u1) = PotRx(u1) + problem.alphaChannels(u1,i)*...
                        pattern(handle_Conf_Array,problem.freq,...
                        problem.phiChannels(u1,i),...
                        problem.thetaChannels(u1,i),...
                        'Type','Power','Normalize',false);
                end
            end
        else
            PotRx(u1) = pattern(handle_Conf_Array,problem.freq,...
                    problem.phiUsers(u1),...
                    problem.thetaUsers(u1),...
                    'Type','Powerdb','Normalize',false);
        end
    end
    
    PotRx = 10*log10(PotRx);
    
    PRx = PotRx(problem.IDUserAssigned);
    PotRx(problem.IDUserAssigned) = [];
    % The interferences vector in the solution files does not contain the
    % interference inflicted to the user being analyzed [U] (it's nonsense)
    % Therefore, we need to shift the read index once we have read the
    % users with ID lower than U, and assign 0 to the interference
    % inflicted to himself.
    shift = 0;
    I = zeros(1,problem.nUsers);
    for m=1:problem.nUsers
        if m == problem.IDUserAssigned
            I(m) = 0;
            shift = -1;
        else
            I(m) = PotRx(m+shift);
        end
    end
    amplitude = x_elem(NAntToBeAssigned+1:2*NAntToBeAssigned);
    phase = x_elem(2*NAntToBeAssigned+1:end);
    W = zeros(1,problem.NxPatch*problem.NyPatch);
    W(x_elem(1:NAntToBeAssigned)) = amplitude.*cos(phase) + 1i*amplitude.*sin(phase);
end

