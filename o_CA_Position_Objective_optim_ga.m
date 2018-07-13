function [x,fval,bestScores,exitflag,output,population,score] = o_CA_Position_Objective_optim_ga(conf,...
    problem,lb,ub,PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,...
    MaxGenerations_Data,MaxStallGenerations_Data,FunctionTolerance_Data)
    %% This is an auto generated MATLAB file from Optimization Tool.
    % Copyright 2017  The MathWorks, Inc.
    % global variable to save all the best scores
    global bestS;
    %% Start with the default options
    options = optimoptions('ga');
    %% Modify options setting
    options = optimoptions(options,'PopulationSize', PopulationSize_Data);
    options = optimoptions(options,'EliteCount', EliteCount_Data);
    options = optimoptions(options,'CrossoverFraction', CrossoverFraction_Data);
    options = optimoptions(options,'OutputFcn',@gaoutfun);
    options = optimoptions(options,'MaxGenerations', MaxGenerations_Data);
    options = optimoptions(options,'FunctionTolerance', FunctionTolerance_Data);
    options = optimoptions(options,'MaxStallGenerations', MaxStallGenerations_Data);
    options = optimoptions(options,'Display', 'off');
    if conf.verbosity >= 1
        options = optimoptions(options,'PlotFcn', {@gaplotbestf,@gaplotstopping});
    end
    options = optimoptions(options,'UseParallel', true);
    nvar = numel(lb); % == numel(ub) == num. variables in the gen
    if strcmp(conf.genStructure, 'nchoosek')
        integervars = 1;
        if conf.multiPath
            [x,fval,exitflag,output,population,score] = ...
            ga(@(x)o_Position_Objective_optim_cost_multipath(x,conf,problem),...
            nvar,[],[],[],[],lb,ub,@(x)o_all_NonLinear_Constraints(x,conf,problem),...
            integervars,options);
        else
            [x,fval,exitflag,output,population,score] = ...
            ga(@(x)o_Position_Objective_optim_cost_singlepath(x,conf,problem),...
            nvar,[],[],[],[],lb,ub,@(x)o_all_NonLinear_Constraints(x,conf,problem),...
            integervars,options);
        end
    elseif strcmp(conf.genStructure, 'allAntennas')
        options = optimoptions(options,'CrossoverFcn',...
            @(p,o,n,z1,z2,t)o_crossoverArrayGA(p,o,n,z1,z2,t,problem,conf));
        options = optimoptions(options,'MutationFcn',...
            @(p,o,n,z1,z2,z3,t)o_mutationArrayGA(p,o,n,z1,z2,z3,t,problem,conf));
        options = optimoptions(options,'CreationFcn',...
            @(g,f,o)o_creationArrayGA(g,f,o,problem,conf));
        if conf.multiPath
            [x,fval,exitflag,output,population,score] = ...
            ga(@(x)o_Position_Objective_optim_cost_multipath(x,conf,problem),...
            nvar,[],[],[],[],[],[],[],[],options);
        else
            [x,fval,exitflag,output,population,score] = ...
            ga(@(x)o_Position_Objective_optim_cost_singlepath(x,conf,problem),...
            nvar,[],[],[],[],[],[],[],[],options);
        end
    else
        %integervars = 1:problem.Nmax;
        % Configure options for GA
        options = optimoptions(options,'CrossoverFcn',...
            @(p,o,n,z1,z2,t)o_crossoverArrayGA(p,o,n,z1,z2,t,problem,conf));
        options = optimoptions(options,'MutationFcn',...
            @(p,o,n,z1,z2,z3,t)o_mutationArrayGA(p,o,n,z1,z2,z3,t,problem,conf));
        options = optimoptions(options,'CreationFcn',...
            @(g,f,o)o_creationArrayGA(g,f,o,problem,conf));
        if isfield(problem,'initialW') && problem.IDUserAssigned==1
            [initPopulation,initScores] = assignInitialScore(problem,conf);
            options.InitialPopulation = initPopulation;
            options.InitialScores = initScores;
        end
        if conf.multiPath
            [x,fval,exitflag,output,population,score] = ...
            ga(@(x)o_Position_Objective_optim_cost_multipath(x,conf,problem),...
            nvar,[],[],[],[],[],[],[],... % @lb,ub,(x)o_all_NonLinear_Constraints(x,conf,problem)
            [],options); %integervars,options);
        else
            [x,fval,exitflag,output,population,score] = ...
            ga(@(x)o_Position_Objective_optim_cost_singlepath(x,conf,problem),...
            nvar,[],[],[],[],[],[],[],... % lb,ub,@(x)o_all_NonLinear_Constraints(x,conf,problem),...
            [],options); %integervars,options);
        end
    end
    bestScores = bestS;
end

function [InitialPopulation,InitialScores] = assignInitialScore(problem,conf)
    % Antenna Selection
    antAlloc   = find(problem.initialW(problem.IDUserAssigned,:)~=0);
    partition  = o_antennas_to_subarrays(antAlloc,problem.Partition);
    % Amplitude Selection
    amplitudes = abs(problem.initialW(problem.IDUserAssigned,antAlloc));
    % Phase Selection
    phases     = angle(problem.initialW(problem.IDUserAssigned,antAlloc));    
    % Baseband Fbb
    Fbb = 1;  % Already accounted for in the other terms
    % Final initial population
    InitialPopulation = [partition,amplitudes,phases,Fbb];
    % Final initial score
    InitialScore = o_Position_Objective_optim_cost_singlepath(InitialPopulation, conf, problem);
    % Post-processing of the results
    InitialPopulation = repmat(InitialPopulation,conf.PopulationSize_Data,1);
    InitialScores = repmat(InitialScore,conf.PopulationSize_Data,1);
    
%     % Alternative method to compare that the scores check (to be removed)
%     mygene = InitialPopulation(1,1:end-1);
%     [handle_Conf_Array_USER,~,~,~] = o_geneToAssignment(mygene,problem,conf);
%     % Extract Rx Power (in dB)
%     id = problem.IDUserAssigned;
%     DirOK = -Inf(problem.nUsers,1);
%     DirNOK = -Inf(problem.nUsers,problem.nUsers);
%     DirOK(id) = patternAzimuth(handle_Conf_Array_USER,problem.freq,problem.thetaUsers(id),'Azimuth',problem.phiUsers(id),'Type','powerdb');
%     % Extract interference generated to others (in dB)
%     for id1 = 1:1:problem.nUsers
%         if id1~=id
%             DirNOK(id,id1) = patternAzimuth(handle_Conf_Array_USER,problem.freq,problem.thetaUsers(id1),'Azimuth',problem.phiUsers(id1),'Type','powerdb');
%         end
%     end
%     % Compute basic parameters for SINR and Capacity computations
%     chLoss_lin = ((problem.lambda ./ (4*pi*problem.dUsers(1:problem.nUsers))).^2).';  % Losses
%     Noise_lin = db2pow(problem.Noise);  % Noise power
%     Noise_lin = repmat(Noise_lin,problem.nUsers,1);
%     % Parse results for specific case
%     DirOK_lin = db2pow(DirOK);
%     DirNOK_lin = db2pow(DirNOK);
%     DirNOK_pcvd_lin = sum(DirNOK_lin,2); % Perceived interference
%     % Compute SINR and Capacities - LCMV
%     SINR_PB_lin = (DirOK_lin.*chLoss_lin) ./(DirNOK_pcvd_lin.*chLoss_lin + Noise_lin);  % Compute SINR Pass-Band (BB)
%     SINR_PB = pow2db(SINR_PB_lin);  % Compute SINR Pass-Band (BB)
%     Cap_lin = log2(1 + SINR_PB_lin);  % Compute final Capacity (bits/Hz/s)
%     if conf.verbosity >= 1
%         fprintf('* Capacity(%d): %.2f (bits/Hz/s)\n',id,Cap_lin(id));
%         fprintf('* SINR(%d): %.2f (dB)\n',id,SINR_PB(id));
%         fprintf('* Directivity IDmax: %.2f (dB)\n',DirOK(id));
%         for id1 = 1:1:problem.nUsers
%             if id1~=id
%                 fprintf('  Directivity IDmin(%d): %.2f (dB)\n',id1,DirNOK(id,id1));
%             end
%         end
%     end
%     % Retrieve initial score for GA
%     scorePrx = o_transferScore(DirOK(id),1);
%     set = (1:1:problem.nUsers);
%     set(id) = [];
%     scoreInt = o_transferScore(pow2db(sum(db2pow(DirNOK(id,set)))),0);
%     conf.Fweights(1) = 1/problem.nUsers;
%     conf.Fweights(2) = 1 - 1/problem.nUsers;
%     initialScore = (-1)*( conf.Fweights(1)*scorePrx + conf.Fweights(2)*scoreInt );
%     % Store initial scores in options variable
%     initialScores = repmat(initialScore,conf.PopulationSize_Data,1);
end