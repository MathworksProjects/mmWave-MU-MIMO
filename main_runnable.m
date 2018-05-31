%% Setup the environment
clear; clc; close all;
addpath('utilities','-end');  % Add utilities folder at the end of search path
%% Define several experiments here and override variable values accordingly
experimentList = 21;
%% Experiment selection

%% EXPERIMENT 1
if any(experimentList(:)==1);    experiment1();   end

%% EXPERIMENT 2
if any(experimentList(:)==2)
    arrRestctList    = {'None','Localized'};
    % Input parameters
    input.nIter               = 2;  % Total number of iterations
    input.nUsers              = 2;  % Number of users deployed
    input.nAntennasList       = [4 8 12 16].^2;  % Number of antennas in array
    input.Maxgenerations_Data = 100;
    input.algorithm           = 'GA';  % Heuristic algorithm
    input.detLocation         = true;  % Deterministic locations if true
    input.useCasesLocation    = true;  % Use-Case locations if true
    input.useCaseLocation     = 1;  % Use-case ID
    fileNameList = cell(length(arrRestctList),1);
    score_tot = zeros(length(arrRestctList),length(input.nAntennasList),input.Maxgenerations_Data);
    for restIdx = 1:length(arrRestctList)
        % Input parameters extra
        input.arrayRestriction = arrRestctList{restIdx};
        % Main
        fileNameList{restIdx} = experiment2(input);
        % Parse results
        load(fileNameList{restIdx},'bestScores');
        score_tot(restIdx,:,:) = bestScores;
    end
    % Save results
    fileName = strcat('temp/exp2_',input.algorithm,'_TOT_',mat2str(input.nUsers),'_',mat2str(input.detLocation),'_',mat2str(input.useCasesLocation),'_',mat2str(input.useCaseLocation));
    nAntennasList = input.nAntennasList; Maxgenerations_Data=input.Maxgenerations_Data;
    save(fileName,'score_tot','nAntennasList','arrRestctList','Maxgenerations_Data');
    % Plot results
    experiment2_plot(fileName);
end

%% EXPERIMENT 2 - PLOTTING
if any(experimentList(:)==21)
    fileName = strcat('temp/exp2_GA_TOT_2_true_true_1');
    experiment2_plot(fileName);
end

%% EXPERIMENT 3
if any(experimentList(:)==3);    experiment3();   end

%% EXPERIMENT 4
if any(experimentList(:)==4)
    nUsers = 2;
    nAntennasList = [4 6 8 10];
    nIter = 10;
    plotFLAG = true;
    experiment4(nIter,nUsers,nAntennasList,plotFLAG);
end

%% EXPERIMENT 5
if any(experimentList(:)==5)
    arrRestctList    = {'None','Localized'};
    % Output parameters
    fileNameList = cell(length(arrRestctList),1);
    SINR_PB_tot = zeros(length(input.nAntennasList),length(arrRestctList));
    SINR_BB_tot = zeros(length(input.nAntennasList),length(arrRestctList));
    Cap_tot = zeros(length(input.nAntennasList),length(arrRestctList));
    for restIdx = 1:length(arrRestctList)
        % Input parameters
        input.nIter            = 5;  % Total number of iterations
        input.nUsers           = 2;  % Number of users deployed
        input.nAntennasList    = [4 8 12 16 20 24 28 32].^2;  % Number of antennas in array
        input.arrayRestriction = arrRestctList{restIdx};
        input.algorithm        = 'GA';  % Heuristic algorithm
        input.detLocation      = true;  % Deterministic locations if true
        input.useCasesLocation = true;  % Use-Case locations if true
        input.useCaseLocation  = 1;  % Use-case ID
        plotFLAG = false;  % Plotting flag
        % Main
        fileNameList{restIdx} = experiment5(input,plotFLAG);
        % Plot
        experiment5_plot(fileNameList{restIdx});
        % Parse results
        load(fileNameList{restIdx});
        SINR_PB_lin = db2pow(SINR_PB);
        SINR_BB_lin = db2pow(SINR_BB);
        Cap_lin = db2pow(Cap);
        SINR_PB_tot(:,restIdx) = pow2db(mean(SINR_PB_lin,1)).';
        SINR_BB_tot(:,restIdx) = pow2db(mean(SINR_BB_lin,1)).';
        Cap_tot(:,restIdx) = pow2db(mean(Cap_lin,1)).';
    end
    % Save results
    fileName = strcat('temp/exp5_GA_TOT_',mat2str(input.nUsers),'_',mat2str(input.detLocation),'_',mat2str(input.useCasesLocation),'_',mat2str(input.useCaseLocation));
    nUsers = input.nUsers;  nAntennasList = input.nAntennasList;
    save(fileName,'Cap_tot','SINR_BB_tot','SINR_PB_tot','nUsers','nAntennasList','arrRestctList');
end

%% EXPERIMENT 5 - PLOTTING
if any(experimentList(:)==51)
    arrRestctList = {'None','Localized'};
    fileNameList = {'temp/exp5_GA_None_2_true_true_1','temp/exp5_GA_Localized_2_true_true_1'};
%     fileNameList = {'temp/exp5-results_GA_None_2_true_true_4','temp/exp5-results_GA_Localized_2_true_true_4'};
%     fileNameList = {'temp/exp5-results_GA_None_2_true_true_5','temp/exp5-results_GA_Localized_2_true_true_5'};
%     fileNameList = {'temp/exp5-results_GA_None_2_true_true_6','temp/exp5-results_GA_Localized_2_true_true_6'};
    for restIdx = 1:length(arrRestctList)
        experiment5_plot(fileNameList{restIdx});
    end
    fileName = strcat('temp/exp5_GA_TOT_2_true_true_1');
    load(fileName,'Cap_tot','SINR_BB_tot','SINR_PB_tot','nUsers','nAntennasList','arrRestctList');
    figure;
    plot(nAntennasList,SINR_PB_tot,'LineWidth',2,'Marker','s');
    title('Average SINR-PB','FontSize',12);
    xlabel('Number of antennas','FontSize',12);
    ylabel('SINR in dB','FontSize',12);
    legend(arrRestctList);
    grid minor;
    figure;
    plot(nAntennasList,SINR_BB_tot,'LineWidth',2,'Marker','s');
    title('Average SINR-BB','FontSize',12);
    xlabel('Number of antennas','FontSize',12);
    ylabel('SINR in dB','FontSize',12);
    legend(arrRestctList);
    grid minor;
    figure;
    plot(nAntennasList,Cap_tot,'LineWidth',2,'Marker','s');
    title('Average capacity','FontSize',12);
    xlabel('Number of antennas','FontSize',12);
    ylabel('Capacity in bits/Hz/s','FontSize',12);
    legend(arrRestctList);
    grid minor;
end

%% EXPERIMENT 7
if any(experimentList(:)==7)
    nUsersList = [2];
%     nAntennasList = [4 5 6 7 8 9 10].^2;
    nAntennasList = [2 3].^2;
    nIter = 2;
    plotFLAG = true;
    experiment7(nIter,nUsersList,nAntennasList,plotFLAG);
end

%% EXPERIMENT 8
if any(experimentList(:)==8)
    nUsersList = [2];
%     nAntennasList = [4 5 6 7 8 9 10].^2;
    nAntennasList = [2 3 4 5].^2;
    nIter = 2;
    plotFLAG = true;
    experiment8(nIter,nUsersList,nAntennasList,plotFLAG);
end

%% --------------------- EXPERIMENT IMPLEMENTATION --------------------- %%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%






%% EXPERIMENT 1
function experiment1(varargin)
    % EXPERIMENT 1 - Capacity offered
    % In this experiment we evaluate the capacity that the heuristics are able
    % to offer to the devices. Heuristics assigns antennas as a function of the
    % priority. The traffic is overloaded. The users location and channel
    % varies across simulations.
    %
    %------------- BEGIN CODE EXPERIMENT 1 --------------
    %
    fprintf('Running experiment 1...\n');
    % Load basic configuration - static and/or default
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override parameters
    problem.iat = 60;
    problem.deadline = 50;
    problem.payload = 1500*8*5e3;
    % Configure the simulation environment
    [problem,~,flows] = f_configuration(conf,problem);
    baseFlows = flows;  % For printing purposes at the end of execution
    % Main function
    [flows,CapTot,TXbitsTot,THTot,lastSlotSimm,lastSelFlow] = main(conf,problem,flows);
    % Report of single execution
    [ratioOK,ratioNOK] = f_generateReport(flows,DEBUG);
    % Plotting of single execution
    % main_plotting(problem,TXbitsTot,THTot,baseFlows,lastSelFlow);
end







%% EXPERIMENT 2 - Convergency analysis
function fileName = experiment2(input)
    fprintf('Running experiment 2...\n');
    % Store input struct in local
    nUsers              = input.nUsers;
    nIter               = input.nIter;
    nAntennasList       = input.nAntennasList;
    arrayRestriction    = input.arrayRestriction;
    Maxgenerations_Data = input.Maxgenerations_Data;
    algorithm           = input.algorithm;
    detLocation         = input.detLocation;
    useCasesLocation    = input.useCasesLocation;
    useCaseLocation     = input.useCaseLocation;
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override (problem) parameters
    problem.nUsers = nUsers;  % Number of users in the simulation
    problem.MinObjFIsSNR = true;  % (arbitrary)
    problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
    problem.arrayRestriction = 'None';  % Possibilities: "None", "Localized", "Interleaved", "DiagInterleaved"
    % Override (conf) parameters
    conf.verbosity = 2;
    conf.algorithm = algorithm;  % Heuristic algorithm
    conf.arrayRestriction = arrayRestriction;
    conf.PopulationSize_Data = 30;
    conf.Maxgenerations_Data = Maxgenerations_Data;
    conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);
    conf.MaxStallgenerations_Data = conf.Maxgenerations_Data;  % Force it to cover all the generations
    conf.FunctionTolerance_Data = 1e-10;  % Heuristics stops when not improving solution by this much
    conf.NumPhaseShifterBits = 60;  % Number of 
    conf.multiPath = false;  % LoS channel (for now)
    conf.detLocation = detLocation;  % Use fixed pre-stored locations
    conf.useCasesLocation = useCasesLocation;  % Use the use-case locations
    conf.useCaseLocation = useCaseLocation;  % Specify the use-case location
    % Output parameters
    bestScores1 = zeros(length(nAntennasList),nIter,conf.Maxgenerations_Data);  % temp
    bestScores = zeros(length(nAntennasList),conf.Maxgenerations_Data);  % final
    fileName = strcat('temp/exp2_',conf.algorithm,'_',problem.arrayRestriction,'_',mat2str(nUsers),'_',mat2str(conf.detLocation),'_',mat2str(conf.useCasesLocation),'_',mat2str(conf.useCaseLocation));
    % For each case we execute ES and the GA
    for idxAnt = 1:length(nAntennasList)
        for idxIter = 1:nIter
            fprintf('Iteration %d with nAntenas %d\n',idxIter,nAntennasList(idxAnt));
            % Configure the simulation environment. Need to place users in new
            % locations (if not fixed) and create new channels 
            % to have statistically meaningful results (if not LoS)
            [problem,~,~] = f_configuration(conf,problem);
            % Select number of antennas
            problem.N_Antennas = nAntennasList(idxAnt);
            % Adjust parameters
            problem.NxPatch = floor(sqrt(problem.N_Antennas));
            problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
            problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
            % Call heuristics
            fprintf('\t** %d Antennas and %d Users...\n',problem.N_Antennas,problem.nUsers);
            % We will paralelize the solution computations: we need (if not already
            % created) a parallelization processes pool
            gcp;
            %Create subarray partition
            problem = o_create_subarray_partition(problem);
            problem.NzPatch = problem.NxPatch;
            problem.dz = problem.dx;
            %Create the antenna handler and the data structure with all possible pos.
            problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
                'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque sí
            handle_ConformalArray = phased.URA([problem.NyPatch,problem.NzPatch],...
                'Lattice','Rectangular','Element',problem.handle_Ant,...
            'ElementSpacing',[problem.dy,problem.dz]);
            problem.possible_locations = handle_ConformalArray.getElementPosition;
            % Boolean flag indicating if we have already found a feasible solution
            problem = o_compute_antennas_per_user(problem,1:nUsers);
            % We will accumulate in the assignments_status var the
            % antennas / subarrays assigned as soon as we assign them
            [~,orderedIndices] = sort(problem.MinObjF,'descend');
            u = orderedIndices(1);
            problem.IDUserAssigned = u;
            % Genetic Algorithm
            fprintf('Solving... (Genetic Algorithm)\n')
            [~,~,~,~,tempBS] = ...
                o_solveSingleNmaxUserInstance(conf,problem,...
                problem.NmaxArray(problem.IDUserAssigned));
            bestScores1(idxAnt,idxIter,1:length(tempBS)) = tempBS;
            % Smoothen out results in tail
            if length(tempBS)<conf.Maxgenerations_Data
                stillCover = length(tempBS)+1:conf.Maxgenerations_Data;
                bestScores1(idxAnt,idxIter,stillCover) = repmat(bestScores1(idxAnt,idxIter,length(tempBS)),1,length(stillCover));
            end
        end
        t = mean(bestScores1(idxAnt,:,:),2);
        bestScores(idxAnt,:) = t(:).';
        save(fileName,'bestScores',...
             'nUsers','nAntennasList','arrayRestriction',...
             'detLocation','useCaseLocation','useCaseLocation');
    end
end

%% EXPERIMENT 2 - Plotting
function experiment2_plot(fileName)
    % Plot results
    load(fileName,'score_tot','nAntennasList','arrRestctList');
    Maxgenerations_Data = size(score_tot,3);  %#ok
    figure; hold on;
    colList = {'r','b','g'};
    colList{1} = [0 0 255]./255;
    colList{2} = [0 51 102]./255;
    colList{3} = [0 128 255]./255;
    colList{4} = [0 180 255]./255;
    colList{5} = [0 220 255]./255;
    lineStyleList = {'-','-.','--',':'};
    leg = cell(length(arrRestctList)*length(nAntennasList),1);  %#ok
    for restIdx = 1:length(arrRestctList)
        col = colList{restIdx};
        for antIdx = 1:length(nAntennasList)
%             lin = lineStyleList{1+mod(antIdx-1,length(nAntennasList))};
            res = score_tot(restIdx,antIdx,:);  %#ok
            score_TOT(:,(restIdx-1)*length(nAntennasList) + antIdx) = res(:);  %#ok
            leg((restIdx-1)*length(nAntennasList) + antIdx) = strcat(arrRestctList{restIdx},{' '},'--',{' '},mat2str(nAntennasList(antIdx)));  %#ok
        end
    end
    plot(1:Maxgenerations_Data,score_TOT,'LineWidth',1,'color',col,'LineStyle',lin,'Marker','.','MarkerSize',10);
    grid minor;
    xlabel('Generations','FontSize',12);
    ylabel('Fitness value','FontSize',12);
    title('Convergency analysis','FontSize',12);
    hleg = legend(leg);
    set(hleg,'FontSize',12);
end







%% EXPERIMENT 3 - Performance comparisson against null-forcing technique
% To-Do. The comparisson needs to be agains a more updated technique such
% as the JSDM. Waiting for reply from Kaushik to write to them and get the
% code.







%% EXPERIMENT 4 - Convergency analysis
function experiment4(nIter,nUsers,nAntennasList,plotFLAG)
    % EXPERIMENT 4 -- 
    % 
    % Aim: Obtain the approximate number of generations (operational cycles)
    % that we need to obtain a certain quality of solution by comparing the
    % Genetic Algorithm convergence with the global optimum found by means
    % of an Exhaustive Search
    % 
	% Assumptions (Fixed):
    %   1. User location: Fixed, from config file.
    %   2. Sub-array geometry: 'None'.
    %   3. Antenna Array geometry: Fixed to URA.
    %   4. Algorithm: GA & ES
    % Variable:
    %   1. Antenna Array size variable: nAntennasList
    %   2. Population size: Fixed with all the rest of GA parameters, in
    %   order to have a unique dependence on #generations
    % 
    % Syntax:  [] =
    % experiment4(nIters,nUsers,nAntennasList,plotFLAG)
    % 
    % Inputs:
    %    nIter - Number of iterations to extract average values
    %    nUsers - Number of users considered
    %    nAntennaList - Number of antenas (set of)
    %    plotFLAG - True for plotting directivity and antenna allocation
    %
    % Outputs: None
    %
    %------------- BEGIN CODE EXPERIMENT 4 --------------
    %
    fprintf('Running experiment 4...\n');
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override (problem) parameters
    % Override (problem) parameters
    problem.nUsers = nUsers;  % Number of users in the simulation
    problem.MinObjFIsSNR = true;  % (arbitrary)
    problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
    problem.arrayRestriction = 'None';  % Possibilities: "None", "Localized", "Interleaved", "DiagInterleaved"
    % Override (conf) parameters
    conf.verbosity = 0;
    conf.NumPhaseShifterBits = 2;  % Number of 
    conf.NbitsAmplitude = 2;
    conf.FunctionTolerance_Data = 1e-6;  % Heuristics stops when not improving solution by this much
    conf.multiPath = false;  % LoS channel (for now)
    
    % Override GA parameters
    conf.PopulationSize_Data = 150;
    conf.Maxgenerations_Data = 100;
    conf.EliteCount_Data = 25;
    conf.MaxStallgenerations_Data = 40;  % Force it to cover all the generations
    %h1 = figure;
    %hold on
    
    globalMin = zeros(1,length(nAntennasList));
    globalMin_Loc = zeros(1,length(nAntennasList));
    bestScores = zeros(length(nAntennasList),nIter,conf.Maxgenerations_Data);
    bestScores_Loc = zeros(length(nAntennasList),nIter,conf.Maxgenerations_Data);
    
    % For each case we execute ES and the GA
    for idxAnt = 1:length(nAntennasList)
        for idxIter = 1:nIter
            fprintf('Iteration %d with nAntenas %d\n',idxIter,nAntennasList(idxAnt));
            % Configure the simulation environment. Need to place users in new
            % locations (if not fixed) and create new channels 
            % to have statistically meaningful results (if not LoS)
            [problem,~,~] = f_configuration(conf,problem);
            % Select number of antennas
            problem.N_Antennas = nAntennasList(idxAnt);
            % Adjust parameters
            problem.NxPatch = floor(sqrt(problem.N_Antennas));
            problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
            problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
            % Call heuristics
            fprintf('\t** %d Antennas and %d Users...\n',problem.N_Antennas,problem.nUsers);
            % We will paralelize the solution computations: we need (if not already
            % created) a parallelization processes pool
            gcp;

            % Create subarray partition
            problem = o_create_subarray_partition(problem);

            problem.NzPatch = problem.NxPatch;
            problem.dz = problem.dx;

            % Create the antenna handler and the data structure with all possible pos.
            problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
                'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque sí
            handle_ConformalArray = phased.URA([problem.NyPatch,problem.NzPatch],...
                'Lattice','Rectangular','Element',problem.handle_Ant,...
            'ElementSpacing',[problem.dy,problem.dz]);

            problem.possible_locations = handle_ConformalArray.getElementPosition;

            % Boolean flag indicating if we have already found a feasible solution
            problem = o_compute_antennas_per_user(problem,1:nUsers);
            % We will accumulate in the assignments_status var the
            % antennas / subarrays assigned as soon as we assign them
            [~,orderedIndices] = sort(problem.MinObjF,'descend');
            u = orderedIndices(1);
            problem.IDUserAssigned = u;
            
            % First we execute the Exhaustive Search
            % Only the first iteration, because the result will be the same
            % in other iterations (if assuming static environment)
            if idxIter == 1
                conf.arrayRestriction = 'None';
                conf.algorithm = 'ES';  % Heuristic algorithm
                fprintf('Solving... (Exhaustive Search - Free subarray allocation)\n')
                [~,~,~,~,globalMin(idxAnt)] = ...
                    o_solveSingleNmaxUserInstance(conf,problem,...
                    problem.NmaxArray(problem.IDUserAssigned));
                fprintf('Solved!\n')
%                 figure(h1)
%                 disp(globalMin)
%                 line(1:conf.Maxgenerations_Data,ones(1,conf.Maxgenerations_Data)*globalMin);
%                 drawnow
                conf.arrayRestriction = 'Localized';  % Heuristic algorithm
                fprintf('Solving... (Exhaustive Search - Localized)\n')
                [~,~,~,~,globalMin_Loc(idxAnt)] = ...
                    o_solveSingleNmaxUserInstance(conf,problem,...
                    problem.NmaxArray(problem.IDUserAssigned));
                fprintf('Solved!\n')
            end
            
            % And secondly using Genetic Algorithm
            conf.arrayRestriction = 'None';
            conf.algorithm = 'GA';  % Heuristic algorithm
            fprintf('Solving... (Genetic Algorithm)\n')
            [~,~,~,~,tempBS] = ...
                o_solveSingleNmaxUserInstance(conf,problem,...
                problem.NmaxArray(problem.IDUserAssigned));
            bestScores(idxAnt,idxIter,1:length(tempBS)) = tempBS;
            fprintf('Solved!\n')
%             figure(h1)
%             disp(bestScores)
%             line(1:length(bestScores),bestScores);
%             drawnow
            % And secondly using Genetic Algorithm
            conf.arrayRestriction = 'Localized';
            conf.algorithm = 'GA';  % Heuristic algorithm
            fprintf('Solving... (Genetic Algorithm)\n')
            [~,~,~,~,tempBS] = ...
                o_solveSingleNmaxUserInstance(conf,problem,...
                problem.NmaxArray(problem.IDUserAssigned));
            bestScores_Loc(idxAnt,idxIter,1:length(tempBS)) = tempBS;
            fprintf('Solved!\n')
            save('temp/exp4-results_so_far','globalMin','globalMin_Loc',...
                'bestScores','bestScores_Loc','nUsers','nAntennasList');
        end
    end
end







%% EXPERIMENT 5
function fileName = experiment5(input,plotFLAG)
    % EXPERIMENT 5 -- 
    % 
    % Aim: Evaluate the average received power (Prx) at the intended users
    % and Interference (Int) generated at other users. The secuential
    % allocation policy leads to an unfair allocation policy, which leads
    % to different Prx and Int values. This experiment analyzes how the
    % size of the antenna array and the place in the priority list impact
    % on the performance of the system.
    % 
    % The experiment stores the relevant variables in a .mat file with the
    % following format:     temp/exp5_A_B_C_D_E_F.mat, 
    % where 'A' is the heuristic algorithm used, 'B' is the array geometry,
    % 'C' is the num. of users, 'D' is the deterministic localization, 'E'
    % is the use case localization ('E' overrides 'D') and 'F' is the use
    % case localization used.
    % 
    % The .mat file contains the info (dim arrays [nUsers x nAntennasList]):
    %    Cap - Capacity in b/Hz/s 
    %    SINR_BB - BaseBand SINR (dB)
    %    SINR_PB - PassBand SINR (dB)
    %    DirOK - Directivity to the intended transmitter (dB)
    %    DirNOK_gntd - Directivity generated to other nodes (dB)
    %    DirNOK_pcvd - Directivity perceived due to interfeering nodes (dB)
    %    nUsers - The number of users inputed (integer)
    %    nAntennasList - The antenna size inputted (array)
    %    arrayRestriction - The array restriction inputted (cell)
    %    detLocation - Use fixed pre-stored locations (boolean)
    %    useCasesLocation - Use the use-case locations (boolean)
    %    useCaseLocation - Specify the use-case location (integer)
    % 
	% Assumptions (Fixed):
    %   1. Number of antennas: Same across users and prop. to Array size.
    %   2. Number of users: nUsers.
    %   3. User location: From config file.
    %   4. Sub-array geometry: Fixed to 'None', 'Localized', etc.
    %   5. Antenna Array geometry: Fixed to URA.
    %   6. Algorithm: Fixed to 'GA', 'PSO' or 'PS'.
    % Variable:
    %   1. Antenna Array size variable: nAntennasList
    % 
    % Syntax:  [CapTot,SINRTot,DirOKTot,DirOKAv,DirNOKTot,DirNOKAv] =
    % experiment5(nIter,nUsers,nAntennasList,plotFLAG)
    % 
    % Inputs (stored in struct 'input'):
    %    nIter - Number of iterations to extract average values
    %    nUsers - Number of users considered
    %    nAntennaList - Number of antenas
    %    arrayRestriction - Defines the sub-array restriction at the BS
    %    algorithm - Defines tre heuristic algorithm to run
    %    detLocation - Determines whether to use static locations
    %    useCasesLocation - Determines wheter to use the loc. use-cases
    %    useCaseLocation  - Determines the use-case ID
    %    plotFLAG - True for plotting directivity and antenna allocation
    %
    % Outputs: 
    %    fileName - The name of the .mat file with the results
    %
    %------------- BEGIN CODE EXPERIMENT 5 --------------
    %
    fprintf('Running experiment 5...\n');
    % Store input struct in local
    nUsers           = input.nUsers;
    nIter            = input.nIter;
    nAntennasList    = input.nAntennasList;
    arrayRestriction = input.arrayRestriction;
    algorithm        = input.algorithm;
    detLocation      = input.detLocation;
    useCasesLocation = input.useCasesLocation;
    useCaseLocation  = input.useCaseLocation;
    fprintf('Input parameters:\n\tnUsers:\t%d\n\tIter:\t%d\n\tAntennasList:\t%s\n\tarrayRestriction:\t%s\n\talgorithm:\t%s\n\tdetLocation:\t%s\n\tuseCasesLocation:\t%s\n\tuseCaseLocation:\t%d\n',nUsers,nIter,mat2str(nAntennasList),arrayRestriction,algorithm,mat2str(detLocation),mat2str(useCasesLocation),useCaseLocation);
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override (problem) parameters
    problem.nUsers = nUsers;  % Number of users in the simulation
    problem.MinObjFIsSNR = true;  % (arbitrary)
    problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
    problem.arrayRestriction = arrayRestriction;  % Possibilities: "None", "Localized", "Interleaved", "DiagInterleaved"
    % Override (conf) parameters
    conf.verbosity = 1;
    conf.algorithm = algorithm;  % Heuristic algorithm
    conf.PopulationSize_Data = 30;
    conf.Maxgenerations_Data = 30;
    conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);
    conf.MaxStallgenerations_Data = ceil(conf.Maxgenerations_Data/10);  % Force it to cover all the generations
    conf.FunctionTolerance_Data = 1e-10;  % Heuristics stops when not improving solution by this much
    conf.NumPhaseShifterBits = 60;  % Number of 
    conf.multiPath = false;  % LoS channel (for now)
    conf.detLocation = detLocation;  % Use fixed pre-stored locations
    conf.useCasesLocation = useCasesLocation;  % Use the use-case locations
    conf.useCaseLocation = useCaseLocation;  % Specify the use-case location
	% Configure basic parameters
    candSet = (1:1:problem.nUsers);  % Set of users to be considered
	% Create output variables
    CapTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    SINRTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    DirOKTot = -Inf(problem.nUsers,length(nAntennasList),nIter);
    DirNOKTot = -Inf(problem.nUsers,problem.nUsers,length(nAntennasList),nIter);
    fileName_temp = strcat('temp/exp5_',problem.arrayRestriction,'_',mat2str(nUsers),'_',conf.algorithm,'_so_far');
    fileName = strcat('temp/exp5_',conf.algorithm,'_',problem.arrayRestriction,'_',mat2str(nUsers),'_',mat2str(conf.detLocation),'_',mat2str(conf.useCasesLocation),'_',mat2str(conf.useCaseLocation));
    % Main execution
    for idxAnt = 1:length(nAntennasList)
        for idxIter = 1:nIter
            fprintf('Iteration %d\n',idxIter);
            % Configure the simulation environment. Need to place users in new
            % locations and create new channels to have statistically
            % meaningful results
            [problem,~,~] = f_configuration(conf,problem);
            % Select number of antennas
            problem.N_Antennas = nAntennasList(idxAnt);
            % Adjust parameters
            problem.NxPatch = floor(sqrt(problem.N_Antennas));
            problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
            problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
            % Call heuristics
            fprintf('\t** %d Antennas and %d Users...\n',problem.N_Antennas,problem.nUsers);
            [~,W,~,estObj] = f_heuristics(problem,conf,candSet);
            % Heuristics - Post Processing
            if conf.MinObjFIsSNR;     CapTot(:,idxAnt,idxIter)  = log2(estObj+1);  % in bps/Hz
                                      SINRTot(:,idxAnt,idxIter) = pow2db(estObj);  % in dB
            else;                     CapTot(:,idxAnt,idxIter)  = estObj;  % in bps/Hz
                                      SINRTot(:,idxAnt,idxIter) = pow2db(2.^(estTH/problem.Bw) - 1);  % in dB
            end
            % Reconstruct array
            % Create handle per user
            problem1 = o_create_subarray_partition(problem);
            problem1.NzPatch = problem1.NxPatch;
            problem1.dz = problem1.dx;
            problem1.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                                    [problem1.freq-(problem1.Bw/2) problem1.freq+(problem1.Bw/2)],...
                                    'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque sí
            handle_ConformalArray = phased.URA([problem1.NyPatch,problem1.NzPatch],...
                                    'Lattice','Rectangular','Element',problem1.handle_Ant,...
                                    'ElementSpacing',[problem1.dy,problem1.dz]);
            problem1.possible_locations = handle_ConformalArray.getElementPosition;
            for id = 1:1:problem1.nUsers
                problem1.ant_elem = sum(W(id,:)~=0);
                relevant_positions = (W(id,:)~=0);
                Taper_user = W(id,relevant_positions);
                handle_Conf_Array = phased.ConformalArray('Element',problem1.handle_Ant,...
                                      'ElementPosition',...
                                      [zeros(1,problem1.ant_elem);...
                                      problem1.possible_locations(2,relevant_positions);...
                                      problem1.possible_locations(3,relevant_positions)],...
                                      'Taper',Taper_user);
                % Extract Rx Power (in dB)
                DirOKTot(id,idxAnt,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id),'Azimuth',problem.phiUsers(id),'Type','powerdb');
                fprintf('* Directivity IDmax: %.2f (dB)\n',DirOKTot(id,idxAnt,idxIter));
                % Extract interference generated to others (in dB)
                for id1 = 1:1:problem1.nUsers
                    if id1~=id
                        DirNOKTot(id,id1,idxAnt,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id1),'Azimuth',problem.phiUsers(id1),'Type','powerdb');
                        fprintf('  Directivity IDmin(%d): %.2f (dB)\n',id1,DirNOKTot(id,id1,idxAnt,idxIter));
                    end
                end
                problem1.IDUserAssigned = id;
                if plotFLAG
                    % Plot beam pattern obtained with assignation and BF configuration
                    o_plotAssignment_mod(problem1, handle_Conf_Array);
                end
            end
            if plotFLAG
                % Plot assignation
                px = problem1.possible_locations(3,:);  % Antenna allocation on x-axis
                py = problem1.possible_locations(2,:);  % Antenna allocation on y-axis
                pz = problem1.possible_locations(1,:);  % Antenna allocation on z-axispatch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
                patch = o_getPatch(problem1.NxPatch,problem1.NyPatch,px,py);
                arrays = o_getArrays(problem1.nUsers,W,px,py,pz);
                o_plot_feasible_comb(problem1,conf,patch,arrays);
            end
            save(fileName_temp,'DirOKTot','DirNOKTot','nUsers','nAntennasList');
        end
    end
    % Convert back to Watts (from dB)
    DirOKTot_lin = db2pow(DirOKTot);
    DirNOKTot_lin = db2pow(DirNOKTot);
    % Compute average Directivities
    DirOK_lin = zeros(nUsers,length(nAntennasList));  % Directivity generated by intended user
    DirNOK_gntd_lin = zeros(nUsers,length(nAntennasList));  % Generated interference by intended user
    DirNOK_pcvd_lin = zeros(nUsers,length(nAntennasList));  % Perceived interference by intended user
    for antIdx = 1:length(nAntennasList)
        DirOK_lin(:,antIdx) = mean(DirOKTot_lin(:,antIdx,:),3);
        DirNOK_gntd_lin(:,antIdx) = sum(mean(DirNOKTot_lin(:,:,antIdx,:),4),2); % Generated interference 
        DirNOK_pcvd_lin(:,antIdx) = sum(mean(DirNOKTot_lin(:,:,antIdx,:),4),1); % Perceived interference
    end
    DirOK = pow2db(DirOK_lin);  %#ok % Directivity generated to intended user
    DirNOK_gntd = pow2db(DirNOK_gntd_lin);  %#ok % Directivity being generated by intended user
    DirNOK_pcvd = pow2db(DirNOK_pcvd_lin);  %#ok % Directivity inflicted to intended user
    % Compute SINR and Capacities
    Ptx_lin = db2pow(problem.Ptx);  % Initial transmit power
    Ptx_lin = repmat(Ptx_lin,1,length(nAntennasList));
    chLoss_lin = (((4*pi*problem.dUsers(1:nUsers)) ./ problem.lambda).^2 ).';  % Losses
    chLoss_lin = repmat(chLoss_lin,1,length(nAntennasList));
    Noise_lin = db2pow(problem.Noise);  % Noise power
    Noise_lin = repmat(Noise_lin,1,length(nAntennasList));
    SINR_PB_lin = (Ptx_lin.*DirOK_lin.*chLoss_lin) ./ (Ptx_lin.*DirNOK_pcvd_lin.*chLoss_lin + Noise_lin);  % SINR
    SINR_PB = pow2db(SINR_PB_lin);  %#ok
    SINR_BB_lin = mean(db2pow(SINRTot),3);  % Compute SINR Base-Band (BB)
    SINR_BB = pow2db(SINR_BB_lin);  %#ok
    Cap = mean(CapTot,3);  %#ok % Compute Average Capacities in the system
    % Store results in .mat file
    save(fileName,'Cap','SINR_BB','SINR_PB',...
             'DirOK','DirNOK_gntd','DirNOK_pcvd','DirOKTot','DirNOKTot',...
             'nUsers','nAntennasList','arrayRestriction',...
             'detLocation','useCaseLocation','useCaseLocation');
end

%% EXPERIMENT 5 - PLOTTING
function experiment5_plot(fileName)
    % EXPERIMENT 5 - Plotting results
    % Load results
    load(fileName);
    % Get figure number
    h = findobj('type','figure');
    figNum = length(h) + 1;
    % Plot Directivities
    figure(figNum);  figNum = figNum + 1;
    leg = cell(nUsers,1);
    for id = 1:nUsers
        subplot(1,3,1); hold on;
        plot(nAntennasList,DirOK(id,:),'LineWidth',2,'Marker','s');
        subplot(1,3,2); hold on;
        plot(nAntennasList,DirNOK_gntd(id,:),'LineWidth',2,'Marker','s');
        subplot(1,3,3); hold on;
        plot(nAntennasList,DirNOK_pcvd(id,:),'LineWidth',2,'Marker','s');
        leg{id} = cell2mat(strcat('user',{' '},num2str(id)));
    end
    subplot(1,3,1);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    title('Directivity to intended user','FontSize',12);
    legend(leg,'FontSize',12);
    subplot(1,3,2);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    title('Interference generated to other users','FontSize',12);
    legend(leg,'FontSize',12);
    subplot(1,3,3);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    title('Interference generated to intended user','FontSize',12);
    suptt = strcat(mat2str(nUsers),{' '},'usr -',{' '},arrayRestriction,{' '},'geometry');
    tit = suptitle(suptt{:});
    set(tit,'FontSize',12)
    legend(leg,'FontSize',12);
    % Plot perceived SINRs
    figure(figNum);  figNum = figNum + 1;
    for id = 1:nUsers
        subplot(1,2,1); hold on;
        plot(nAntennasList,SINR_BB(id,:),'LineWidth',2,'Marker','s');
        subplot(1,2,2); hold on;
        plot(nAntennasList,SINR_PB(id,:),'LineWidth',2,'Marker','s');
    end
    subplot(1,2,1);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    tit = strcat('BB-SINR -',{' '},mat2str(nUsers),{' '},'usr -',{' '},arrayRestriction,{' '},'geometry');
    title(tit,'FontSize',12);
    legend(leg,'FontSize',12);
    subplot(1,2,2);
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Power in dB','FontSize',12);
    tit = strcat('PB-SINR -',{' '},mat2str(nUsers),{' '},'usr -',{' '},arrayRestriction,{' '},'geometry');
    title(tit,'FontSize',12);
    legend(leg,'FontSize',12);
    % Plot perceived Capacities
    figure(figNum);  figNum = figNum + 1;
    for id = 1:nUsers
        hold on;
        plot(nAntennasList,Cap(id,:),'LineWidth',2,'Marker','s');
    end
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Capacity in bits/Hz/s','FontSize',12);
    tit = strcat('Capacities  -',{' '},mat2str(nUsers),{' '},'usr -',{' '},arrayRestriction,{' '},'geometry');
    title(tit,'FontSize',12);
    legend(leg,'FontSize',12);
    % Plot average Capacities
    figure(figNum);  figNum = figNum + 1;
    Cap_lin = db2pow(Cap);
    Cap_av = pow2db(mean(Cap_lin,1));
	plot(nAntennasList,Cap_av,'LineWidth',2,'Marker','s');
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('Average Capacity in bits/Hz/s','FontSize',12);
    tit = strcat('Capacity  -',{' '},mat2str(nUsers),{' '},'usr -',{' '},arrayRestriction,{' '},'geometry');
    title(tit,'FontSize',12);
    % Plot average SINR (BB)
    figure(figNum);  figNum = figNum + 1;
    SINR_BB_lin = db2pow(SINR_BB);
    SINR_BB_av = pow2db(mean(SINR_BB_lin,1));
	plot(nAntennasList,SINR_BB_av,'LineWidth',2,'Marker','s');
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('SINR in dB','FontSize',12);
    tit = strcat('Average BB-SINR -',{' '},mat2str(nUsers),{' '},'usr -',{' '},arrayRestriction,{' '},'geometry');
    title(tit,'FontSize',12);
    % Plot average SINR (PB)
    figure(figNum);  figNum = figNum + 1;                              %#ok
    SINR_PB_lin = db2pow(SINR_PB);
    SINR_PB_av = pow2db(mean(SINR_PB_lin,1));
	plot(nAntennasList,SINR_PB_av,'LineWidth',2,'Marker','s');
    grid minor;
    xlabel('Number of available antennas','FontSize',12);
    ylabel('SINR in dB','FontSize',12);
    tit = strcat('Average PB-SINR -',{' '},mat2str(nUsers),{' '},'usr -',{' '},arrayRestriction,{' '},'geometry');
    title(tit,'FontSize',12);
end







%% EXPERIMENT 6
function [Cap,SINR_BB,SINR_PB,DirOK,DirNOK_gntd,DirNOK_pcvd] = experiment6(nIter,nUsers,nAntennasList,plotFLAG)
    %
    fprintf('Running experiment 6...\n');
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override (problem) parameters
    problem.nUsers = nUsers;  % Number of users in the simulation
    problem.MinObjFIsSNR = true;  % (arbitrary)
    problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
    problem.arrayRestriction = 'Localized';  % Possibilities: "None", "Localized", "Interleaved", "DiagInterleaved"
    % Override (conf) parameters
    conf.verbosity = 0;
    conf.algorithm = 'GA';  % Heuristic algorithm
    conf.NumPhaseShifterBits = 60;  % Number of 
    conf.FunctionTolerance_Data = 1e-10;  % Heuristics stops when not improving solution by this much
    conf.multiPath = false;  % LoS channel (for now)
	% Configure basic parameters
    candSet = (1:1:problem.nUsers);  % Set of users to be considered
	% Create output variables
    CapTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    SINRTot = zeros(problem.nUsers,length(nAntennasList),nIter);
    DirOKTot = -Inf(problem.nUsers,length(nAntennasList),nIter);
    DirNOKTot = -Inf(problem.nUsers,problem.nUsers,length(nAntennasList),nIter);
    fileName_temp = strcat('temp/exp5-results_',problem.arrayRestriction,'_',mat2str(nUsers),'_',conf.algorithm,'_so_far');
    fileName = strcat('temp/exp5-results_',problem.arrayRestriction,'_',mat2str(nUsers),'_',conf.algorithm);
    % Linearize combinations and asign Population size (To be replaced with
    % convergency analysis values)
%     totComb = log10(problem.nUsers.*factorial(ceil(nAntennasList/problem.nUsers)));
%     maxPop = 70;  % Maximum population size
%     minPop = 40;  % Minimum population size
%     slope = (maxPop - minPop) / (totComb(end)-totComb(1));
%     ordIdx = minPop - slope*totComb(1);
%     PopSizeList = ceil(slope*totComb + ordIdx);
    PopSizeList = 150*ones(length(nAntennasList),1);
    % Main execution
    for idxAnt = 1:length(nAntennasList)
        conf.PopulationSize_Data = PopSizeList(idxAnt);
        conf.Maxgenerations_Data = 150;
        conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);
        conf.MaxStallgenerations_Data = ceil(conf.Maxgenerations_Data/10);  % Force it to cover all the generations
        for idxIter = 1:nIter
            fprintf('Iteration %d with PopSize %d\n',idxIter,PopSizeList(idxAnt));
            % Configure the simulation environment. Need to place users in new
            % locations and create new channels to have statistically
            % meaningful results
            [problem,~,~] = f_configuration(conf,problem);
            % Select number of antennas
            problem.N_Antennas = nAntennasList(idxAnt);
            % Adjust parameters
            problem.NxPatch = floor(sqrt(problem.N_Antennas));
            problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
            problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
            % Call heuristics
            fprintf('\t** %d Antennas and %d Users...\n',problem.N_Antennas,problem.nUsers);
            [~,W,~,estObj] = f_heuristics(problem,conf,candSet);
            % Heuristics - Post Processing
            if conf.MinObjFIsSNR;     CapTot(:,idxAnt,idxIter)  = log2(estObj+1);  % in bps/Hz
                                      SINRTot(:,idxAnt,idxIter) = 10*log10(estObj);  % in dB
            else;                     CapTot(:,idxAnt,idxIter)  = estObj;  % in bps/Hz
                                      SINRTot(:,idxAnt,idxIter) = 10*log10(2.^(estTH/problem.Bw) - 1);  % in dB
            end
            % Reconstruct array
            % Create handle per user
            problem1 = o_create_subarray_partition(problem);
            problem1.NzPatch = problem1.NxPatch;
            problem1.dz = problem1.dx;
            problem1.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                                    [problem1.freq-(problem1.Bw/2) problem1.freq+(problem1.Bw/2)],...
                                    'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque sí
            handle_ConformalArray = phased.URA([problem1.NyPatch,problem1.NzPatch],...
                                    'Lattice','Rectangular','Element',problem1.handle_Ant,...
                                    'ElementSpacing',[problem1.dy,problem1.dz]);
            problem1.possible_locations = handle_ConformalArray.getElementPosition;
            for id = 1:1:problem1.nUsers
                problem1.ant_elem = sum(W(id,:)~=0);
                relevant_positions = (W(id,:)~=0);
                Taper_user = W(id,relevant_positions);
                handle_Conf_Array = phased.ConformalArray('Element',problem1.handle_Ant,...
                                      'ElementPosition',...
                                      [zeros(1,problem1.ant_elem);...
                                      problem1.possible_locations(2,relevant_positions);...
                                      problem1.possible_locations(3,relevant_positions)],...
                                      'Taper',Taper_user);
                % Extract Rx Power (in dB)
                DirOKTot(id,idxAnt,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id),'Azimuth',problem.phiUsers(id),'Type','powerdb');
                fprintf('* Directivity IDmax: %.2f (dB)\n',DirOKTot(id,idxAnt,idxIter));
                % Extract interference generated to others (in dB)
                for id1 = 1:1:problem1.nUsers
                    if id1~=id
                        DirNOKTot(id,id1,idxAnt,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id1),'Azimuth',problem.phiUsers(id1),'Type','powerdb');
                        fprintf('  Directivity IDmin(%d): %.2f (dB)\n',id1,DirNOKTot(id,id1,idxAnt,idxIter));
                    end
                end
                % Compare performance with other Beamforming mechanisms
                % Simulate a test signal using a simple rectangular pulse
                t = linspace(0,0.3,300)';
                testsig = zeros(size(t));
                testsig(201:205) = 1;
                % Incident signal
                angle_of_arrival = [problem.phiUsers(1);problem.thetaUsers(1)];
                x = collectPlaneWave(handle_ConformalArray,testsig,angle_of_arrival,problem.freq);
                % Add AWGN to signal
                rng default
                npower = 0.5;
                x = x + sqrt(npower/2)*(randn(size(x)) + 1i*randn(size(x)));
                % Create interference
                jammer_angle = [problem.phiUsers(2);problem.thetaUsers(2)];
                jamsig = collectPlaneWave(handle_ConformalArray,x,jammer_angle,problem.freq);
                % Add AWGN to jamming signal
                noise = sqrt(noisePwr/2)*...
                    (randn(size(jamsig)) + 1j*randn(size(jamsig)));
                jamsig = jamsig + noise;
                rxsig = x + jamsig;
                % Create conventional Beamformer
                convbeamformer = phased.PhaseShiftBeamformer('SensorArray',handle_Conf_Array,...
                                'OperatingFrequency',problem.freq,'Direction',angle_of_arrival,...
                                'WeightsOutputPort',true);
                [~,W_convent] = convbeamformer(rxsig);
                % Create LCMV Beamformer
                steeringvector = phased.SteeringVector('SensorArray',handle_Conf_Array,...
                                 'PropagationSpeed',physconst('LightSpeed'));
                LCMVbeamformer = phased.LCMVBeamformer('DesiredResponse',1,...
                                 'TrainingInputPort',true,'WeightsOutputPort',true);
                LCMVbeamformer.Constraint = steeringvector(problem.freq,angle_of_arrival);
                LCMVbeamformer.DesiredResponse = 1;
                [~,wLCMV] = LCMVbeamformer(rxsig,jamsig);
                problem1.IDUserAssigned = id;
                if plotFLAG
                    % Plot beam pattern obtained with assignation and BF configuration
                    o_plotAssignment_mod(problem1, handle_Conf_Array);
                    % Plot beam pattern obtained with LCMV Beamforming
                    figure;
                    subplot(211)
                    pattern(handle_Conf_Array,problem.freq,(-180:180),0,'PropagationSpeed',physconst('LightSpeed'),...
                                'CoordinateSystem','rectangular','Type','powerdb','Normalize',true,...
                                'Weights',W_convent)
                    title('Array Response with Conventional Beamforming Weights');
                    subplot(212)
                    pattern(handle_Conf_Array,problem.freq,(-180:180),0,'PropagationSpeed',physconst('LightSpeed'),...)
                                'CoordinateSystem','rectangular','Type','powerdb','Normalize',true,...
                                'Weights',wLCMV)
                    title('Array Response with LCMV Beamforming Weights');
                end
            end
            if plotFLAG
                % Plot assignation
                px = problem1.possible_locations(3,:);  % Antenna allocation on x-axis
                py = problem1.possible_locations(2,:);  % Antenna allocation on y-axis
                pz = problem1.possible_locations(1,:);  % Antenna allocation on z-axispatch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
                patch = o_getPatch(problem1.NxPatch,problem1.NyPatch,px,py);
                arrays = o_getArrays(problem1.nUsers,W,px,py,pz);
                o_plot_feasible_comb(problem1,conf,patch,arrays);
            end
            save(fileName_temp,'DirOKTot','DirNOKTot','nUsers','nAntennasList');
        end
    end
    % Convert back to Watts (from dB)
    DirOKTot_lin = db2pow(DirOKTot);
    DirNOKTot_lin = db2pow(DirNOKTot);
    % Compute average Directivities
    DirOK_lin = zeros(nUsers,length(nAntennasList));  % Directivity generated by intended user
    DirNOK_gntd_lin = zeros(nUsers,length(nAntennasList));  % Generated interference by intended user
    DirNOK_pcvd_lin = zeros(nUsers,length(nAntennasList));  % Perceived interference by intended user
    for antIdx = 1:length(nAntennasList)
        DirOK_lin(:,antIdx) = mean(DirOKTot_lin(:,antIdx,:),3);
        DirNOK_gntd_lin(:,antIdx) = sum(mean(DirNOKTot_lin(:,:,antIdx,:),4),1); % Generated interference 
        DirNOK_pcvd_lin(:,antIdx) = sum(mean(DirNOKTot_lin(:,:,antIdx,:),4),2); % Perceived interference
    end
    DirOK = pow2db(DirOK_lin);  % Directivity generated to intended user
    DirNOK_gntd = pow2db(DirNOK_gntd_lin);  % Directivity being generated by intended user
    DirNOK_pcvd = pow2db(DirNOK_pcvd_lin);  % Directivity inflicted to intended user
    % Compute SINR and Capacities
    chLoss = pow2db( ((4*pi*problem.dUsers(1:nUsers)) ./ problem.lambda).^2 ).';  % Losses
    chLoss = repmat(chLoss,1,length(nAntennasList));
%     Ptx = repmat(problem.Ptx,1,length(nAntennasList));  % Initial transmit power
    Noise_lin = repmat(db2pow(problem.Noise),1,length(nAntennasList));  % Noise power
%     SINR = Ptx + DirOK - DirNOK_pcvd - chLoss - problem.Noise;  % Compute SINR Pass-Band (PB)
    SINR_PB = (DirOK_lin*db2pow(chLoss)) ./(DirNOK_gntd_lin*db2pow(chLoss) + Noise_lin);
    SINR_BB = mean(SINRTot,3);  % Compute SINR Base-Band (BB)
    Cap = mean(CapTot,3);  % Compute Average Capacities in the system
    save(fileName,'Cap','SINR_BB','SINR_PB','DirOK','DirNOK_gntd','DirNOK_pcvd','DirOKTot','DirNOKTot','nUsers','nAntennasList');
end







%% EXPERIMENT 7
function experiment7(nIter,nUsersList,nAntennasList,plotFLAG)
    % EXPERIMENT 7 -- 
    % 
    % Aim: Compare the performance of the different subarray allocation
    % techniques in terms of capacity offered. We need to do long GA
    % executions in order to be sure of offering the best solution
    % 
	% Assumptions (Fixed):
    %   1. User location: Fixed, from config file.
    %   2. Antenna Array geometry: Fixed to URA.
    %   3. Algorithm: GA
    % Variable:
    %   1. Antenna Array size variable: nAntennasList
    %   2. Users in system: nUsersList
    %   3. Subarray allocation policy
    % 
    % Syntax:  [] =
    % experiment7(nIters,nUsers,nAntennasList,nUsersList)
    % 
    % Inputs:
    %    nIter - Number of iterations to extract average values
    %    nUsersList - Number of users (set of)
    %    nAntennaList - Number of antenas (set of)
    %
    % Outputs: None
    %
    %------------- BEGIN CODE EXPERIMENT 6 --------------
    %
    fprintf('Running experiment 7...\n');
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');    
    % Override (problem) parameters
    problem.MinObjFIsSNR = true;  % (arbitrary)
    % Override (conf) parameters
    conf.verbosity = 1;
    conf.NumPhaseShifterBits = 0;  % Number of 
    conf.NbitsAmplitude = 0;
    conf.FunctionTolerance_Data = 1e-6;  % Heuristics stops when not improving solution by this much
    conf.multiPath = false;  % LoS channel (for now)
    
    % Override GA parameters
    conf.algorithm = 'GA';  % Heuristic algorithm
    conf.PopulationSize_Data = 150;
    conf.Maxgenerations_Data = 100;
    conf.EliteCount_Data = 25;
    conf.MaxStallgenerations_Data = 40;  % Force it to cover all the generations
    %h1 = figure;
    %hold on
    
    estObj_Free = cell(length(nAntennasList),nIter,length(nUsersList));
    estObj_Loc = cell(length(nAntennasList),nIter,length(nUsersList));
    estObj_Interl = cell(length(nAntennasList),nIter,length(nUsersList));
    estObj_Diag = cell(length(nAntennasList),nIter,length(nUsersList));
    
    % For each case we execute ES and the GA
    for idxUsers = 1:length(nUsersList)
        problem.nUsers = nUsersList(idxUsers);  % Number of users in the simulation
        % Configure basic parameters
        candSet = (1:1:problem.nUsers);  % Set of users to be considered
        problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
        for idxAnt = 1:length(nAntennasList)
            for idxIter = 1:nIter
                fprintf('Iteration %d with nUsers %d and nAntenas %d\n',...
                    idxIter,nUsersList(idxUsers),nAntennasList(idxAnt));
                % Configure the simulation environment. Need to place users in new
                % locations (if not fixed) and create new channels 
                % to have statistically meaningful results (if not LoS)
                [problem_temp,~,~] = f_configuration(conf,problem);
                % Select number of antennas
                problem_temp.N_Antennas = nAntennasList(idxAnt);
                % Adjust parameters
                problem_temp.NxPatch = floor(sqrt(problem_temp.N_Antennas));
                problem_temp.NyPatch = floor(problem_temp.N_Antennas./problem_temp.NxPatch);
                problem_temp.N_Antennas = problem_temp.NxPatch.*problem_temp.NyPatch;
                % Call heuristics
                fprintf('\t** %d Antennas and %d Users...\n',problem_temp.N_Antennas,problem_temp.nUsers);
                for aRestr = {'None','Localized','Interleaved','DiagInterleaved'}
                    conf.arrayRestriction = aRestr;
                    [~,~,~,estObj] = f_heuristics(problem_temp,conf,candSet);
                    switch cell2mat(aRestr)
                        case 'None'
                            estObj_Free(idxAnt,idxIter,idxUsers) = mat2cell(estObj,1);
                        case 'Localized'
                            estObj_Loc(idxAnt,idxIter,idxUsers) = mat2cell(estObj,1);
                        case 'Interleaved'
                            estObj_Interl(idxAnt,idxIter,idxUsers) = mat2cell(estObj,1);
                        case 'DiagInterleaved'
                            estObj_Diag(idxAnt,idxIter,idxUsers) = mat2cell(estObj,1);
                    end
                end
                save('temp/exp7-results_so_far','estObj_Free','estObj_Loc',...
                            'estObj_Interl','estObj_Diag','nUsersList','nAntennasList');
            end
        end
    end
end







%% EXPERIMENT 8
function experiment8(nIter,nUsersList,nAntennasList,plotFLAG)
    % EXPERIMENT 8 -- 
    % 
    % Aim: Compare the performance of our heuristic approach with a random 
    % modification of itself, aiming at analyzing the performance of our
    % result.
    % 
	% Assumptions (Fixed):
    %   1. User location: Fixed, from config file.
    %   2. Antenna Array geometry: Fixed to URA.
    %   3. Algorithm: GA
    %   3. Subarray allocation policy: Free
    % Variable:
    %   1. Antenna Array size variable: nAntennasList
    %   2. Users in system: nUsersList
    % 
    % Syntax:  [] =
    % experiment8(nIters,nUsers,nAntennasList,nUsersList)
    % 
    % Inputs:
    %    nIter - Number of iterations to extract average values
    %    nUsersList - Number of users (set of)
    %    nAntennaList - Number of antenas (set of)
    %
    % Outputs: None
    %
    %------------- BEGIN CODE EXPERIMENT 8 --------------
    %
    fprintf('Running experiment 8...\n');
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');    
    % Override (problem) parameters
    problem.MinObjFIsSNR = true;  % (arbitrary)
    % Override (conf) parameters
    conf.verbosity = 1;
    conf.NumPhaseShifterBits = 0;  % Number of 
    conf.NbitsAmplitude = 0;
    conf.FunctionTolerance_Data = 1e-6;  % Heuristics stops when not improving solution by this much
    conf.multiPath = false;  % LoS channel (for now)
    
    % Override GA parameters
    conf.PopulationSize_Data = 15;
    conf.Maxgenerations_Data = 10;
    conf.EliteCount_Data = 2;
    conf.MaxStallgenerations_Data = 40;  % Force it to cover all the generations
    %h1 = figure;
    %hold on
    
    estObj_heur = cell(length(nAntennasList),nIter,length(nUsersList));
    estObj_rnd = cell(length(nAntennasList),nIter,length(nUsersList));
    
    % For each case we execute ES and the GA
    for idxUsers = 1:length(nUsersList)
        problem.nUsers = nUsersList(idxUsers);  % Number of users in the simulation
        % Configure basic parameters
        candSet = (1:1:problem.nUsers);  % Set of users to be considered
        problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
        for idxAnt = 1:length(nAntennasList)
            for idxIter = 1:nIter
                fprintf('Iteration %d with nUsers %d and nAntenas %d\n',...
                    idxIter,nUsersList(idxUsers),nAntennasList(idxAnt));
                % Configure the simulation environment. Need to place users in new
                % locations (if not fixed) and create new channels 
                % to have statistically meaningful results (if not LoS)
                [problem_temp,~,~] = f_configuration(conf,problem);
                % Select number of antennas
                problem_temp.N_Antennas = nAntennasList(idxAnt);
                % Adjust parameters
                problem_temp.NxPatch = floor(sqrt(problem_temp.N_Antennas));
                problem_temp.NyPatch = floor(problem_temp.N_Antennas./problem_temp.NxPatch);
                problem_temp.N_Antennas = problem_temp.NxPatch.*problem_temp.NyPatch;
                % Call heuristics
                fprintf('\t** %d Antennas and %d Users...\n',problem_temp.N_Antennas,problem_temp.nUsers);
                conf.algorithm = 'GA';  % Heuristic algorithm
                [~,~,~,estObj] = f_heuristics(problem_temp,conf,candSet);
                estObj_heur(idxAnt,idxIter,idxUsers) = mat2cell(estObj,1);
                conf.algorithm = 'GA-rnd';  % Heuristic algorithm
                [~,~,~,estObj] = f_heuristics(problem_temp,conf,candSet);
                estObj_rnd(idxAnt,idxIter,idxUsers) = mat2cell(estObj,1);
                save('temp/exp8-results_so_far','estObj_heur',...
                    'estObj_rnd','nUsersList','nAntennasList');
            end
        end
    end
end