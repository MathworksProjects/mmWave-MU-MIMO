%% Setup the environment
clear; clc; close all;
addpath('utilities','-end');  % Add utilities folder at the end of search path

%% Define several experiments here and override variable values accordingly
experimentList = [31];

%% Experiment selection

%% EXPERIMENT 1
if any(experimentList(:)==1)
    input.nUsersList           = [2 4 6 8 10];
    input.payloadList          = (1/8).*(50/10).*[385e6 962.5e6 1155e6 1540e6 1925e6 2502e6 3080e6 4620e6].*10e-3;
    % Input parameters
    input.nIter                = 1;
    input.nAntennas            = 4.^2;
    input.arrayRestriction     = 'None';
    input.algorithm            = 'GA';
    input.detLocation          = true;
    input.useCasesLocation     = true;
    input.useCaseLocation      = 4;
    plotFLAG                   = false;
    for payload = input.payloadList
        input.payload = payload;
        for nUsers = input.nUsersList
            input.nUsers = nUsers;
            experiment1(input);
        end
    end
end

%% EXPERIMENT 1 - PLOTTING
if any(experimentList(:)==11)
    input.nUsersList           = [2 4 6 8 10];
    input.payloadList          = (1/8).*(50/10).*[385e6 962.5e6 1155e6 1540e6 1925e6 2502e6 3080e6 4620e6].*10e-3;
    input.nAntennas            = 4.^2;
    input.arrayRestriction     = 'None';
    input.algorithm            = 'GA';
    input.detLocation          = true;
    input.useCasesLocation     = true;
    input.useCaseLocation	   = 4;
    experiment1_plot(input);
end

%% EXPERIMENT 2
if any(experimentList(:)==2)
%     arrRestctList    = {'None','Localized'};
    arrRestctList = {'None'};
    % Input parameters
    input.nIter               = 10;  % Total number of iterations
    input.nUsers              = 2;  % Number of users deployed
    input.nAntennasList       = [4 8 12 16].^2;  % Number of antennas in array
    input.Maxgenerations_Data = 100;
    input.algorithm           = 'GA';  % Heuristic algorithm
    input.detLocation         = true;  % Deterministic locations if true
    input.useCasesLocation    = true;  % Use-Case locations if true
    useCaseLocationList = 1;
    for locID = 1:length(useCaseLocationList)
        input.useCaseLocation     = useCaseLocationList(locID);  % Use-case ID
        fileNameList = cell(length(arrRestctList),1);
        score_tot = zeros(length(arrRestctList),length(input.nAntennasList),input.Maxgenerations_Data);
        DirOK_tot = zeros(length(arrRestctList),length(input.nAntennasList));
        DirNOK_tot = zeros(length(arrRestctList),length(input.nAntennasList));
        for restIdx = 1:length(arrRestctList)
            % Input parameters extra
            input.arrayRestriction = arrRestctList{restIdx};
            % Main
            fileNameList{restIdx} = experiment2(input);
            % Parse results
            load(fileNameList{restIdx},'bestScores','DirOKTot','DirNOKTot');
            score_tot(restIdx,:,:) = bestScores;
            DirOK_tot(restIdx,:) = DirOKTot;
            DirNOK_tot(restIdx,:) = DirNOKTot;
        end
        % Save results
        fileName = strcat('temp/exp2_',input.algorithm,'_TOT_',mat2str(input.nUsers),'_',mat2str(input.detLocation),'_',mat2str(input.useCasesLocation),'_',mat2str(input.useCaseLocation));
        nAntennasList = input.nAntennasList; Maxgenerations_Data=input.Maxgenerations_Data;
        save(fileName,'score_tot','DirOK_tot','DirNOK_tot','nAntennasList','arrRestctList','Maxgenerations_Data');
        % Plot results
        experiment2_plot(fileName);
    end
end

%% EXPERIMENT 2 - PLOTTING
if any(experimentList(:)==21)
    % Get results by parameters - individual
    algorithm = 'GA';
    nUsers = 2;
    useCaseLocation = 4;
%     fileName = strcat('temp/exp2_',algorithm,'_TOT_',mat2str(nUsers),'_true_true_',mat2str(useCaseLocation));
%     experiment2_plot(fileName);
    % Call plot
    useCaseLocationList = [1 2 3 4 5 6];
	experiment2_plot(algorithm,nUsers,useCaseLocationList);
end

%% EXPERIMENT 3
if any(experimentList(:)==3)
    input.nUsersList           = [2 4 6 8 10];
    input.deadlineList         = [20 30 40 50 60 70 80 90 100];
    % Input parameters
    input.nIter                = 1;
    input.nAntennas            = 4.^2;
    input.arrayRestriction     = 'None';
    input.algorithm            = 'GA';
    input.detLocation          = true;
    input.useCasesLocation     = true;
    input.useCaseLocation      = 4;
    plotFLAG                   = false;
    for deadline = input.deadlineList
        input.deadline = deadline;
        input.payload = (1/8).*(deadline/10).*1925e6.*10e-3;
        for nUsers = input.nUsersList
            input.nUsers = nUsers;
            experiment3(input);
        end
    end
end

%% EXPERIMENT 3 - PLOTTING
if any(experimentList(:)==31)
    input.nUsersList           = [2 4 6 8 10];
    input.deadlineList         = [20 30 40 50 60 70 80 90 100];
    input.nAntennas            = 4.^2;
    input.payload              = (1/8).*(50/10).*1925e6.*10e-3;
    input.arrayRestriction     = 'None';
    input.algorithm            = 'GA';
    input.detLocation          = true;
    input.useCasesLocation     = true;
    input.useCaseLocation	   = 4;
    experiment3_plot(input);
end

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
    for restIdx = 1:length(arrRestctList)
        % Input parameters
        input.nIter            = 5;  % Total number of iterations
        input.nUsers           = 2;  % Number of users deployed
        input.nAntennasList    = [4 8 12 16 20 24 28 32].^2;  % Number of antennas in array
        input.arrayRestriction = arrRestctList{restIdx};
        input.algorithm        = 'GA';  % Heuristic algorithm
        input.detLocation      = true;  % Deterministic locations if true
        input.useCasesLocation = true;  % Use-Case locations if true
        input.useCaseLocation  = 3;  % Use-case ID
        plotFLAG = false;  % Plotting flag
        % Main
        fileNameList{restIdx} = experiment5(input,plotFLAG);
        % Plot
        experiment5_plot(fileNameList{restIdx});
    end
    % Parse results along array geometry
    SINR_PB_tot = zeros(length(input.nAntennasList),length(arrRestctList));
    SINR_BB_tot = zeros(length(input.nAntennasList),length(arrRestctList));
    Cap_tot = zeros(length(input.nAntennasList),length(arrRestctList));
    for restIdx = 1:length(arrRestctList)
        load(fileNameList{restIdx});
        % Parse results
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
    nUsers = 2;  % Number of users considered in the simulation
    arrRestctList = {'None','Localized'};  % List of sub-array restrictions to consider
%     my_nAntennasList = [4 8 12 16].^2;  % List of antennas to plot over
    my_nAntennasList    = [4 8 12 16 20 24 28 32].^2;  % Number of antennas in array
    locList = [1 2 3 4 5 6];  % List of locations to plot over
    % Call plot
    experiment5_plot(nUsers,arrRestctList,locList,my_nAntennasList);
end

%% EXPERIMENT 6
if any(experimentList(:)==6)
    nUsersList = [2 4 6];
    % Input parameters
    input.nIter                = 1;
    input.nAntennas            = 12.^2;
    input.arrayRestriction     = 'None';
    input.algorithm            = 'GA';
    input.detLocation          = true;
    input.useCasesLocation     = true;
    input.useCaseLocationList  = 1;  % List of locations to plot over;
    plotFLAG                   = false;  % Plotting flag
    for nUsers = nUsersList
        input.nUsers = nUsers;
        experiment6(input,plotFLAG);
    end
end

%% EXPERIMENT 6 - PLOTTING
if any(experimentList(:)==61)
    nUsersList = [2 4 6];
    input.nAntennas            = 12.^2;
    input.arrRestct            = 'None';
    input.algorithm            = 'GA';
    input.detLocation          = true;
    input.useCasesLocation     = true;
    input.useCaseLocationList  = 1;  % List of locations to plot over;
    % Call plot
    experiment6_plot(nUsersList,input);
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

%% EXPERIMENT 9
if any(experimentList(:)==9)
%     nUsersList = [2 4 6 8 10 12];
    nUsersList = [2 4 6];
    % Input parameters
    input.nIter                = 1;
%     input.nAntennasList        = [12 14 16 18 20 22].^2;
    input.nAntennasList        = [6].^2;
    input.arrayRestriction     = 'None';
    input.algorithm            = 'GA';
    input.detLocation          = true;
    input.useCasesLocation     = true;
    input.useCaseLocation      = 4;  % List of locations to plot over;
    plotFLAG                   = false;  % Plotting flag
    for nUsers = nUsersList
        input.nUsers = nUsers;
        experiment9(input,plotFLAG);
    end
end

%% EXPERIMENT 9 - PLOTTING
if any(experimentList(:)==91)
%     nUsersList = [2 4 6 8 10 12];
%     nUsersList = [2 4 6 8];
    nUsersList = [2 4 6];
%     input.nAntennasList        = [12 14 16 18 20 22].^2;
%     input.nAntennasList        = [8 10].^2;
    input.nAntennasList        = [6].^2;
    input.arrRestct            = 'None';
    input.algorithm            = 'GA';
    input.detLocation          = true;
    input.useCasesLocation     = true;
    input.useCaseLocation      = 4;  % List of locations to plot over;
    % Call plot
    experiment9_plot(nUsersList,input);
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
function experiment1(input)
    % EXPERIMENT 1 - Capacity offered
    % In this experiment we evaluate the capacity that the heuristics are able
    % to offer to the devices. Heuristics assigns antennas as a function of the
    % priority. The traffic is overloaded. The users location and channel
    % varies across simulations.
    %
    %------------- BEGIN CODE EXPERIMENT 1 --------------
    %
%     fprintf('Running experiment 1...\n');
%     fprintf('Input parameters:\n\tnUsers:\t%d\n\tIter:\t%d\n\tAntennas:\t%d\n\tarrayRestriction:\t%s\n\talgorithm:\t%s\n\tdetLocation:\t%s\n\tuseCasesLocation:\t%s\n\tuseCaseLocation:\t%d\n',input.nUsers,input.nIter,input.nAntennas,input.arrayRestriction,input.algorithm,mat2str(input.detLocation),mat2str(input.useCasesLocation),input.useCaseLocation);
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override parameters (config)
    conf.verbosity            = 0;
    conf.algorithm            = input.algorithm;
	conf.detLocation          = input.detLocation;
    conf.useCasesLocation     = input.useCasesLocation;
    conf.useCaseLocation      = input.useCaseLocation;
    % Override parameters (problem)
    problem.nUsers            = input.nUsers;
    problem.N_Antennas        = input.nAntennas;
    problem.payload           = input.payload;
    problem.arrayRestriction  = input.arrayRestriction;
    problem.DEBUG             = true;  % no ouput messages (Scheduling)
    problem.PLOT_DEBUG        = false;  % No plotting
    problem.numPkts           = 2;
    % Output variables
    ratioOKTot = zeros(input.nUsers,input.nIter);
    ratioNOKTot = zeros(input.nUsers,input.nIter);
    for iter = 1:input.nIter
        % Configure the simulation environment
        [problem,~,flows] = f_configuration(conf,problem);
        % Main function
        [flows,~,~,~,~,~] = main(conf,problem,flows);
        % Report of single execution
        [ratioOKTot(:,iter),ratioNOKTot(:,iter)] = f_generateReport(flows,problem.DEBUG);
        fprintf('\tNusers %d -> OK=%.2f(%%)\n',input.nUsers,mean(ratioOKTot(:,iter)));
    end
    fprintf('Nusers %d -> OK = %.2f(%%)\n',input.nUsers,mean(mean(ratioOKTot,2)));
    % Store results in temp file (to be retrieved by plotting function)
    fileName = strcat('temp/exp1_',conf.algorithm,'_',problem.arrayRestriction,'_',mat2str(problem.nUsers),'_',mat2str(problem.payload),'_',mat2str(problem.N_Antennas),'_',mat2str(conf.detLocation),'_',mat2str(conf.useCasesLocation),'_',mat2str(conf.useCaseLocation));
    % Store in local to store in file
    ratioOK          = mean(ratioOKTot,2);  %#ok
    ratioNOK         = mean(ratioNOKTot,2);  %#ok
    nUsers           = problem.nUsers;  %#ok
    nAntennas        = problem.N_Antennas;  %#ok
    payload          = problem.payload;  %#ok
    arrayRestriction = problem.arrayRestriction;  %#ok
    detLocation      = conf.detLocation;  %#ok
    useCaseLocation  = conf.useCaseLocation;  %#ok
    useCaseLocation  = conf.useCaseLocation;  %#ok
    save(fileName,'ratioOK','ratioNOK',...
         'nUsers','payload','arrayRestriction',...
         'detLocation','useCaseLocation','useCaseLocation');
end



%% EXPERIMENT 1 - PLOTTING
function experiment1_plot(input)
    h =  findobj('type','figure');
    figIdx = length(h) + 1;
    ratioOK_av = zeros(length(input.nUsersList),1);
    ratioNOK_av = zeros(length(input.nUsersList),1);
    leg = cell(length(input.payloadList),1);
    LineStyleList = {'-.','-','--',':'};
    MarkerList = {'s','*','d','+'};
    for payload = input.payloadList
        for nUsers = input.nUsersList
            fileName = strcat('temp/exp1_',input.algorithm,'_',...
                                           input.arrayRestriction,'_',...
                                           mat2str(nUsers),'_',...
                                           mat2str(payload),'_',...
                                           mat2str(input.nAntennas),'_',...
                                           mat2str(input.detLocation),'_',...
                                           mat2str(input.useCasesLocation),'_',...
                                           mat2str(input.useCaseLocation));
            load(fileName,'ratioOK','ratioNOK');
            ratioOK_av(nUsers==input.nUsersList) = mean(ratioOK);
            ratioNOK_av(nUsers==input.nUsersList) = mean(ratioNOK);
        end
        idxMarkerStyle= mod(find(payload==input.payloadList),4)+1;
        idxLineStyle = ceil(find(payload==input.payloadList)/4);
        % No interpolation
        figure(figIdx); hold on; grid minor;
        plot(input.nUsersList,ratioOK_av,'LineWidth',1,'LineStyle',LineStyleList{idxLineStyle},'Marker',MarkerList{idxMarkerStyle},'Color','k');
        % With interpolation
        figure(figIdx+1); hold on; grid minor;
        xq = (1:1:input.nUsersList(end));
        ratioOK_av_int = interp1(input.nUsersList,ratioOK_av,xq,'pchip');
        plot(xq,ratioOK_av_int,'LineWidth',1,'LineStyle',LineStyleList{idxLineStyle},'Marker',MarkerList{idxMarkerStyle},'Color','k');
        % Legend
        leg(payload==input.payloadList) = strcat('Network sat.',{' '},sprintf('%.2f',payload/max(input.payloadList)),'x');
    end
    figure(figIdx); grid minor;
    xlabel('Number of users','FontSize',12);
    ylabel('Ratio OK (%)','FontSize',12);
    title('Ratio of application that meet deadline (%)','FontSize',12);
    hleg = legend(leg);
    set(hleg,'FontSize',10);
    figure(figIdx+1); grid minor;
    xlabel('Number of users','FontSize',12);
    ylabel('Ratio OK (%)','FontSize',12);
    title('Ratio of application that meet deadline (%)','FontSize',12);
    hleg = legend(leg);
    set(hleg,'FontSize',10);
end







%% EXPERIMENT 2 - Convergency analysis
% EXPERIMENT 2 - Chances of achieving the demanded throughput
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
    conf.verbosity = 0;
    conf.algorithm = algorithm;  % Heuristic algorithm
    conf.PopulationSize_Data = 30;
    conf.Maxgenerations_Data = Maxgenerations_Data;
    conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);
    conf.MaxStallgenerations_Data = conf.Maxgenerations_Data;  % Force it to cover all the generations
    conf.FunctionTolerance_Data = 1e-10;  % Heuristics stops when not improving solution by this much
    conf.NumPhaseShifterBits = 60;  % Number of bits to control heuristic solution
    conf.multiPath = false;  % LoS channel (for now)
    conf.detLocation = detLocation;  % Use fixed pre-stored locations
    conf.useCasesLocation = useCasesLocation;  % Use the use-case locations
    conf.useCaseLocation = useCaseLocation;  % Specify the use-case location
    % Output parameters
    bestScores1 = zeros(length(nAntennasList),nIter,conf.Maxgenerations_Data);  % temp
    bestScores = zeros(length(nAntennasList),conf.Maxgenerations_Data);  % final
    DirOKTot = -Inf(1,length(nAntennasList));  % final
    DirNOKTot = -Inf(nUsers-1,length(nAntennasList));  % final
    fileName = strcat('temp/exp2_',conf.algorithm,'_',problem.arrayRestriction,'_',mat2str(nUsers),'_',mat2str(conf.detLocation),'_',mat2str(conf.useCasesLocation),'_',mat2str(conf.useCaseLocation));
    % For each case we execute ES and the GA
    for idxAnt = 1:length(nAntennasList)
        DirOK1 = -Inf(1,nIter);
        DirNOK1 = -Inf(nUsers,nIter);
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
                'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque s\ED
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
            [~,W,~,~,tempBS] = ...
                o_solveSingleNmaxUserInstance(conf,problem,...
                problem.NmaxArray(problem.IDUserAssigned));
            bestScores1(idxAnt,idxIter,1:length(tempBS)) = tempBS;
            % Smoothen out results in tail
            if length(tempBS)<conf.Maxgenerations_Data
                stillCover = length(tempBS)+1:conf.Maxgenerations_Data;
                bestScores1(idxAnt,idxIter,stillCover) = repmat(bestScores1(idxAnt,idxIter,length(tempBS)),1,length(stillCover));
            end
            % Extract Directivity (intended user and towards others)
            id = 1;  % We only evaluate user 1
            problem.ant_elem = sum(W(id,:)~=0);
            relevant_positions = (W(id,:)~=0);
            Taper_user = W(id,relevant_positions);
            handle_Conf_Array = phased.ConformalArray('Element',problem.handle_Ant,...
                                  'ElementPosition',...
                                  [zeros(1,problem.ant_elem);...
                                  problem.possible_locations(2,relevant_positions);...
                                  problem.possible_locations(3,relevant_positions)],...
                                  'Taper',Taper_user);
            % Extract Rx Power (in dB)
            DirOK1(1,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id),'Azimuth',problem.phiUsers(id),'Type','powerdb');
            fprintf('* Directivity IDmax: %.2f (dB)\n',DirOK1(1,idxIter));
            % Extract interference generated to others (in dB)
            for id1 = 1:1:problem.nUsers
                if id1~=id
                    DirNOK1(id1,idxIter) = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(id1),'Azimuth',problem.phiUsers(id1),'Type','powerdb');
                    fprintf('  Directivity IDmin(%d): %.2f (dB)\n',id1,DirNOK1(id1,idxIter));
                end
            end
        end
        t = mean(bestScores1(idxAnt,:,:),2);
        bestScores(idxAnt,:) = t(:).';
        DirOKTot(1,idxAnt) = pow2db(mean(db2pow(DirOK1)));
        DirNOKTot(:,idxAnt) = pow2db(mean(db2pow(DirNOK1(2:end,:)),2));
    end
    save(fileName,'bestScores','DirOKTot','DirNOKTot',...
             'nUsers','nAntennasList','arrayRestriction',...
             'detLocation','useCaseLocation','useCaseLocation');
end

%% EXPERIMENT 2 - Plotting
function experiment2_plot(varargin)
    h =  findobj('type','figure');
    figIdx = length(h) + 1;
    FigSize = [530 260];
    colNoneList = [[139 0 0]./255 ; [178 34 34]./255 ; [220 20 60]./255 ; [205 92 92]./255 ; [250 128 114]./255 ; [255 160 122]./255];
    colLocList  = [[0 0 128]./255 ; [0 0 205]./255 ; [0 0 255]./255 ; [65 105 225]./255 ; [100 149 237]./255 ; [0 191 255]./255];
    MarkerList = {'s','*','o','x','d','p','h','^','v','>','<'};
    if nargin == 1
        % Plot individual
        fileName = varargin{1};
        % Load variables
        load(fileName,'score_tot','DirOK_tot','DirNOK_tot','nAntennasList','arrRestctList');
        parse = strsplit(fileName,'_');
        userLocation = parse{end};
        nUsers = parse{4};
        Maxgenerations_Data = size(score_tot,3);  %#ok
        % Configure plots 1
        colorList = [[0 0 255]./255 ; [0 51 102]./255];
    %     lineStyleList = {'-s','-.o','--x',':*'};
        lineStyleList = {'-','-.','--',':'};
        legendList = cell(length(arrRestctList)*length(nAntennasList),1);  %#ok
        % Parse aditional results
        for antIdx = 1:length(nAntennasList)
            for restIdx = 1:length(arrRestctList)
                idx = (antIdx-1)*length(arrRestctList) + restIdx;
                res = score_tot(restIdx,antIdx,:);  %#ok
                score_TOT(:,idx) = res(:);  %#ok
                legendList(idx) = strcat(arrRestctList{restIdx},{' '},'-',{' '},mat2str(nAntennasList(antIdx)));  %#ok
            end
        end
        % Plot 1
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        set(gca, 'ColorOrder', colorList, 'linestyleorder', lineStyleList, 'NextPlot', 'replacechildren');  % Change to new colors.
        plot(1:Maxgenerations_Data,score_TOT,'LineWidth',1.5,'MarkerSize',4);
        ylim([-1 -0.05]);
        grid minor;
        xlabel('Generations','FontSize',12);
        ylabel('Fitness value','FontSize',12);
        titl = strcat('Convergency analysis for',{' '},nUsers,{' '},'users - Location',{' '},userLocation);
        title(titl,'FontSize',12);
        hList = findobj('Type', 'line');  % Returns lines in inverse order
        hList = hList(end:-1:1);
        ah1 = gca;
        leg = legend(ah1,hList(1:2:2*length(nAntennasList)),legendList(1:2:2*length(nAntennasList)));
        set(leg,'FontSize',8);
        ah2 = axes('position',get(gca,'position'),'visible','off');
        leg = legend(ah2,hList(2:2:2*length(nAntennasList)),legendList(2:2:2*length(nAntennasList)));
        set(leg,'FontSize',8);
        set(fHandle,'Position',[fHandle.Position(1) fHandle.Position(2) FigSize(1) FigSize(2)]);

        % Configure plot Directivities 2 (dB)
        colorList = [[0 0 255]./255;[0 51 102]./255;[255 0 0]./255;[139 0 0]./255];
        FigSize = [530 260];
        % Plot Directivities (dB)
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        set(gca, 'ColorOrder', colorList, 'NextPlot', 'replacechildren');  % Change to new colors.
        data(:,1,1) = DirOK_tot(1,:).';  %#ok
        data(:,2,1) = DirOK_tot(2,:).';
        data(:,1,2) = -DirNOK_tot(1,:).';  %#ok
        data(:,2,2) = -DirNOK_tot(2,:).';
        Xneg = data;
        Xneg(Xneg>0) = 0;
        Xpos = data;
        Xpos(Xpos<0) = 0;
        groupLabels = nAntennasList;
        plotBarStackGroups(Xpos, groupLabels);
        plotBarStackGroups(Xneg, groupLabels);
        ylim([-10 90]);
        xlim([1-0.5 length(nAntennasList)+0.5])
        xticklabels(nAntennasList);
        grid minor;
        titl = strcat(nUsers,{' '},'users - Location',{' '},userLocation);
        title(titl,'FontSize',11)
        ylabel('Directivity (dB)','FontSize',11)
        xlabel('Number of antennas','FontSize',11)
        leg = legend('No Restr. - target','No Restr. - others','Restr. - target','Restr. - others');
        set(leg,'FontSize',8,'Location','NorthWest');
        set(fHandle,'Position',[fHandle.Position(1) fHandle.Position(2) FigSize(1) FigSize(2)]);
        nUsers = str2double(nUsers);  % Parse it for future plots
    elseif nargin==3
        % Plot list
        algorithm = varargin{1};
        nUsers = varargin{2};
        useCaseLocationList = varargin{3};
        legendList = cell(4,1);
        
        for idx = 1:length(useCaseLocationList)
            userLocation = useCaseLocationList(idx);
            fileName = strcat('temp/exp2_',algorithm,'_TOT_',mat2str(nUsers),'_true_true_',mat2str(userLocation));
            load(fileName,'score_tot','DirOK_tot','DirNOK_tot','nAntennasList','arrRestctList');
            DirOK_None(idx,:) = DirOK_tot(1,:);  %#ok
            DirOK_Localized(idx,:) = DirOK_tot(2,:);  %#ok
            DirNOK_None(idx,:) = DirNOK_tot(1,:);  %#ok
            DirNOK_Localized(idx,:) = DirNOK_tot(2,:);  %#ok
            data_OK_None(:,idx,1) = DirOK_tot(1,:).';  %#ok
            data_NOK_None(:,idx,1) = DirNOK_tot(1,:).';  %#ok
            data_OK_Localized(:,idx,1) = DirOK_tot(2,:).';  %#ok
            data_NOK_Localized(:,idx,1) = DirNOK_tot(2,:).';  %#ok
            % Correction for higher values
            DirNOK_None(DirNOK_None<-100) = -100;  %#ok
            DirNOK_Localized(DirNOK_Localized<-100) = -100;  %#ok
            data_NOK_None(data_NOK_None<-100) = -100;  %#ok
            data_NOK_Localized(data_NOK_Localized<-100) = -100;  %#ok
            legendList(idx) = strcat('Localization -',{' '},mat2str(userLocation));
            
            % Plot
            fHandle = figure(figIdx); hold on;
            plot(nAntennasList,DirOK_tot(1,:),'Marker',MarkerList{idx},'LineStyle','-','Color','k');
            plot(nAntennasList,DirOK_tot(2,:),'Marker',MarkerList{idx},'LineStyle','--','Color','k');
            fHandle = figure(figIdx+1); hold on; 
            plot(nAntennasList,DirNOK_tot(1,:),'Marker',MarkerList{idx},'LineStyle','-','Color','k');
            plot(nAntennasList,DirNOK_tot(2,:),'Marker',MarkerList{idx},'LineStyle','--','Color','k');
        end
        fHandle = figure(figIdx); grid minor
        xlabel('Number of antennas','FontSize',12);
        ylabel('Directivity (dB)','FontSize',12);
        title('Directivity towards intended user','FontSize',12);
        leg = legend(legendList);
        set(leg,'FontSize',10,'Location','NorthWest');
        xlim([min(nAntennasList)-5 max(nAntennasList)+5]);
        fHandle = figure(figIdx+1); grid minor
        xlabel('Number of antennas','FontSize',12);
        ylabel('Directivity (dB)','FontSize',12);
        title('Directivity towards other user','FontSize',12);
        leg = legend(legendList);
        set(leg,'FontSize',10,'Location','NorthWest');
        xlim([min(nAntennasList)-5 max(nAntennasList)+5]);
        figIdx = figIdx + 2;
        % Plot No geometry restriction - OK
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        groupLabels = nAntennasList;
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        set(gca, 'ColorOrder', colNoneList, 'NextPlot', 'replacechildren');  % Change to new colors
        plotBarStackGroups(data_OK_None, groupLabels);
        grid minor;
        ylabel('Directivity (dB)','FontSize',12)
        xlabel('Number of antennas','FontSize',12)
        title('No restriction - target','FontSize',14)
        leg = legend(legendList);
        set(leg,'FontSize',10,'Location','NorthWest');
        % Plot No geometry restriction - NOK
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        set(gca, 'ColorOrder', colNoneList, 'NextPlot', 'replacechildren');  % Change to new colors
        plotBarStackGroups(-data_NOK_None, groupLabels);
        grid minor;
        ylabel('Directivity (dB)','FontSize',12)
        xlabel('Number of antennas','FontSize',12)
        title('No restriction - others','FontSize',14)
        leg = legend(legendList);
        set(leg,'FontSize',10,'Location','NorthWest');
        yt = yticklabels;
        yt = -str2double(yt);
        yt = num2cell(yt);
        yticklabels(yt);
        % Plot Localized geometry - OK
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        groupLabels = nAntennasList;
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        set(gca, 'ColorOrder', colLocList, 'NextPlot', 'replacechildren');  % Change to new colors
        plotBarStackGroups(data_OK_Localized, groupLabels);
        grid minor;
        ylabel('Directivity (dB)','FontSize',11)
        xlabel('Number of antennas','FontSize',11)
        title('Localized restriction - target','FontSize',12)
        leg = legend(legendList);
        set(leg,'FontSize',10,'Location','NorthWest');
        % Plot Localized geometry - NOK
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        set(gca, 'ColorOrder', colLocList, 'NextPlot', 'replacechildren');  % Change to new colors
        plotBarStackGroups(-data_NOK_Localized, groupLabels);
        grid minor;
        ylabel('Directivity (dB)','FontSize',12)
        xlabel('Number of antennas','FontSize',12)
        title('Localized restriction - others','FontSize',14)
        leg = legend(legendList);
        set(leg,'FontSize',10,'Location','NorthWest');
        yt = yticklabels;
        yt = -str2double(yt);
        yt = num2cell(yt);
        yticklabels(yt);
        
        MaxGenerations = size(score_tot,3);  %#ok
        tNon = zeros(MaxGenerations,4);
        tLoc = zeros(MaxGenerations,4);
        fHandle = figure(figIdx); hold on; figIdx = figIdx + 1;
        MarkerList = {'s','*','+','o'};
        for p = 1:4
            t = score_tot(1,p,:);
            tNon(:,p) = t(:);
            p1(p) = line_fewer_markers(1:MaxGenerations,tNon(:,p),20,'Color','k','lineStyle','-','Marker',MarkerList{p},'MarkerSize',7);
            t = score_tot(2,p,:);
            tLoc(:,p) = t(:);
            p2(p) = line_fewer_markers(1:MaxGenerations,tLoc(:,p),20,'Color','k','lineStyle',':','Marker',MarkerList{p},'MarkerSize',7);
        end
        xlabel('Generations','FontSize',12);
        ylabel('Fitness value','FontSize',12);
        title('Convergency analysis in Heuristics-BF','FontSize',12);
        ah1 = gca;
        legend1 = {'No restr. - 16 antennas','No restr. - 64 antennas','No restr. - 144 antennas','No restr. - 256 antennas'};
        legend2 = {'Localized - 16 antennas','Localized - 64 antennas','Localized - 144 antennas','Localized - 256 antennas'};
        leg = legend(ah1,p1,legend1);
        set(leg,'FontSize',10);
        ah2 = axes('position',get(gca,'position'),'visible','off');
        leg = legend(ah2,p2,legend2);
        set(leg,'FontSize',10);
        
    else
        error('wrong number of parameters');
    end

    % Print use case locations
    nUsers = 50;  % For better visualization
    problem.nUsers = nUsers;
    Delta = 90/problem.nUsers;
    if mod(problem.nUsers,2)==0  % even
        ini = Delta/2;
        vect = [(-problem.nUsers/2:1:-1) (1:1:problem.nUsers/2)];
    else                         % odd
        ini = Delta;
        vect = [(-floor(problem.nUsers/2):1:-1) 0 (1:1:floor(problem.nUsers/2))];
    end
    vect = ini.*vect;
    % UC 1 - Located horizontally (no elevation)
    uc_el(1,:) = zeros(1,problem.nUsers);
    uc_az(1,:) = vect;
    uc_dist(1,:) = 5.*ones(1,problem.nUsers);
    % UC 2 - Located horizontally (no elevation - a bit more separation)
    uc_el(2,:) = zeros(1,problem.nUsers);
    uc_az(2,:) = 1.5.*vect;
    uc_dist(2,:) = 5.*ones(1,problem.nUsers);
    % UC 3 - Located vertically (no azymuth)
    uc_el(3,:) = vect;
    uc_az(3,:) = zeros(1,problem.nUsers);
    uc_dist(3,:) = 5.*ones(1,problem.nUsers);
    % UC 4 - Located diagonaly
    uc_el(4,:) = vect;
    uc_az(4,:) = vect;
    uc_dist(4,:) = 5.*ones(1,problem.nUsers);
    % UC 5 - Located vertically (15 deg azymuth)
    uc_el(5,:) = vect;
    uc_az(5,:) = 15.*ones(1,problem.nUsers);
    uc_dist(5,:) = 5.*ones(1,problem.nUsers);
    % UC 6 - Located horizontally (15 deg elevation)
    uc_el(6,:) = 15.*ones(1,problem.nUsers);
    uc_az(6,:) = vect;
    uc_dist(6,:) = 5.*ones(1,problem.nUsers);
    
    % Visualize user location/channels
    dx = 0.1;   dy = 0.1;   dz = 0.3;
    colorList = {'r','b','g','k','c','m'};
    IDmax = 1;  % selected user for exp2 is 1
    legendList = cell(6,1);
    figure(figIdx); hold on; figIdx = figIdx + 1;
    grid minor;
    for location = 1:1:6
        problem.phiUsers = uc_az(location,:);
        problem.thetaUsers = uc_el(location,:);
        [x_u,y_u,z_u] = sph2cart(problem.phiUsers/360*2*pi,...
                                 problem.thetaUsers/360*2*pi,...
                                 3+50*ones(size(problem.phiUsers)));
        for i = 1:nUsers
            hF(location) = line([0,x_u(i)],[0,y_u(i)],[0,z_u(i)],'color',colorList{location},'LineWidth',2);
            if i == IDmax;   scatter3(x_u(i),y_u(i),z_u(i),'MarkerEdgeColor',colorList{location},'MarkerFaceColor',colorList{location});
            else;            scatter3(x_u(i),y_u(i),z_u(i),'MarkerEdgeColor',colorList{location},'MarkerFaceColor',colorList{location});
            end
%             text(x_u(i)+dx, y_u(i)+dy, z_u(i)+dz, strcat('User [',num2str(problem.phiUsers(i)),{' '},num2str(problem.thetaUsers(i)),']'));
        end
        legendList(location) = strcat('Use-case',{' '},mat2str(location));
    end
    % Draw array reference
    problem = o_read_input_problem('data/metaproblem_test.dat');
    problem.NxPatch = floor(sqrt(problem.N_Antennas));
    problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
    problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
    problem = o_create_subarray_partition(problem);
    problem.NzPatch = problem.NxPatch;
    problem.dz = problem.dx;
    problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                         [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
                         'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque s\ED
    handle_ConformalArray = phased.URA([problem.NyPatch,problem.NzPatch],...
                         'Lattice','Rectangular','Element',problem.handle_Ant,...
                         'ElementSpacing',[problem.dy,problem.dz]);
    possible_locations = handle_ConformalArray.getElementPosition;
    possible_locations = possible_locations.*1000;
    for i = 1:size(possible_locations,2)
        scatter3(possible_locations(1,i),possible_locations(2,i),possible_locations(3,i),'s','MarkerEdgeColor','k','MarkerFaceColor','k');
    end
    % Draw edges references
    a = [max(xlim)-5 0 0]; b = [0 min(ylim)+5 0]; c = [0 0 max(zlim)];
    starts = zeros(3,3);
    ends = [a;b;c];
    quiver3(starts(:,1), starts(:,2), starts(:,3), ends(:,1), ends(:,2), ends(:,3),'color','k','lineStyle','-.','LineWidth',1.5);
    hl = legend(hF,legendList);
    set(hl,'FontSize',9,'Location','SouthWest');
    view([50,30]);
    title('User use-cases fixed location','FontSize',12);
    xlabel('x-plane','FontSize',12);
    ylabel('y-plane','FontSize',12);
    zlabel('z-plane','FontSize',12);
end







%% EXPERIMENT 3
function experiment3(input)
%     fprintf('Running experiment 1...\n');
%     fprintf('Input parameters:\n\tnUsers:\t%d\n\tIter:\t%d\n\tAntennas:\t%d\n\tarrayRestriction:\t%s\n\talgorithm:\t%s\n\tdetLocation:\t%s\n\tuseCasesLocation:\t%s\n\tuseCaseLocation:\t%d\n',input.nUsers,input.nIter,input.nAntennas,input.arrayRestriction,input.algorithm,mat2str(input.detLocation),mat2str(input.useCasesLocation),input.useCaseLocation);
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override parameters (config)
    conf.verbosity            = 0;
    conf.algorithm            = input.algorithm;
	conf.detLocation          = input.detLocation;
    conf.useCasesLocation     = input.useCasesLocation;
    conf.useCaseLocation      = input.useCaseLocation;
    % Override parameters (problem)
    problem.nUsers            = input.nUsers;
    problem.N_Antennas        = input.nAntennas;
    problem.payload           = input.payload;
    problem.deadline          = input.deadline;
    problem.arrayRestriction  = input.arrayRestriction;
    problem.DEBUG             = true;  % no ouput messages (Scheduling)
    problem.PLOT_DEBUG        = false;  % No plotting
    problem.numPkts           = 2;
    % Output variables
    ratioOKTot = zeros(input.nUsers,input.nIter);
    ratioNOKTot = zeros(input.nUsers,input.nIter);
    for iter = 1:input.nIter
        % Configure the simulation environment
        [problem,~,flows] = f_configuration(conf,problem);
        % Main function
        [flows,~,~,~,~,~] = main(conf,problem,flows);
        % Report of single execution
        [ratioOKTot(:,iter),ratioNOKTot(:,iter)] = f_generateReport(flows,problem.DEBUG);
        fprintf('\tNusers %d -> OK=%.2f(%%)\n',input.nUsers,mean(ratioOKTot(:,iter)));
    end
    % Store results in temp file (to be retrieved by plotting function)
    fileName = strcat('temp/exp3_',conf.algorithm,'_',problem.arrayRestriction,'_',mat2str(problem.nUsers),'_',mat2str(problem.deadline),'_',mat2str(problem.N_Antennas),'_',mat2str(conf.detLocation),'_',mat2str(conf.useCasesLocation),'_',mat2str(conf.useCaseLocation));
    % Store in local to store in file
    ratioOK          = mean(ratioOKTot,2);  %#ok
    ratioNOK         = mean(ratioNOKTot,2);  %#ok
    nUsers           = problem.nUsers;  %#ok
    nAntennas        = problem.N_Antennas;  %#ok
    deadline         = problem.deadline;  %#ok
    arrayRestriction = problem.arrayRestriction;  %#ok
    detLocation      = conf.detLocation;  %#ok
    useCaseLocation  = conf.useCaseLocation;  %#ok
    useCaseLocation  = conf.useCaseLocation;  %#ok
    save(fileName,'ratioOK','ratioNOK',...
         'nUsers','deadline','arrayRestriction',...
         'detLocation','useCaseLocation','useCaseLocation');
end



%% EXPERIMENT 3 - PLOTTING
function experiment3_plot(input)
    h =  findobj('type','figure');
    figIdx = length(h) + 1;
    ratioOK_av = zeros(length(input.nUsersList),1);
    ratioNOK_av = zeros(length(input.nUsersList),1);
    leg = cell(length(input.deadlineList),1);
    LineStyleList = {'-.','-','--',':'};
    MarkerList = {'s','*','d','+'};
    for deadline = input.deadlineList
        for nUsers = input.nUsersList
            fileName = strcat('temp/exp3_',input.algorithm,'_',...
                                           input.arrayRestriction,'_',...
                                           mat2str(nUsers),'_',...
                                           mat2str(deadline),'_',...
                                           mat2str(input.nAntennas),'_',...
                                           mat2str(input.detLocation),'_',...
                                           mat2str(input.useCasesLocation),'_',...
                                           mat2str(input.useCaseLocation));
            load(fileName,'ratioOK','ratioNOK');
            ratioOK_av(nUsers==input.nUsersList) = mean(ratioOK);
            ratioNOK_av(nUsers==input.nUsersList) = mean(ratioNOK);
        end
        idxMarkerStyle= mod(find(deadline==input.deadlineList),4)+1;
        idxLineStyle = ceil(find(deadline==input.deadlineList)/4);
        % No interpolation
        figure(figIdx); hold on; grid minor;
        plot(input.nUsersList,ratioOK_av,'LineWidth',1,'LineStyle',LineStyleList{idxLineStyle},'Marker',MarkerList{idxMarkerStyle},'Color','k');
        % With interpolation
        figure(figIdx+1); hold on; grid minor;
        xq = (1:1:input.nUsersList(end));
        ratioOK_av_int = interp1(input.nUsersList,ratioOK_av,xq,'pchip');
        plot(xq,ratioOK_av_int,'LineWidth',1,'LineStyle',LineStyleList{idxLineStyle},'Marker',MarkerList{idxMarkerStyle},'Color','k');
        % Legend
        leg(deadline==input.deadlineList) = strcat('Network sat.',{' '},sprintf('%.2f',deadline/max(input.deadlineList)),'x');
    end
    figure(figIdx); grid minor;
    xlabel('Number of users','FontSize',12);
    ylabel('Ratio OK (%)','FontSize',12);
    title('Ratio of application that meet deadline (%)','FontSize',12);
    hleg = legend(leg);
    set(hleg,'FontSize',10);
    figure(figIdx+1); grid minor;
    xlabel('Number of users','FontSize',12);
    ylabel('Ratio OK (%)','FontSize',12);
    title('Ratio of application that meet deadline (%)','FontSize',12);
    hleg = legend(leg);
    set(hleg,'FontSize',10);
end






% EXPERIMENT 4 - Convergency analysis
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
    conf.NumPhaseShifterBits = 2;  % Number of bits to control heuristic solution
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

            %% Create subarray partition
            problem = o_create_subarray_partition(problem);

            problem.NzPatch = problem.NxPatch;
            problem.dz = problem.dx;

            %% Create the antenna handler and the data structure with all possible pos.
            problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
                'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque si
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
    conf.verbosity = 0;
    conf.algorithm = algorithm;  % Heuristic algorithm
    conf.PopulationSize_Data = 30;
    conf.Maxgenerations_Data = 150;
    conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);
    conf.MaxStallgenerations_Data = ceil(conf.Maxgenerations_Data/10);  % Force it to cover all the generations
    conf.FunctionTolerance_Data = 1e-10;  % Heuristics stops when not improving solution by this much
    conf.NumPhaseShifterBits = 60;  % Number of bits to control heuristic solution
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
            [~,W,arrayHandle,estObj] = f_heuristics(problem,conf,candSet);
            % Heuristics - Post Processing
            if conf.MinObjFIsSNR;     CapTot(:,idxAnt,idxIter)  = log2(estObj+1);  % in bps/Hz
                                      SINRTot(:,idxAnt,idxIter) = pow2db(estObj);  % in dB
            else;                     CapTot(:,idxAnt,idxIter)  = estObj;  % in bps/Hz
                                      SINRTot(:,idxAnt,idxIter) = pow2db(2.^(estTH/problem.Bw) - 1);  % in dB
            end
            % Process results after Beamforming
            [DirOKTot(:,idxAnt,idxIter),DirNOKTot(:,:,idxAnt,idxIter)]  = ...
                f_BF_results(W,arrayHandle,problem,conf,false);
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
function experiment5_plot(varargin)
    h =  findobj('type','figure');
    figIdx = length(h) + 1;
    
    if nargin ==1
        fileName = varargin{1};
        % EXPERIMENT 5 - Plotting results
        % Load results
        load(fileName,'DirOK','DirNOK_gntd','DirNOK_pcvd','nUsers','nAntennasList','arrayRestriction');
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
        
    elseif nargin == 4
        nUsers = varargin{1};
        arrRestctList = varargin{2};
        locList = varargin{3};
        my_nAntennasList = varargin{4};

        fileNameList = cell(length(arrRestctList),1);
        legendList = cell(length(arrRestctList),1);
        % colorList from: https://en.wikipedia.org/wiki/Web_colors
        colNoneList = [[139 0 0]./255 ; [178 34 34]./255 ; [220 20 60]./255 ; [205 92 92]./255 ; [250 128 114]./255 ; [255 160 122]./255];
        colLocList  = [[0 0 128]./255 ; [0 0 205]./255 ; [0 0 255]./255 ; [65 105 225]./255 ; [100 149 237]./255 ; [0 191 255]./255];
        colListTot = [colNoneList(1:length(locList),:);colLocList(1:length(locList),:)];
        MarkerList = {'s','*','o','x','d','p','h','^','v','>','<'};

        for locIdx = 1:length(locList)
            fileNameList{1} = strcat('temp/exp5_GA_',arrRestctList{1},'_2_true_true_',mat2str(locList(locIdx)));
            fileNameList{2} = strcat('temp/exp5_GA_',arrRestctList{2},'_2_true_true_',mat2str(locList(locIdx)));
            fileName = strcat('temp/exp5_GA_TOT_',mat2str(nUsers),'_true_true_',mat2str(locList(locIdx)));
            % Store global results
            load(fileName,'Cap_tot','SINR_BB_tot','SINR_PB_tot','nUsers','nAntennasList','arrRestctList');
            Cap_TOT(:,locIdx) = Cap_tot(:,1);  %#ok
            Cap_TOT(:,locIdx + length(locList)) = Cap_tot(:,2);  %#ok
            legendList(locIdx)                   = strcat(' No restr. - Loc. ',{' '},mat2str(locList(locIdx)));
            legendList(locIdx + length(locList)) = strcat('Sq. restr. - Loc. ',{' '},mat2str(locList(locIdx)));
            figure(figIdx); hold on;
            p(locIdx) = plot(my_nAntennasList,Cap_tot(:,2),'Color','k','LineStyle',':','Marker',MarkerList{locIdx});
            p(locIdx + length(locList)) = plot(my_nAntennasList,Cap_tot(:,1),'Color','k','LineStyle','-','Marker',MarkerList{locIdx});
        end
        title('Average capacity','FontSize',12);
        xlabel('Number of antennas','FontSize',12);
        ylabel('Capacity in bits/Hz/s','FontSize',12);
        ah1 = gca;
        leg = legend(ah1,p(1:length(locList)),legendList(length(locList)+1:2*length(locList)));
        set(leg,'FontSize',10);
        ah2 = axes('position',get(gca,'position'),'visible','off');
        leg = legend(ah2,p(length(locList)+1:2*length(locList)),legendList(1:length(locList)));
        set(leg,'FontSize',10);
        grid minor;
        
        figIdx = figIdx + 1;
        figure(figIdx); hold on; figIdx = figIdx + 1;
        set(gca, 'ColorOrder', colListTot, 'NextPlot', 'replacechildren');  % Change to new colors.
        [~,ia,~] = intersect(nAntennasList,my_nAntennasList);  %#ok
        nAntennasList = nAntennasList(ia);
        plot(nAntennasList,Cap_TOT(1:length(nAntennasList),:),'LineWidth',1.3,'Marker','x','MarkerSize',3.5);
        title('Average capacity','FontSize',12);
        xlabel('Number of antennas','FontSize',12);
        ylabel('Capacity in bits/Hz/s','FontSize',12);
        hList = findobj('Type', 'line');  % Returns lines in inverse order
        ah1 = gca;
        leg = legend(ah1,hList(1:length(locList)),legendList(length(locList)+1:2*length(locList)));
        set(leg,'FontSize',10);
        ah2 = axes('position',get(gca,'position'),'visible','off');
        leg = legend(ah2,hList(length(locList)+1:2*length(locList)),legendList(1:length(locList)));
        set(leg,'FontSize',10);
        figure(figIdx - 1);
        xlim([min(nAntennasList)-50 max(nAntennasList)+50]);
        grid minor;

        data(:,:,1) = Cap_TOT(1:length(nAntennasList),:);

        groupLabels = nAntennasList;
        figure(figIdx); hold on;
        xlim([0.5 length(nAntennasList)+0.5]);
        grid minor;
        set(gca, 'ColorOrder', colListTot, 'NextPlot', 'replacechildren');  % Change to new colors
        plotBarStackGroups(data, groupLabels);
        title('Average capacity','FontSize',14);
        xlabel('Number of antennas','FontSize',12);
        ylabel('Capacity in bits/Hz/s','FontSize',12);
        hList = findobj('Type', 'line');  % Returns lines in inverse order
        ah1 = gca;
        leg = legend(ah1,hList(1:length(locList)),legendList(length(locList)+1:2*length(locList)));
        set(leg,'FontSize',10,'Location','SouthEast');
        ah2 = axes('position',get(gca,'position'),'visible','off');
        leg = legend(ah2,hList(length(locList)+1:2*length(locList)),legendList(1:length(locList)));
        set(leg,'FontSize',10,'Location','SouthWest');
        figure(figIdx);
    else
        error('wrong number of parameters');
    end
end







%% EXPERIMENT 6
function experiment6(input,plotFLAG)
    %
    fprintf('Running experiment 6...\n');
	% Store input struct in local
    nUsers               = input.nUsers;
    nIter                = input.nIter;
    nAntennas            = input.nAntennas;
    arrayRestriction     = input.arrayRestriction;
    algorithm            = input.algorithm;
    detLocation          = input.detLocation;
    useCasesLocation     = input.useCasesLocation;
    useCaseLocationList  = input.useCaseLocationList;
    fprintf('Input parameters:\n\tnUsers:\t%d\n\tIter:\t%d\n\tAntennas:\t%d\n\tarrayRestriction:\t%s\n\talgorithm:\t%s\n\tdetLocation:\t%s\n\tuseCasesLocation:\t%s\n\tuseCaseLocation:\t%s\n',nUsers,nIter,nAntennas,arrayRestriction,algorithm,mat2str(detLocation),mat2str(useCasesLocation),mat2str(useCaseLocationList));
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override (problem) parameters
    problem.nUsers = nUsers;  % Number of users in the simulation
    problem.MinObjFIsSNR = true;  % (arbitrary)
    problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
    problem.arrayRestriction = arrayRestriction;
    % Override (conf) parameters
    conf.verbosity = 1;
    conf.algorithm = algorithm;  % Heuristic algorithm
    conf.NumPhaseShifterBits = 60;  % Number of bits to control heuristic solution
    conf.FunctionTolerance_Data = 1e-10;  % Heuristics stops when not improving solution by this much
    conf.multiPath = false;  % LoS channel (for now)
	% Configure basic parameters
    candSet = (1:1:problem.nUsers);  % Set of users to be considered
	% Create output variables
    CapHEUTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % Capacity (heuristics)
    SINRHEUTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % SINR (heuristics)
    CapHEU_CBFTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % Capacity (heuristics)
    SINRHEU_CBFTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % SINR (heuristics)
    CapHEU_LCMVTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % Capacity (heuristics)
    SINRHEU_LCMVTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % SINR (heuristics)
    CapLCMVTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % Capacity (LCMV)
    SINRLCMVTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % SINR (LCMV)
    CapCBFTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % Capacity (Conventional)
    SINRCBFTot = zeros(problem.nUsers,length(useCaseLocationList),nIter);  % SINR (Conventional)
    DirOKHEUTot = -Inf(problem.nUsers,length(useCaseLocationList),nIter);  % Directivity target (heuristics)
    DirNOKHEUTot = -Inf(problem.nUsers,problem.nUsers,length(useCaseLocationList),nIter);  % Directivity others (heuristics)
    DirOKHEU_LCMVTot = -Inf(problem.nUsers,length(useCaseLocationList),nIter);  % Directivity target (heuristics)
    DirNOKHEU_LCMVTot = -Inf(problem.nUsers,problem.nUsers,length(useCaseLocationList),nIter);  % Directivity others (heuristics)
    DirOKHEU_CBFTot = -Inf(problem.nUsers,length(useCaseLocationList),nIter);  % Directivity target (heuristics)
    DirNOKHEU_CBFTot = -Inf(problem.nUsers,problem.nUsers,length(useCaseLocationList),nIter);  % Directivity others (heuristics)
    DirOKLCMVTot = -Inf(problem.nUsers,length(useCaseLocationList),nIter);  % Directivity target (LCMV)
    DirNOKLCMVTot = -Inf(problem.nUsers,problem.nUsers,length(useCaseLocationList),nIter);  % Directivity others (LCMV)
    DirOKCBFTot = -Inf(problem.nUsers,length(useCaseLocationList),nIter);  % Directivity target (Conventional)
    DirNOKCBFTot = -Inf(problem.nUsers,problem.nUsers,length(useCaseLocationList),nIter);  % Directivity others (Conventional)
    % Linearize combinations and asign Population size (To be replaced with
    % convergency analysis values)
%     totComb = log10(problem.nUsers.*factorial(ceil(nAntennasList/problem.nUsers)));
%     maxPop = 70;  % Maximum population size
%     minPop = 40;  % Minimum population size
%     slope = (maxPop - minPop) / (totComb(end)-totComb(1));
%     ordIdx = minPop - slope*totComb(1);
%     PopSizeList = ceil(slope*totComb + ordIdx);
    PopSizeList = 30*ones(length(useCaseLocationList),1);
    % Main execution
    for idxLoc = 1:length(useCaseLocationList)
        conf.PopulationSize_Data = PopSizeList(idxLoc);
        conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);
%         conf.Maxgenerations_Data = 150;
%         conf.MaxStallgenerations_Data = ceil(conf.Maxgenerations_Data/4);
        conf.Maxgenerations_Data = 30;
        conf.MaxStallgenerations_Data = 10;
        % Select the localization
        conf.detLocation = detLocation;
        conf.useCasesLocation = useCasesLocation;
        conf.useCaseLocation = useCaseLocationList(idxLoc);
        fprintf('Nusers=%d - Nant=%d - rest=%s - Loc=%d\n',problem.nUsers,nAntennas,arrayRestriction,conf.useCaseLocation);
        for idxIter = 1:nIter
            fprintf('Iteration %d with PopSize %d\n',idxIter,PopSizeList(idxLoc));
            % Configure the simulation environment. Need to place users in new
            % locations and create new channels to have statistically
            % meaningful results
            [problem,~,~] = f_configuration(conf,problem);
            % Select number of antennas
            problem.N_Antennas = nAntennas;
            % Adjust parameters
            problem.NxPatch = floor(sqrt(problem.N_Antennas));
            problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
            problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
            % Conventional Beamforming + LCMV
            [W_LCMV,W_CBF,handle_ConformalArray,~,~] = f_conventionalBF(problem,conf,candSet);
            fprintf('SOLVING LCMV\n');
            % Compute directivity for LCMV
            [DirOKLCMVTot(:,idxLoc,idxIter),DirNOKLCMVTot(:,:,idxLoc,idxIter), ...
             CapLCMVTot(:,idxLoc,idxIter),SINRLCMVTot(:,idxLoc,idxIter)]     = ...
               f_BF_results(W_LCMV,handle_ConformalArray,candSet,problem,conf,plotFLAG);
            % Compute directivity for Conventional
            fprintf('SOLVING CBF\n');
            [DirOKCBFTot(:,idxLoc,idxIter),DirNOKCBFTot(:,:,idxLoc,idxIter), ...
             CapCBFTot(:,idxLoc,idxIter),SINRCBFTot(:,idxLoc,idxIter)]     = ...
               f_BF_results(W_CBF,handle_ConformalArray,candSet,problem,conf,plotFLAG);
            % Call heuristics - No bootstrapping
            fprintf('SOLVING HEURISTICS\n');
%             [~,W,handle_ConformalArray,~] = f_heuristics(problem,conf,candSet);
            [~,W,handle_ConformalArray,~] = CBG_solveit(problem,conf,candSet);
            % Compute directivity for Heuristics
            [DirOKHEUTot(:,idxLoc,idxIter),DirNOKHEUTot(:,:,idxLoc,idxIter), ...
             CapHEUTot(:,idxLoc,idxIter),SINRHEUTot(:,idxLoc,idxIter)]     = ...
               f_BF_results(W,handle_ConformalArray,candSet,problem,conf,plotFLAG);
           % Call heuristics - bootstrapping from LCMV
            fprintf('SOLVING HEURISTICS - BOOTSTRAPPING LCMV\n');
            problem.initialW = W_LCMV;
%             [~,W,handle_ConformalArray,~] = f_heuristics(problem,conf,candSet);
            [~,W,handle_ConformalArray,~] = CBG_solveit(problem,conf,candSet);
            % Compute directivity for Heuristics + bootstrapped LCMV
            [DirOKHEU_LCMVTot(:,idxLoc,idxIter),DirNOKHEU_LCMVTot(:,:,idxLoc,idxIter), ...
             CapHEU_LCMVTot(:,idxLoc,idxIter),SINRHEU_LCMVTot(:,idxLoc,idxIter)]     = ...
               f_BF_results(W,handle_ConformalArray,candSet,problem,conf,plotFLAG);
           % Call heuristics - bootstrapping from CBF
            fprintf('SOLVING HEURISTICS - BOOTSTRAPPING CBF\n');
            problem.initialW = W_CBF;
%             [~,W,handle_ConformalArray,~] = f_heuristics(problem,conf,candSet);
            [~,W,handle_ConformalArray,~,~,~] = CBG_solveit(problem,conf,candSet);
            % Compute directivity for Heuristics + bootstrapped CBF
            [DirOKHEU_CBFTot(:,idxLoc,idxIter),DirNOKHEU_CBFTot(:,:,idxLoc,idxIter), ...
             CapHEU_CBFTot(:,idxLoc,idxIter),SINRHEU_CBFTot(:,idxLoc,idxIter)]     = ...
               f_BF_results(W,handle_ConformalArray,candSet,problem,conf,plotFLAG);
        end
        
        % Parse results for specific case - Heuristics
        DirOKHEUTot_lin = db2pow(DirOKHEUTot);
        DirNOKHEUTot_lin = db2pow(DirNOKHEUTot);
        DirOKHEU_lin = mean(DirOKHEUTot_lin(:,idxLoc,:),3);
        DirNOKHEU_gntd_lin = sum(mean(DirNOKHEUTot_lin(:,:,idxLoc,:),4),1).'; % Generated interference 
        DirNOKHEU_pcvd_lin = sum(mean(DirNOKHEUTot_lin(:,:,idxLoc,:),4),2); % Perceived interference
        DirOKHEU = pow2db(DirOKHEU_lin);  %#ok % Directivity generated to intended user (heuristics)
        DirNOKHEU_gntd = pow2db(DirNOKHEU_gntd_lin);  %#ok  % Directivity being generated by intended user (heuristics)
        DirNOKHEU_pcvd = pow2db(DirNOKHEU_pcvd_lin);  %#ok  % Directivity inflicted to intended user (heuristics)
        % Parse results for specific case - Heuristics - Bootstrap from LCMV
        DirOKHEU_LCMVTot_lin = db2pow(DirOKHEU_LCMVTot);
        DirNOKHEU_LCMVTot_lin = db2pow(DirNOKHEU_LCMVTot);
        DirOKHEU_LCMV_lin = mean(DirOKHEU_LCMVTot_lin(:,idxLoc,:),3);
        DirNOKHEU_LCMV_gntd_lin = sum(mean(DirNOKHEU_LCMVTot_lin(:,:,idxLoc,:),4),1).'; % Generated interference 
        DirNOKHEU_LCMV_pcvd_lin = sum(mean(DirNOKHEU_LCMVTot_lin(:,:,idxLoc,:),4),2); % Perceived interference
        DirOKHEU_LCMV = pow2db(DirOKHEU_LCMV_lin);  %#ok % Directivity generated to intended user (heuristics)
        DirNOKHEU_LCMV_gntd = pow2db(DirNOKHEU_LCMV_gntd_lin);  %#ok  % Directivity being generated by intended user (heuristics)
        DirNOKHEU_LCMV_pcvd = pow2db(DirNOKHEU_LCMV_pcvd_lin);  %#ok  % Directivity inflicted to intended user (heuristics)
        % Parse results for specific case - Heuristics - Bootstrap from CBF
        DirOKHEU_CBFTot_lin = db2pow(DirOKHEU_CBFTot);
        DirNOKHEU_CBFTot_lin = db2pow(DirNOKHEU_CBFTot);
        DirOKHEU_CBF_lin = mean(DirOKHEU_CBFTot_lin(:,idxLoc,:),3);
        DirNOKHEU_CBF_gntd_lin = sum(mean(DirNOKHEU_CBFTot_lin(:,:,idxLoc,:),4),1).'; % Generated interference 
        DirNOKHEU_CBF_pcvd_lin = sum(mean(DirNOKHEU_CBFTot_lin(:,:,idxLoc,:),4),2); % Perceived interference
        DirOKHEU_CBF = pow2db(DirOKHEU_CBF_lin);  %#ok % Directivity generated to intended user (heuristics)
        DirNOKHEU_CBF_gntd = pow2db(DirNOKHEU_CBF_gntd_lin);  %#ok  % Directivity being generated by intended user (heuristics)
        DirNOKHEU_CBF_pcvd = pow2db(DirNOKHEU_CBF_pcvd_lin);  %#ok  % Directivity inflicted to intended user (heuristics)
        % Parse results for specific case - LCMV
        DirOKLCMVTot_lin = db2pow(DirOKLCMVTot);
        DirNOKLCMVTot_lin = db2pow(DirNOKLCMVTot);
        DirOKLCMV_lin = mean(DirOKLCMVTot_lin(:,idxLoc,:),3);
        DirNOKLCMV_gntd_lin = sum(mean(DirNOKLCMVTot_lin(:,:,idxLoc,:),4),1).'; % Generated interference 
        DirNOKLCMV_pcvd_lin = sum(mean(DirNOKLCMVTot_lin(:,:,idxLoc,:),4),2); % Perceived interference
        DirOKLCMV = pow2db(DirOKLCMV_lin);  %#ok  % Directivity generated to intended user (LCMV)
        DirNOKLCMV_gntd = pow2db(DirNOKLCMV_gntd_lin);  %#ok  % Directivity being generated by intended user (LCMV)
        DirNOKLCMV_pcvd = pow2db(DirNOKLCMV_pcvd_lin);  %#ok  % Directivity inflicted to intended user (LCMV)
        % Parse results for specific case - CBF (Conventional)
        DirOKCBFTot_lin = db2pow(DirOKCBFTot);
        DirNOKCBFTot_lin = db2pow(DirNOKCBFTot);
        DirOKCBF_lin = mean(DirOKCBFTot_lin(:,idxLoc,:),3);
        DirNOKCBF_gntd_lin = sum(mean(DirNOKCBFTot_lin(:,:,idxLoc,:),4),1).'; % Generated interference 
        DirNOKCBF_pcvd_lin = sum(mean(DirNOKCBFTot_lin(:,:,idxLoc,:),4),2); % Perceived interference
        DirOKCBF = pow2db(DirOKCBF_lin);  %#ok  % Directivity generated to intended user (Conventional)
        DirNOKCBF_gntd = pow2db(DirNOKCBF_gntd_lin);  %#ok  % Directivity being generated by intended user (Conventional)
        DirNOKCBF_pcvd = pow2db(DirNOKCBF_pcvd_lin);  %#ok  % Directivity inflicted to intended user (Conventional)
        % Compute SINR and Capacities - Heuristics
        SINRHEU_PB_lin = mean(SINRHEUTot(:,idxLoc,:),3);  % Compute SINR Pass-Band (Linear)
        SINRHEU_PB = pow2db(SINRHEU_PB_lin);  %#ok  % Compute SINR Pass-Band (dB)
        CapHEU = log2(1 + SINRHEU_PB_lin);  % Compute final Average Capacity (bits/Hz/s)
        CapSumHEU = sum(CapHEU);  %#ok  % Compute final Total Capacity (bits/Hz/s)
        % Compute SINR and Capacities - Heuristics + Bootstrapping LCMV
        SINRHEU_LCMV_PB_lin = mean(SINRHEU_LCMVTot(:,idxLoc,:),3);  % Compute SINR Pass-Band (Linear)
        SINRHEU_LCMV_PB = pow2db(SINRHEU_LCMV_PB_lin);  %#ok  % Compute SINR Pass-Band (dB)
        CapHEU_LCMV = log2(1 + SINRHEU_LCMV_PB_lin);  % Compute final Average Capacity (bits/Hz/s)
        CapSumHEU_LCMV = sum(CapHEU_LCMV);  %#ok  % Compute final Total Capacity (bits/Hz/s)
        % Compute SINR and Capacities - Heuristics + Bootstrapping CBF
        SINRHEU_CBF_PB_lin = mean(SINRHEU_CBFTot(:,idxLoc,:),3);  % Compute SINR Pass-Band (Linear)
        SINRHEU_CBF_PB = pow2db(SINRHEU_CBF_PB_lin);  %#ok  % Compute SINR Pass-Band (dB)
        CapHEU_CBF = log2(1 + SINRHEU_CBF_PB_lin);  % Compute final Average Capacity (bits/Hz/s)
        CapSumHEU_CBF = sum(CapHEU_CBF);  %#ok  % Compute final Total Capacity (bits/Hz/s)
        % Compute SINR and Capacities - LCMV
        SINRLCMV_PB_lin = mean(SINRLCMVTot(:,idxLoc,:),3);  % Compute SINR Pass-Band (Linear)
        SINRLCMV_PB = pow2db(SINRLCMV_PB_lin);  %#ok  % Compute SINR Pass-Band (dB)
        CapLCMV = log2(1 + SINRLCMV_PB_lin);  % Compute final Average Capacity (bits/Hz/s)
        CapLCMVSum = sum(CapLCMV);  %#ok  % Compute final Total Capacity (bits/Hz/s)
        % Compute SINR and Capacities - LCMV
        SINRCBF_PB_lin = mean(SINRCBFTot(:,idxLoc,:),3);  % Compute SINR Pass-Band (Linear)
        SINRCBF_PB = pow2db(SINRCBF_PB_lin);  %#ok  % Compute SINR Pass-Band (dB)
        CapCBF = log2(1 + SINRCBF_PB_lin);  % Compute final Average Capacity (bits/Hz/s)
        CapCBFSum = sum(CapCBF);  %#ok  % Compute final Total Capacity (bits/Hz/s)
        % Store results in mat file
        fileName = strcat('temp/exp6_',conf.algorithm,'_',problem.arrayRestriction,'_',mat2str(nUsers),'_',mat2str(nAntennas),'_',mat2str(conf.detLocation),'_',mat2str(conf.useCasesLocation),'_',mat2str(conf.useCaseLocation));
        save(fileName,'DirOKHEU','DirNOKHEU_gntd','DirNOKHEU_pcvd','SINRHEU_PB','CapHEU','CapSumHEU',...
                      'DirOKHEU_LCMV','DirNOKHEU_LCMV_gntd','DirNOKHEU_LCMV_pcvd','SINRHEU_LCMV_PB','CapHEU_LCMV','CapSumHEU_LCMV',...
                      'DirOKHEU_CBF','DirNOKHEU_CBF_gntd','DirNOKHEU_CBF_pcvd','SINRHEU_CBF_PB','CapHEU_CBF','CapSumHEU_CBF',...
                      'DirOKLCMV','DirNOKLCMV_gntd','DirNOKLCMV_pcvd','SINRLCMV_PB','CapLCMV','CapLCMVSum',...
                      'DirOKCBF','DirNOKCBF_gntd','DirNOKCBF_pcvd','SINRCBF_PB','CapCBF','CapCBFSum',...
                      'nUsers','nAntennas','arrayRestriction');
    end
end

%% EXPERIMENT 6 - PLOTTING
% function experiment6_plot(nUsersList,arrRestct,locList,nAntennas)
function experiment6_plot(nUsersList,input)
    arrRestct = input.arrRestct;
    algorithm = input.algorithm;
    nAntennas = input.nAntennas;
    detLocation = input.detLocation;
    useCasesLocation = input.useCasesLocation;
    useCaseLocationList = input.useCaseLocationList;
    MarkerList = {'s','*','o','x','d','p','h','^','v','>','<'};
    for idxLoc = 1:length(useCaseLocationList)
        CapHEU_final = zeros(length(nUsersList),1);
        CapSumHEU_final = zeros(length(nUsersList),1);
        SINRHEU_PB_final = zeros(length(nUsersList),1);
        CapHEU_LCMV_final = zeros(length(nUsersList),1);
        CapSumHEU_LCMV_final = zeros(length(nUsersList),1);
        SINRHEU_LCMV_PB_final = zeros(length(nUsersList),1);
        CapHEU_CBF_final = zeros(length(nUsersList),1);
        CapSumHEU_CBF_final = zeros(length(nUsersList),1);
        SINRHEU_CBF_PB_final = zeros(length(nUsersList),1);
        CapLCMV_final = zeros(length(nUsersList),1);
        CapSumLCMV_final = zeros(length(nUsersList),1);
        SINRLCMV_PB_final = zeros(length(nUsersList),1);
        CapCBF_final = zeros(length(nUsersList),1);
        CapSumCBF_final = zeros(length(nUsersList),1);
        SINRCBF_PB_final = zeros(length(nUsersList),1);
        for idxUsers = 1:length(nUsersList)
            nUsers = nUsersList(idxUsers);
            loc = useCaseLocationList(idxLoc);
            fileName = strcat('temp/exp6_',algorithm,'_',arrRestct,'_',mat2str(nUsers),'_',mat2str(nAntennas),'_',mat2str(detLocation),'_',...
                mat2str(useCasesLocation),'_',mat2str(loc));
            load(fileName,'SINRHEU_PB','CapHEU','CapSumHEU',...
                          'SINRHEU_LCMV_PB','CapHEU_LCMV','CapSumHEU_LCMV',...
                          'SINRHEU_CBF_PB','CapHEU_CBF','CapSumHEU_CBF',...
                          'SINRLCMV_PB','CapLCMV','CapLCMVSum',...
                          'SINRCBF_PB','CapCBF','CapCBFSum',...
                          'nUsers','nAntennas','arrayRestriction');
            CapHEU_final(idxUsers) = pow2db(mean(db2pow(CapHEU)));
            CapSumHEU_final(idxUsers) = pow2db(mean(db2pow(CapSumHEU)));
            SINRHEU_PB_final(idxUsers) = pow2db(mean(db2pow(SINRHEU_PB)));
            CapHEU_LCMV_final(idxUsers) = pow2db(mean(db2pow(CapHEU_LCMV)));
            CapSumHEU_LCMV_final(idxUsers) = pow2db(mean(db2pow(CapSumHEU_LCMV)));
            SINRHEU_LCMV_PB_final(idxUsers) = pow2db(mean(db2pow(SINRHEU_LCMV_PB)));
            CapHEU_CBF_final(idxUsers) = pow2db(mean(db2pow(CapHEU_CBF)));
            CapSumHEU_CBF_final(idxUsers) = pow2db(mean(db2pow(CapSumHEU_CBF)));
            SINRHEU_CBF_PB_final(idxUsers) = pow2db(mean(db2pow(SINRHEU_CBF_PB)));
            CapLCMV_final(idxUsers) = pow2db(mean(db2pow(CapLCMV)));
            CapSumLCMV_final(idxUsers) = pow2db(mean(db2pow(CapLCMVSum)));
            SINRLCMV_PB_final(idxUsers) = pow2db(mean(db2pow(SINRLCMV_PB)));
            CapCBF_final(idxUsers) = pow2db(mean(db2pow(CapCBF)));
            CapSumCBF_final(idxUsers) = pow2db(mean(db2pow(CapCBFSum)));
            SINRCBF_PB_final(idxUsers) = pow2db(mean(db2pow(SINRCBF_PB)));
        end
        figure; hold on; grid minor;
        plot(nUsersList,CapHEU_final,'Color','k','Marker',MarkerList{4},'LineStyle','-');
        plot(nUsersList,CapHEU_LCMV_final,'Color','k','Marker',MarkerList{4},'LineStyle',':');
        plot(nUsersList,CapHEU_CBF_final,'Color','k','Marker',MarkerList{4},'LineStyle','--');
        plot(nUsersList,CapLCMV_final,'Color','k','Marker',MarkerList{2},'LineStyle','-.');
%         plot(nUsersList,CapCBF_final,'Color','k','Marker',MarkerList{3},'LineStyle','-.');
        xlabel('Number of users','FontSize',12);
        ylabel('Capacity (bits/Hz/s)','FontSize',12);
        title('Average Capacity achieved','FontSize',12);
        hleg = legend('Heuristics','Heuristics + LCMV','Heuristics + CBF','LCMV','CBF');
        set(hleg,'FontSize',10,'Location','NorthEast');
        figure; hold on; grid minor;
        plot(nUsersList,CapSumHEU_final,'Color','k','Marker',MarkerList{4},'LineStyle','-');
        plot(nUsersList,CapSumHEU_LCMV_final,'Color','k','Marker',MarkerList{4},'LineStyle',':');
        plot(nUsersList,CapSumHEU_CBF_final,'Color','k','Marker',MarkerList{4},'LineStyle','--');
        plot(nUsersList,CapSumLCMV_final,'Color','k','Marker',MarkerList{2},'LineStyle','-.');
%         plot(nUsersList,CapSumCBF_final,'Color','k','Marker',MarkerList{3},'LineStyle','-.');
        xlabel('Number of users','FontSize',12);
        ylabel('Capacity (bits/Hz/s)','FontSize',12);
        title('Total Capacity achieved','FontSize',12);
        hleg = legend('Heuristics','Heuristics + LCMV','Heuristics + CBF','LCMV','CBF');
        set(hleg,'FontSize',10,'Location','NorthEast');
        figure; hold on; grid minor;
        plot(nUsersList,SINRHEU_PB_final,'Color','k','Marker',MarkerList{4},'LineStyle','-');
        plot(nUsersList,SINRHEU_LCMV_PB_final,'Color','k','Marker',MarkerList{4},'LineStyle',':');
        plot(nUsersList,SINRHEU_CBF_PB_final,'Color','k','Marker',MarkerList{4},'LineStyle','--');
        plot(nUsersList,SINRLCMV_PB_final,'Color','k','Marker',MarkerList{2},'LineStyle','-.');
%         plot(nUsersList,SINRCBF_PB_final,'Color','k','Marker',MarkerList{3},'LineStyle','-.');
        xlabel('Number of users','FontSize',12);
        ylabel('SINR (dB)','FontSize',12);
        title('Average SINR achieved','FontSize',12);
        hleg = legend('Heuristics','Heuristics + LCMV','Heuristics + CBF','LCMV','CBF');
        set(hleg,'FontSize',10,'Location','NorthEast');
    end
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
    conf.NumPhaseShifterBits = 0;  % Number of bits to control heuristic solution
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
    conf.NumPhaseShifterBits = 0;  % Number of bits to control heuristic solution
    conf.NbitsAmplitude = 0;
    conf.FunctionTolerance_Data = 1e-6;  % Heuristics stops when not improving solution by this much
    conf.multiPath = false;  % LoS channel (for now)
    
    % Override GA parameters
    conf.PopulationSize_Data = 150;
    conf.Maxgenerations_Data = 100;
    conf.EliteCount_Data = 25;
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

%% EXPERIMENT 9
function experiment9(input,plotFLAG)
    %
    fprintf('Running experiment 9...\n');
	% Store input struct in local
    nUsers               = input.nUsers;
    nIter                = input.nIter;
    nAntennasList        = input.nAntennasList;
    arrayRestriction     = input.arrayRestriction;
    algorithm            = input.algorithm;
    detLocation          = input.detLocation;
    useCasesLocation     = input.useCasesLocation;
    useCaseLocation      = input.useCaseLocation;
    fprintf('Input parameters:\n\tnUsers:\t%d\n\tIter:\t%d\n\tAntennasList:\t%s\n\tarrayRestriction:\t%s\n\talgorithm:\t%s\n\tdetLocation:\t%s\n\tuseCasesLocation:\t%s\n\tuseCaseLocation:\t%s\n',nUsers,nIter,mat2str(nAntennasList),arrayRestriction,algorithm,mat2str(detLocation),mat2str(useCasesLocation),mat2str(useCaseLocation));
    % Load basic parameters
    problem = o_read_input_problem('data/metaproblem_test.dat');
    conf = o_read_config('data/config_test.dat');
    % Override (problem) parameters
    problem.nUsers = nUsers;  % Number of users in the simulation
    problem.MinObjFIsSNR = true;  % (arbitrary)
    problem.MinObjF = 100.*ones(1,problem.nUsers);  % Same #ant per user. Random SNR (30dB)
    problem.arrayRestriction = arrayRestriction;
    % Override (conf) parameters
    conf.verbosity = 1;
    conf.algorithm = algorithm;  % Heuristic algorithm
    conf.NumPhaseShifterBits = 60;  % Number of bits to control heuristic solution
    conf.FunctionTolerance_Data = 1e-10;  % Heuristics stops when not improving solution by this much
    conf.multiPath = false;  % LoS channel (for now)
    conf.detLocation = detLocation;
    conf.useCasesLocation = useCasesLocation;
    conf.useCaseLocation  = useCaseLocation;
    % Configure basic parameters
    candSet = (1:1:problem.nUsers);  % Set of users to be considered
	% Create output variables
    CapLCMVTot = zeros(problem.nUsers,length(nAntennasList),nIter);  % Capacity (LCMV)
    SINRLCMVTot = zeros(problem.nUsers,length(nAntennasList),nIter);  % SINR (LCMV)
    DirOKLCMVTot = -Inf(problem.nUsers,length(nAntennasList),nIter);  % Directivity target (LCMV)
    DirNOKLCMVTot = -Inf(problem.nUsers,problem.nUsers,length(nAntennasList),nIter);  % Directivity others (LCMV)
    CapCBFTot = zeros(problem.nUsers,length(nAntennasList),nIter);  % Capacity (Conventional)
    SINRCBFTot = zeros(problem.nUsers,length(nAntennasList),nIter);  % SINR (Conventional)
    DirOKCBFTot = -Inf(problem.nUsers,length(nAntennasList),nIter);  % Directivity target (Conventional)
    DirNOKCBFTot = -Inf(problem.nUsers,problem.nUsers,length(nAntennasList),nIter);  % Directivity others (Conventional)
    CapHEUTot = zeros(problem.nUsers,length(nAntennasList),nIter);  % Capacity (Conventional)
    SINRHEUTot = zeros(problem.nUsers,length(nAntennasList),nIter);  % SINR (Conventional)
    DirOKHEUTot = -Inf(problem.nUsers,length(nAntennasList),nIter);  % Directivity target (Conventional)
    DirNOKHEUTot = -Inf(problem.nUsers,problem.nUsers,length(nAntennasList),nIter);  % Directivity others (Conventional)
    % Linearize combinations and asign Population size (To be replaced with
    % convergency analysis values)
%     totComb = log10(problem.nUsers.*factorial(ceil(nAntennasList/problem.nUsers)));
%     maxPop = 70;  % Maximum population size
%     minPop = 40;  % Minimum population size
%     slope = (maxPop - minPop) / (totComb(end)-totComb(1));
%     ordIdx = minPop - slope*totComb(1);
%     PopSizeList = ceil(slope*totComb + ordIdx);
    PopSizeList = 30*ones(length(nAntennasList),1);
    % Main execution
    for idxAnt = 1:length(nAntennasList)
        conf.PopulationSize_Data = PopSizeList(idxAnt);
        conf.Maxgenerations_Data = 10;
        conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);
        conf.MaxStallgenerations_Data = ceil(conf.Maxgenerations_Data/4);
        fprintf('Nusers=%d - Nant=%d - rest=%s - Loc=%d\n',problem.nUsers,nAntennasList(idxAnt),arrayRestriction,conf.useCaseLocation);
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
            % Conventional Beamforming + LCMV
            [W_LCMV,W_CBF,handle_ConformalArray,~,~] = f_conventionalBF(problem,conf,candSet);
            % Heuristics
            problem.initialW = W_LCMV;
            [~,W_HEU,~,~,~,~] = CBG_solveit(problem,conf,candSet);
            % Normalize weights
            for id = 1:problem.nUsers; W_HEU(id,:) = (1/sqrt(W_HEU(id,:)*W_HEU(id,:)'))*W_HEU(id,:); end
            fprintf('SOLVING LCMV\n');
            % Compute directivity for LCMV
            [DirOKLCMVTot(:,idxAnt,idxIter),DirNOKLCMVTot(:,:,idxAnt,idxIter), ...
             CapLCMVTot(:,idxAnt,idxIter),SINRLCMVTot(:,idxAnt,idxIter)]     = ...
               f_BF_results(W_LCMV,handle_ConformalArray,candSet,problem,conf,plotFLAG);
            % Compute directivity for Conventional
            fprintf('SOLVING CBF\n');
            [DirOKCBFTot(:,idxAnt,idxIter),DirNOKCBFTot(:,:,idxAnt,idxIter), ...
             CapCBFTot(:,idxAnt,idxIter),SINRCBFTot(:,idxAnt,idxIter)]     = ...
               f_BF_results(W_CBF,handle_ConformalArray,candSet,problem,conf,plotFLAG);
           % Compute directivity for Conventional
            fprintf('SOLVING CBF\n');
            [DirOKHEUTot(:,idxAnt,idxIter),DirNOKHEUTot(:,:,idxAnt,idxIter), ...
             CapHEUTot(:,idxAnt,idxIter),SINRHEUTot(:,idxAnt,idxIter)]     = ...
               f_BF_results(W_HEU,handle_ConformalArray,candSet,problem,conf,plotFLAG);
        end

        % Parse results for specific case - LCMV
        DirOKLCMVTot_lin = db2pow(DirOKLCMVTot);
        DirNOKLCMVTot_lin = db2pow(DirNOKLCMVTot);
        DirOKLCMV_lin = mean(DirOKLCMVTot_lin(:,idxAnt,:),3);
        DirNOKLCMV_gntd_lin = sum(mean(DirNOKLCMVTot_lin(:,:,idxAnt,:),4),1).'; % Generated interference 
        DirNOKLCMV_pcvd_lin = sum(mean(DirNOKLCMVTot_lin(:,:,idxAnt,:),4),2); % Perceived interference
        DirOKLCMV = pow2db(DirOKLCMV_lin);  %#ok  % Directivity generated to intended user (LCMV)
        DirNOKLCMV_gntd = pow2db(DirNOKLCMV_gntd_lin);  %#ok  % Directivity being generated by intended user (LCMV)
        DirNOKLCMV_pcvd = pow2db(DirNOKLCMV_pcvd_lin);  %#ok  % Directivity inflicted to intended user (LCMV)
        % Parse results for specific case - CBF (Conventional)
        DirOKCBFTot_lin = db2pow(DirOKCBFTot);
        DirNOKCBFTot_lin = db2pow(DirNOKCBFTot);
        DirOKCBF_lin = mean(DirOKCBFTot_lin(:,idxAnt,:),3);
        DirNOKCBF_gntd_lin = sum(mean(DirNOKCBFTot_lin(:,:,idxAnt,:),4),1).'; % Generated interference 
        DirNOKCBF_pcvd_lin = sum(mean(DirNOKCBFTot_lin(:,:,idxAnt,:),4),2); % Perceived interference
        DirOKCBF = pow2db(DirOKCBF_lin);  %#ok  % Directivity generated to intended user (Conventional)
        DirNOKCBF_gntd = pow2db(DirNOKCBF_gntd_lin);  %#ok  % Directivity being generated by intended user (Conventional)
        DirNOKCBF_pcvd = pow2db(DirNOKCBF_pcvd_lin);  %#ok  % Directivity inflicted to intended user (Conventional)
        % Parse results for specific case - HEU
        DirOKHEUTot_lin = db2pow(DirOKHEUTot);
        DirNOKHEUTot_lin = db2pow(DirNOKHEUTot);
        DirOKHEU_lin = mean(DirOKHEUTot_lin(:,idxAnt,:),3);
        DirNOKHEU_gntd_lin = sum(mean(DirNOKHEUTot_lin(:,:,idxAnt,:),4),1).'; % Generated interference 
        DirNOKHEU_pcvd_lin = sum(mean(DirNOKHEUTot_lin(:,:,idxAnt,:),4),2); % Perceived interference
        DirOKHEU = pow2db(DirOKHEU_lin);  %#ok  % Directivity generated to intended user (Conventional)
        DirNOKHEU_gntd = pow2db(DirNOKHEU_gntd_lin);  %#ok  % Directivity being generated by intended user (Conventional)
        DirNOKHEU_pcvd = pow2db(DirNOKHEU_pcvd_lin);  %#ok  % Directivity inflicted to intended user (Conventional)
        % Compute SINR and Capacities - LCMV
        SINRLCMV_PB_lin = mean(SINRLCMVTot(:,idxAnt,:),3);  % Compute SINR Pass-Band (Linear)
        SINRLCMV_PB = pow2db(SINRLCMV_PB_lin);  %#ok  % Compute SINR Pass-Band (dB)
        CapLCMV = log2(1 + SINRLCMV_PB_lin);  % Compute final Average Capacity (bits/Hz/s)
        CapLCMVSum = sum(CapLCMV);  %#ok  % Compute final Total Capacity (bits/Hz/s)
        % Compute SINR and Capacities - CBF (Conventional)
        SINRCBF_PB_lin = mean(SINRCBFTot(:,idxAnt,:),3);  % Compute SINR Pass-Band (Linear)
        SINRCBF_PB = pow2db(SINRCBF_PB_lin);  %#ok  % Compute SINR Pass-Band (dB)
        CapCBF = log2(1 + SINRCBF_PB_lin);  % Compute final Average Capacity (bits/Hz/s)
        CapCBFSum = sum(CapCBF);  %#ok  % Compute final Total Capacity (bits/Hz/s)
        % Compute SINR and Capacities - HEU
        SINRHEU_PB_lin = mean(SINRHEUTot(:,idxAnt,:),3);  % Compute SINR Pass-Band (Linear)
        SINRHEU_PB = pow2db(SINRHEU_PB_lin);  %#ok  % Compute SINR Pass-Band (dB)
        CapHEU = log2(1 + SINRHEU_PB_lin);  % Compute final Average Capacity (bits/Hz/s)
        CapHEUSum = sum(CapHEU);  %#ok  % Compute final Total Capacity (bits/Hz/s)
        % Store results in mat file
        fileName = strcat('temp/exp9_',conf.algorithm,'_',problem.arrayRestriction,'_',mat2str(nUsers),'_',mat2str(nAntennasList(idxAnt)),'_',mat2str(conf.detLocation),'_',mat2str(conf.useCasesLocation),'_',mat2str(conf.useCaseLocation));
        save(fileName,'DirOKLCMV','DirNOKLCMV_gntd','DirNOKLCMV_pcvd','SINRLCMV_PB','CapLCMV','CapLCMVSum',...
                      'DirOKCBF','DirNOKCBF_gntd','DirNOKCBF_pcvd','SINRCBF_PB','CapCBF','CapCBFSum',...
                      'DirOKHEU','DirNOKHEU_gntd','DirNOKHEU_pcvd','SINRHEU_PB','CapHEU','CapHEUSum',...
                      'nUsers','nAntennasList','arrayRestriction','useCaseLocation');
    end
end

%% EXPERIMENT 9 - PLOTTING
function experiment9_plot(nUsersList,input)
    arrRestct = input.arrRestct;
    algorithm = input.algorithm;
    nAntennasList = input.nAntennasList;
    detLocation = input.detLocation;
    useCasesLocation = input.useCasesLocation;
    useCaseLocation = input.useCaseLocation;
    MarkerList = {'s','*','o','x','d','p','h','^','v','>','<'};
    legList = cell(length(nAntennasList),1);
    CapLCMV_final = zeros(length(nUsersList),length(nAntennasList));
    CapSumLCMV_final = zeros(length(nUsersList),length(nAntennasList));
    SINRLCMV_PB_final = zeros(length(nUsersList),length(nAntennasList));
    CapCBF_final = zeros(length(nUsersList),length(nAntennasList));
    CapSumCBF_final = zeros(length(nUsersList),length(nAntennasList));
    SINRCBF_PB_final = zeros(length(nUsersList),length(nAntennasList));
    CapHEU_final = zeros(length(nUsersList),length(nAntennasList));
    CapSumHEU_final = zeros(length(nUsersList),length(nAntennasList));
    SINRHEU_PB_final = zeros(length(nUsersList),length(nAntennasList));
    for idxAnt = 1:length(nAntennasList)
        for idxUsers = 1:length(nUsersList)
            nUsers = nUsersList(idxUsers);
            nAntennas = nAntennasList(idxAnt);
            fileName = strcat('temp/exp9_',algorithm,'_',arrRestct,'_',mat2str(nUsers),'_',mat2str(nAntennas),'_',mat2str(detLocation),'_',...
                mat2str(useCasesLocation),'_',mat2str(useCaseLocation));
            load(fileName,'SINRLCMV_PB','CapLCMV','CapLCMVSum',...
                          'SINRCBF_PB','CapCBF','CapCBFSum',...
                          'SINRHEU_PB','CapHEU','CapHEUSum',...
                          'nUsers','nAntennas','arrayRestriction');
            CapLCMV_final(idxUsers,idxAnt) = pow2db(mean(db2pow(CapLCMV)));
            CapSumLCMV_final(idxUsers,idxAnt) = pow2db(mean(db2pow(CapLCMVSum)));
            SINRLCMV_PB_final(idxUsers,idxAnt) = pow2db(mean(db2pow(SINRLCMV_PB)));
            CapCBF_final(idxUsers,idxAnt) = pow2db(mean(db2pow(CapCBF)));
            CapSumCBF_final(idxUsers,idxAnt) = pow2db(mean(db2pow(CapCBFSum)));
            SINRCBF_PB_final(idxUsers,idxAnt) = pow2db(mean(db2pow(SINRCBF_PB)));
            CapHEU_final(idxUsers,idxAnt) = pow2db(mean(db2pow(CapHEU)));
            CapSumHEU_final(idxUsers,idxAnt) = pow2db(mean(db2pow(CapHEUSum)));
            SINRHEU_PB_final(idxUsers,idxAnt) = pow2db(mean(db2pow(SINRHEU_PB)));
        end
        figure(1); hold on; grid minor;
        p1(2*idxAnt - 1) = plot(nUsersList,CapHEU_final(:,idxAnt),'Color','k','Marker',MarkerList{idxAnt},'LineStyle','-');
        p1(2*idxAnt) = plot(nUsersList,CapLCMV_final(:,idxAnt),'Color','k','Marker',MarkerList{idxAnt},'LineStyle',':');
        figure(2); hold on; grid minor;
        p2(2*idxAnt - 1) = plot(nUsersList,CapSumHEU_final(:,idxAnt),'Color','k','Marker',MarkerList{idxAnt},'LineStyle','-');
        p2(2*idxAnt) = plot(nUsersList,CapSumLCMV_final(:,idxAnt),'Color','k','Marker',MarkerList{idxAnt},'LineStyle',':');
        figure(3); hold on; grid minor;
        p3(2*idxAnt - 1) = plot(nUsersList,SINRHEU_PB_final(:,idxAnt),'Color','k','Marker',MarkerList{idxAnt},'LineStyle','-');
        p3(2*idxAnt) = plot(nUsersList,SINRLCMV_PB_final(:,idxAnt),'Color','k','Marker',MarkerList{idxAnt},'LineStyle',':');
        legList(2*idxAnt - 1) = strcat('HEU -',{' '},mat2str(nAntennas));        
        legList(2*idxAnt) = strcat('LCMV -',{' '},mat2str(nAntennas));
        legList(idxAnt) = strcat(mat2str(nAntennas),{' '},'antennas');
    end
    list = 1:length(p1);
    list = list(1:2:length(p1));
    figure(1); grid minor;
    xlabel('Number of users','FontSize',12);
    ylabel('Capacity (bits/Hz/s)','FontSize',12);
    title('Average Capacity achieved','FontSize',12);
    hleg = legend(p1(list),legList);
    set(hleg,'FontSize',10,'Location','NorthEast');
    figure(2); grid minor;
    xlabel('Number of users','FontSize',12);
    ylabel('Capacity (bits/Hz/s)','FontSize',12);
    title('Total Capacity achieved','FontSize',12);
    hleg = legend(p2(list),legList);
    set(hleg,'FontSize',10,'Location','NorthEast');
    figure(3); grid minor;
    xlabel('Number of users','FontSize',12);
    ylabel('SINR (dB)','FontSize',12);
    title('Average SINR achieved','FontSize',12);
    hleg = legend(p3(list),legList);
    set(hleg,'FontSize',10,'Location','NorthEast');
end