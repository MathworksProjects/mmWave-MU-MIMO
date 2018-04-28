function problem = f_configuration(problem)
    % To be included in data/metaproblem_test.dat once fixed o_read_in*.m
    problem.appPresence=[0.3,0.3,0.05,0.3,0.05];  % Presence of the application types specified in 'appNameInteretList' across users 
    problem.appNameInteretList={'Youtube','Justin TV','Facebook','Web Browsing','Twitter'};  % Applications to be considered in the simulations
    problem.targetPER = 1e-2;  % The target PER specified linearly (i.e. 1e-2 is PER of 1%)
    problem.numPkts = 100;

    % Define Application types and parametrize them - 5G applications Times
    % are defined following the specification "TS 23.203 standard Policy
    % and Charging Control Architecture" by the 3GPP community
    problem.Nclasses = 1;  % Number of classes (applications running)
    deadlineApps = zeros(length(problem.appNameInteretList),1);
    for appType = problem.appNameInteretList
        nUserstoApps = strcmp(appType,problem.appNameInteretList);
        if strcmp(appType,'Youtube') || strcmp(appType,'Justin TV') || strcmp(appType,'Vimeo')
            % Application type is video
            deadlineApps(nUserstoApps) = 50;  % 50ms deadline time for app type
        elseif strcmp(appType,'Facebook')
            % Application type is Facebook
            deadlineApps(nUserstoApps) = 100;  % 100ms deadline time for app type
            % Application type is other
        elseif strcmp(appType,'Web Browsing') || strcmp(appType,'Twitter')
            deadlineApps(nUserstoApps) = 500;  % 500ms deadline time for app type
        end
    end
    
    % Determine application presency based on DataSet
    if ~problem.manuallyAssignApp && strcmp(problem.trafficType,'dataSet')
        if ~exist('traffic','var') || ~exist('appNameList','var') || ~exist('appColorList','var')
            load('data/TrafficDataSetUPC','traffic','appNameList','appColorList');
        end
        % Allocate app types to users based on the population
        % distribution
        problem.appPresence = zeros(1,length(problem.appNameInteretList));
        nflowsTot = 0;
        % Extract the total number of flows
        for appName = problem.appNameInteretList
            nflowsTot = nflowsTot + traffic{strcmp(appName,appNameList)}.numFlows;  %#ok
        end
        % Determine the presence of a certain application based on the
        % number of flows vs the total
        for appName = problem.appNameInteretList
            problem.appPresence(1,strcmp(appName,problem.appNameInteretList)) = traffic{strcmp(appName,appNameList)}.numFlows/nflowsTot;
        end
    else
        % Check that the distribution is consisten across applications
        % We assume that the varibale probel.appPresence has already been
        % initialized previously (supposidly on data/meta*.m)
        if sum(problem.appPresence)~=1
            error('check appPresence in config\n');
        elseif ~isfield(problem,'appPresence')
            error('** ERROR: "appPresence" not manually set **\n');
        end
    end
    
    % Assign traffic type to users
    if problem.DEBUG
        % Select one application type across users
%         appTypePerUser = [1 2 3 4].';
        appTypePerUser = [1 1 1 1].';
    else
        % Distribute application proportionaly to the presence selected
        % and configured in 'problem.appPresence'
        appTypePerUser = gendist(problem.appPresence,problem.nUsers,1);
    end
    
    % Allocate flows of traffic to users assuming only one application type
    % is running at each device
    if strcmp(problem.trafficType,'synthetic')
        % Assign traffic type to users
        for u = 1:problem.nUsers
            appType = appTypePerUser(u);
            problem.class(u).deadline = deadlineApps(appType);  % Deadline in ms
            problem.class(u).name = problem.appNameInteretList{appType};
            problem.class(u).numPkts = problem.numPkts;  % Number of packets for the class
%             problem.class(u).iat = randi([50 200]);
            problem.class(u).iat = 40;  % Deterministic (constant) inter-arrival times (iat)
            problem.class(u).appColor = [0.9290 0.6940 0.1250];
            problem.class(u).payload = 1500*8;  % Default payload in bits
        end
    elseif strcmp(problem.trafficType,'dataSet')
        if ~exist('traffic','var') || ~exist('appNameList','var') || ~exist('appColorList','var')
            load('data/TrafficDataSetUPC','traffic','appNameList','appColorList');
        end
        for u = 1:problem.nUsers
            appType = appTypePerUser(u);
            problem.class(u).deadline = deadlineApps(appType);  % Deadline in ms
            problem.class(u).name = problem.appNameInteretList{appType};
            problem.class(u).numPkts = problem.numPkts;  % Orientative number of packets for the class
            appIdx = strcmp(problem.appNameInteretList{appType},appNameList);
            problem.class(u).trafficFlows = traffic{appIdx};  % Store flow information (inter-packet arrival times and payload)
            problem.class(u).appColor = appColorList(appIdx,:);  %#ok . % Application color to have consistency in plots
        end
    else
        error('ERROR - Traffic type not valid. Please introduce "synthetic" or "dataSet"\n');
    end
    
    % Load the PER-MCS table in the problem struct
    load('data/MCSPERTable.mat','mcsTable');
    problem.MCSPER = mcsTable;
end