%% Sets the configuration parameters (filepath, debug boolean and others)
%%%%  and read the input file to obtain the initial parameters 
%%%%  used for the Algorithm 1.
%   After this script is executed, several variables are available
%    in a struct that can be used for further processing

problem.dir = 'C:\Users\Santi\tfm\Heuristics\BRKGA\LabHeuristics\';
problem.inputDataFile = 'data\newtest_nus2.10_nant36_nmax2.8.dat';
problem.lengthFweights = 3;
problem.Noise = -104.5; % In dBs!!!
% Light speed constant
problem.c = 299792458;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% AUTOMATIC CONFIG PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Arguments that will be read by the readDATInputData function
% We read first all those that are not related to the number of users
problem.inputArgList = { struct('name','NxPatch','type','int');
            struct('name','NyPatch','type','int');
            struct('name','NxSubarrays','type','int');
            struct('name','NySubarrays','type','int');
            struct('name','arrayRestriction','type','char');
            struct('name','invdAntennas','type','float');
            struct('name','freq','type','float');
            struct('name','gamma','type','float');
            struct('name','nUsers','type','int');
            struct('name','Fweights','type','floatArray','size',problem.lengthFweights);
            struct('name','maxnChannelPaths','type','int');
            struct('name','anglesInRadians','type','char')};
input = readDATInputData([problem.dir,problem.inputDataFile],problem.inputArgList);

% Now we merge conf and input structs
f = fieldnames(input);
for i = 1:length(f)
  problem.(f{i}) = input.(f{i});
end

% Now, we will read the arguments that depend on the number of users.

problem.inputArgList = {
            struct('name','thetaUsers','type','floatArray','size',problem.nUsers);
            struct('name','phiUsers','type','floatArray','size',problem.nUsers);
            struct('name','dUsers','type','floatArray','size',problem.nUsers);
            struct('name','MaxThr','type','floatArray','size',problem.nUsers);
            struct('name','MinThr','type','floatArray','size',problem.nUsers);
            struct('name','thetaChannels','type','floatMatrix','size',...
                    {{problem.nUsers,problem.maxnChannelPaths}});
            struct('name','phiChannels','type','floatMatrix','size',...
                    {{problem.nUsers,problem.maxnChannelPaths}});
            struct('name','alphaChannels','type','floatMatrix','size',...
                    {{problem.nUsers,problem.maxnChannelPaths}})};
input = readDATInputData([problem.dir,problem.inputDataFile],problem.inputArgList);

% Now we merge conf and input structs
f = fieldnames(input);
for i = 1:length(f)
  problem.(f{i}) = input.(f{i});
end

% And remove input, as we do not need it anymore
clear input f i
% Conversion from char array to boolean
problem.anglesInRadians = strcmp(problem.anglesInRadians,'true') || ...
    strcmp(problem.anglesInRadians,'True');
% Computation of the lambda for the frequency read from the input file
problem.lambda = problem.c/problem.freq;
% Spacing between antennas x-axis
problem.dx = problem.lambda/(problem.invdAntennas*problem.gamma);
% Spacing between antennas y-axis
problem.dy = problem.lambda/(problem.invdAntennas*problem.gamma);
