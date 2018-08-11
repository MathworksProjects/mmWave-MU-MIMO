function [problem] = o_read_input_problem(problem_file)
    %% Sets the configuration parameters (filepath, debug boolean and others)
    %%%%  and read the input file to obtain the initial parameters 
    %%%%  used for the Algorithm 1.
    %   After this script is executed, several variables are available
    %    in a struct that can be used for further processing

    % Arguments that will be read by the readDATInputData function
    % We read first all those that are not related to the number of users
    inputArgList = { struct('name','N_Antennas','type','interval');
                struct('name','NxSubarrays','type','int');
                struct('name','NySubarrays','type','int');
                struct('name','NmaxArray','type','interval');
                struct('name','arrayRestriction','type','char');
                struct('name','invdAntennas','type','float');
                struct('name','freq','type','float');
                struct('name','gamma','type','float');
                struct('name','nUsers','type','interval');
                struct('name','maxnUsers','type','float');
                struct('name','maxnChannelPaths','type','interval');
                struct('name','thetaUsers','type','floatArray','size','maxnUsers');
                struct('name','phiUsers','type','floatArray','size','maxnUsers');
                struct('name','dUsers','type','floatArray','size','maxnUsers');
                struct('name','anglesInRadians','type','bool');
                struct('name','mindUsers','type','float');
                struct('name','maxdUsers','type','float');
                struct('name','MaxObjF','type','float');
                struct('name','maxMinObjF','type','float');
                struct('name','minMinObjF','type','float');
                struct('name','Noise','type','float');
                struct('name','Ptx','type','float');
                struct('name','Bw','type','float');
                struct('name','c','type','float');
                struct('name','mcsPolicy','type','char');
                struct('name','realizations','type','int');
                struct('name','trafficType','type','char');
                struct('name','iat','type','float');
                struct('name','payload','type','float');
                struct('name','manuallyAssignApp','type','bool');
                struct('name','FLAGagg','type','bool');
                struct('name','DEBUG','type','bool');
                struct('name','PLOT_DEBUG','type','bool');
                struct('name','numPkts','type','float');
                struct('name','Tslot','type','float');
                struct('name','Tsym','type','float');
                struct('name','targetPER','type','float');
                struct('name','BFalgorithm','type','char');
                };
    problem = o_readDATInputData(problem_file,inputArgList);
    
    problem.NxPatch = floor(sqrt(problem.N_Antennas));
    problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);
    % Adjusted number of antennas
    problem.N_Antennas = problem.NxPatch.*problem.NyPatch;
    % Computation of the lambda for the frequency read from the input file
    problem.lambda = problem.c/problem.freq;
    % Spacing between antennas x-axis
    problem.dx = problem.lambda/(problem.invdAntennas*problem.gamma);
    % Spacing between antennas y-axis
    problem.dy = problem.lambda/(problem.invdAntennas*problem.gamma);
    % Number of subarrays
    problem.N_Subarrays = problem.NxSubarrays*problem.NySubarrays;
end
