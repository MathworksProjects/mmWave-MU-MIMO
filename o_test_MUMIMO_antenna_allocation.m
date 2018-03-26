function [] = o_test_MUMIMO_antenna_allocation(conf_file,problem_file)
%TEST_MUMIMO_ANTENNA_ALLOCATION Function to test MUMIMO antenna allocation
%   conf_file is a path to a file in DAT format with the configuration 
%   parameters for the execution
%   problem_file is a path to a file in DAT format with the problem 
%   parameters for the execution
    conf = o_read_config(conf_file);
    % Let's process the input parameters and read the input file (problem). Edit 
    % this script to change the parameters
    meta_problem = o_read_input_problem(problem_file);
    problems = cell(length(meta_problem.nUsers),...
                    length(meta_problem.N_Antennas),...
                    length(meta_problem.maxnChannelPaths),...
                    meta_problem.realizations);
    solutions = problems; % Equally sized cell
    
    for r=1:meta_problem.realizations
        for u=1:length(meta_problem.nUsers)
            [thetaPos, phiPos, dPos] = o_generate_positions(conf,meta_problem.nUsers(u),...
                meta_problem.maxdUsers,meta_problem.mindUsers);
            [MaxThr,MinThr] = o_generate_maxmin_throughput(meta_problem.nUsers(u),...
                meta_problem.maxMinThr,meta_problem.minMinThr,...
                meta_problem.MaxThr);
            for c=1:length(meta_problem.maxnChannelPaths)
                [thetaChannels, phiChannels, alphaChannels] = ...
                    o_generate_channels(conf,meta_problem.nUsers(u),...
                    meta_problem.maxnChannelPaths(c));
                for a=1:length(meta_problem.N_Antennas)
                    fprintf('u:%d - a:%d - c:%d - r:%d\n',...
                        meta_problem.nUsers(u),...
                        meta_problem.NxPatch(a)*meta_problem.NyPatch(a),...
                        meta_problem.maxnChannelPaths(c),...
                        r);
                    problems{u,a,c,r} = meta_problem;
                    problems{u,a,c,r}.nUsers = meta_problem.nUsers(u);
                    problems{u,a,c,r}.N_Antennas = ...
                        meta_problem.NxPatch(a)*meta_problem.NyPatch(a);
                    % NmaxArray relative to N_Antennas
                    problems{u,a,c,r}.NmaxArray = ...
                        floor(problems{u,a,c,r}.N_Antennas * ...
                        meta_problem.NmaxArray);
                    problems{u,a,c,r}.NxPatch=meta_problem.NxPatch(a);
                    problems{u,a,c,r}.NyPatch=meta_problem.NyPatch(a);
                    problems{u,a,c,r}.maxnChannelPaths = meta_problem.maxnChannelPaths(c);
                    problems{u,a,c,r}.thetaUsers = thetaPos;
                    problems{u,a,c,r}.phiUsers = phiPos;
                    problems{u,a,c,r}.dUsers = dPos;
                    problems{u,a,c,r}.MinThr = MinThr;
                    problems{u,a,c,r}.MaxThr = MaxThr;
                    problems{u,a,c,r}.thetaChannels = thetaChannels;
                    problems{u,a,c,r}.phiChannels = phiChannels;
                    problems{u,a,c,r}.alphaChannels = alphaChannels;
                    [averageCap,totTime] = o_MUMIMO_antenna_allocation(conf,problems{u,a,c,r});
                    solutions{u,a,c,r}.averageCap = averageCap;
                    solutions{u,a,c,r}.totTime = totTime;
                end
            end
        end
    end
end

