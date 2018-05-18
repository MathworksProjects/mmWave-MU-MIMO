function [solutions] = o_test_MUMIMO_antenna_allocation(conf_file,problem_file)
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
    solutions = cell(length(meta_problem.nUsers),...
                    length(meta_problem.N_Antennas),...
                    length(meta_problem.maxnChannelPaths),...
                    meta_problem.realizations,...
                    4); % Equally sized matrix, different geometries
    
    for r=1:meta_problem.realizations
        for u=1:length(meta_problem.nUsers)
            tmp_problem = meta_problem;
            tmp_problem.nUsers = meta_problem.nUsers(u);
            [thetaPos, phiPos, dPos] = o_generate_positions(conf,tmp_problem);
            [MaxObjF,MinObjF] = o_generate_maxmin_throughput(meta_problem.nUsers(u),...
                meta_problem.maxMinObjF,meta_problem.minMinObjF,...
                meta_problem.MaxObjF);
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
                    problems{u,a,c,r}.MinObjF = MinObjF;
                    problems{u,a,c,r}.MaxObjF = MaxObjF;
                    problems{u,a,c,r}.thetaChannels = thetaChannels;
                    problems{u,a,c,r}.phiChannels = phiChannels;
                    problems{u,a,c,r}.alphaChannels = alphaChannels;
                    % arrayRestriction = 'None'
                    [sol_found,W,averageCap,totTime,usersAssigned] = ...
                        o_MUMIMO_antenna_allocation(conf,problems{u,a,c,r});
                    solutions{u,a,c,r,1}.sol_found = sol_found;
                    solutions{u,a,c,r,1}.W = W;
                    solutions{u,a,c,r,1}.averageCap = averageCap;
                    solutions{u,a,c,r,1}.totTime = totTime;
                    solutions{u,a,c,r,1}.usersAssigned = ...
                        length(usersAssigned);
                    if sol_found
                        fprintf('Solution found!\n');
                        fprintf('Av. cap = %f\n',averageCap);
                        fprintf('Exec. time = %f\n',totTime);
                        fprintf('Users assigned = %f\n\n',...
                            length(usersAssigned));
                    else
                        fprintf('No solution found!\n');
                    end
%                     % arrayRestriction = 'Localized'
%                     problems{u,a,c,r}.arrayRestriction = 'Localized';
%                     [sol_found,W,averageCap,totTime,usersAssigned] = ...
%                         o_MUMIMO_antenna_allocation(conf,problems{u,a,c,r});
%                     solutions{u,a,c,r,2}.sol_found = sol_found;
%                     solutions{u,a,c,r,2}.W = W;
%                     solutions{u,a,c,r,2}.averageCap = averageCap;
%                     solutions{u,a,c,r,2}.totTime = totTime;
%                     solutions{u,a,c,r,2}.usersAssigned = ...
%                         length(usersAssigned);
%                     if sol_found
%                         fprintf('Solution found!\n');
%                         fprintf('Av. cap = %f\n',averageCap);
%                         fprintf('Exec. time = %f\n',totTime);
%                         fprintf('Users assigned = %f\n\n',...
%                             length(usersAssigned));
%                     else
%                         fprintf('No solution found!\n');
%                     end
%                     % arrayRestriction = 'Interleaved'
%                     problems{u,a,c,r}.arrayRestriction = 'Interleaved';
%                     [sol_found,W,averageCap,totTime,usersAssigned] = ...
%                         o_MUMIMO_antenna_allocation(conf,problems{u,a,c,r});
%                     solutions{u,a,c,r,3}.sol_found = sol_found;
%                     solutions{u,a,c,r,3}.W = W;
%                     solutions{u,a,c,r,3}.averageCap = averageCap;
%                     solutions{u,a,c,r,3}.totTime = totTime;
%                     solutions{u,a,c,r,3}.usersAssigned = ...
%                         length(usersAssigned);
%                     if sol_found
%                         fprintf('Solution found!\n');
%                         fprintf('Av. cap = %f\n',averageCap);
%                         fprintf('Exec. time = %f\n',totTime);
%                         fprintf('Users assigned = %f\n\n',...
%                             length(usersAssigned));
%                     else
%                         fprintf('No solution found!\n');
%                     end
%                     % arrayRestriction = 'DiagInterleaved'
%                     problems{u,a,c,r}.NxSubarrays = meta_problem.nUsers(u);
%                     problems{u,a,c,r}.NySubarrays = meta_problem.nUsers(u);
%                     problems{u,a,c,r}.arrayRestriction = 'DiagInterleaved';
%                     [sol_found,W,averageCap,totTime,usersAssigned] = ...
%                         o_MUMIMO_antenna_allocation(conf,problems{u,a,c,r});
%                     solutions{u,a,c,r,4}.sol_found = sol_found;
%                     solutions{u,a,c,r,4}.W = W;
%                     solutions{u,a,c,r,4}.averageCap = averageCap;
%                     solutions{u,a,c,r,4}.totTime = totTime;
%                     solutions{u,a,c,r,4}.usersAssigned = ...
%                         length(usersAssigned);
                    if sol_found
                        fprintf('Solution found!\n');
                        fprintf('Av. cap = %f\n',averageCap);
                        fprintf('Exec. time = %f\n',totTime);
                        fprintf('Users assigned = %f\n\n',...
                            length(usersAssigned));
                    else
                        fprintf('No solution found!\n');
                    end
                end
            end
        end
    end
end

