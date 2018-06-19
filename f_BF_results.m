function [DirOK,DirNOK,Cap_lin,SINR_PB_lin]  = f_BF_results(W,handle_ConformalArray,problem,conf,plotFLAG)
% f_BF_RESULTS - 
%
% Syntax:  [DirOK,DirNOK,Cap,SINR_PB,estObj]  = ...
%               f_BF_results(W,handle_ConformalArray,problem,conf,plotFLAG)
% Inputs:
%    W - Matrix [nUser x nAntennas] containg the weights for LCMV
%    handle_ConformalArray - Initial conformal array
%    problem - struct containint configuration in data/metaproblem_test.dat
%    conf - Struct containing configuration in data/config_test.dat
%    plotFLAG - True for plotting the beam pattern. False otherwise.
%
% Outputs:
%    DirOK - Vector [nUser x 1] containg Directivity to target user (dB)
%    DirNOK - Matrix [nUser x nUser] containg the Interference generated to
%             the other users (dB)
%    Cap_lin - Capacity achieved per user (in bits/Hz/s)
%    SINR_PB_lin - SINR achieved per user (in Linear scale)
%
% Example: 
%   problem = o_read_input_problem('data/metaproblem_test.dat');
%   conf = o_read_config('data/config_test.dat');
%   conf.verbosity = 1;  % To visualize metrics on command line
%   problem.nUsers = 5;  % Fix number of users manually for example
%   [problem,~,~] = f_configuration(conf,problem);
%   problem.N_Antennas = nAntennas;  % Select number of antennas
%   problem.NxPatch = floor(sqrt(problem.N_Antennas));  % Adjust
%   problem.NyPatch = floor(problem.N_Antennas./problem.NxPatch);  % Adjust
%   problem.N_Antennas = problem.NxPatch.*problem.NyPatch;  % Adjust
%   candSet = [2 4 5];  % Users 1 and 3 are left out
%   [W_LCMV,W_CBF,handle_ConformalArray,estObj_LCMV,estObj_CBF] = ...
%                                    f_conventionalBF(problem,conf,candSet)
%   [DirOK,DirNOK,Cap,SINR_PB]  = ...
%                   f_BF_results(W,handle_ConformalArray,problem,conf,true)
%
% Other m-files required: f_conventionalBF or f_heuristics
% Subfunctions: none
% MAT-files required: none
% DAT-files required: none
%
% See also: main.m,  f_conventionalBF.m,  f_heuristics.m,  main_runnable.m

%------------- BEGIN CODE --------------

% Output parameters
DirOK = -Inf(problem.nUsers,1);  % Directivity target (heuristics)
DirNOK = -Inf(problem.nUsers,problem.nUsers);  % Directivity others (heuristics)

% Antenna location in the array
possible_locations = handle_ConformalArray.getElementPosition;

for id = 1:problem.nUsers
    relevant_positions = (W(id,:)~=0);
    Taper_user = W(id,relevant_positions);

    handle_Conf_Array_USER = phased.ConformalArray(...
                          'Element',handle_ConformalArray.Element,...
                          'ElementPosition', [possible_locations(1,relevant_positions);...
                                              possible_locations(2,relevant_positions);...
                                              possible_locations(3,relevant_positions)],...
                          'Taper',Taper_user);

    % Extract Rx Power (in dB)
    DirOK(id) = patternAzimuth(handle_Conf_Array_USER,problem.freq,problem.thetaUsers(id),'Azimuth',problem.phiUsers(id),'Type','powerdb');
    % Extract interference generated to others (in dB)
    for id1 = 1:1:problem.nUsers
        if id1~=id
            DirNOK(id,id1) = patternAzimuth(handle_Conf_Array_USER,problem.freq,problem.thetaUsers(id1),'Azimuth',problem.phiUsers(id1),'Type','powerdb');
        end
    end

    if plotFLAG
        % Plot beamforming per user
        problem.IDUserAssigned = id;
        o_plotAssignment_mod(problem, handle_Conf_Array_USER);
        % Plot assignation
        px = possible_locations(3,:);  % Antenna allocation on x-axis
        py = possible_locations(2,:);  % Antenna allocation on y-axis
        pz = possible_locations(1,:);  % Antenna allocation on z-axispatch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
        patch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
        arrays = o_getArrays(problem.nUsers,W,px,py,pz);
        o_plot_feasible_comb(problem,conf,patch,arrays);
    end
end

% Compute basic parameters for SINR and Capacity computations
chLoss_lin = ((problem.lambda ./ (4*pi*problem.dUsers(1:problem.nUsers))).^2).';  % Losses
Noise_lin = db2pow(problem.Noise);  % Noise power
Noise_lin = repmat(Noise_lin,problem.nUsers,1);
% Parse results for specific case
DirOK_lin = db2pow(DirOK);
DirNOK_lin = db2pow(DirNOK);
DirNOK_pcvd_lin = sum(DirNOK_lin,2); % Perceived interference
% Compute SINR and Capacities - LCMV
SINR_PB_lin = (DirOK_lin.*chLoss_lin) ./(DirNOK_pcvd_lin.*chLoss_lin + Noise_lin);  % Compute SINR Pass-Band (BB)
SINR_PB = pow2db(SINR_PB_lin);  % Compute SINR Pass-Band (BB)
Cap_lin = log2(1 + SINR_PB_lin);  % Compute final Capacity (bits/Hz/s)

if conf.verbosity >= 1
    for id = 1:problem.nUsers
        fprintf('* Capacity(%d): %.2f (bits/Hz/s)\n',id,Cap_lin(id));
        fprintf('* SINR(%d): %.2f (dB)\n',id,SINR_PB(id));
        fprintf('* Directivity IDmax: %.2f (dB)\n',DirOK(id));
        for id1 = 1:1:problem.nUsers
            if id1~=id
                fprintf('  Directivity IDmin(%d): %.2f (dB)\n',id1,DirNOK(id,id1));
            end
        end
    end
end
    
% EOF
