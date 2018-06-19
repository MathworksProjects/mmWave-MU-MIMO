function [Cap,SINR,estObj] = f_heuristicsDummy(MinObjF,MinObjFIsSNR,snrRange)
% F_HEURISICSDUMMY - Dummy function that emulates real heuristics.
%
% Syntax:  [Cap,SINR,estObj] = ...
%                          f_heuristicsDummy(MinObjF,MinObjFIsSNR,snrRange)
%
% Inputs:
%    MinObjF - Minimum value required
%    MinObjFIsSNR - True if optimization is based on SINR. False if Cap.
%    snrRange - Possible SINR range
%
% Outputs:
%    Cap - Vector [nUsers x 1] containing the capacity values (Linear)
%    SINR - Vector [nUsers x 1] containing the SINR values (Linear)
%    estObj - Vector [nUsers x 1] containing either the Capacity or SINR 
%             values (Linear), dpending upon MinObjFIsSNR
%
% Example: 
%   problem = o_read_input_problem('data/metaproblem_test.dat');
%   conf = o_read_config('data/config_test.dat');
%   conf.verbosity = 1;  % To visualize metrics on command line
%   problem.nUsers = 5;  % Fix number of users manually for example
%   [problem,~,~] = f_configuration(conf,problem);
%   [Cap,SNRList,~] = f_heuristicsDummy(problem.MinObjF,...
%                               conf.MinObjFIsSNR,problem.MCSPER.snrRange);
%
% Other m-files required: f_configuration
% Subfunctions: none
% MAT-files required: none
% DAT-files required: data/metaproblem_test.dat and data/config_test.dat
%
% See also: main.m,  f_conventionalBF.m,  f_heuristics.m,  main_runnable.m

%------------- BEGIN CODE --------------
    
if MinObjFIsSNR
%         MinObjFdB = max(10*log10(MinObjF),min(snrRange));  % in dB
%         MaxObjFdB = 30;  % in dB
%         estObjdB  = rand(1,length(MinObjFdB))*range([MinObjFdB MaxObjFdB]) + MinObjFdB;  % SNR (dB)
%         estObj    = 10.^(estObjdB/10);  % SNR (linear)
    % Uncomment line below for nice plotting of transmission bits
    estObj = MinObjF;  % SINR (linear) 

    SINR = estObj;  % SINR (Linear)
    Cap = log2(estObj + 1);  % Capacity (Linear)
else
    estObj = MinObjF;  % Capacity in bps/s/Hz

    Cap = estObj;  % Capacity (Linear)
    SINR = 2^estObj - 1;  % SINR (Linear)
end


% EOF