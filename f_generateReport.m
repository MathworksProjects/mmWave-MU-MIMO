function [ratioOK,ratioNOK] = f_generateReport(flows, DEBUG)
%F_GENERATEREPORT - Generate a report containing the number of PHY-flows
%that delivered the bits to the receiver before the application deadline.
%
% Syntax:  [ratioOK,ratioNOK] = f_generateReport(flows)
%
% Inputs:
%    flows - Description
%
% Outputs:
%    ratioOK - Description
%    ratioNOK - Description
%
% Example: 
%    [flows,TXbitsTot,THTot,lastSlotSimm,lastSelFlow] = main([]);
%    [ratioOK,ratioNOK] = f_generateReport(flows);
%
% Other m-files required: main.m
% Subfunctions: none
% MAT-files required: none
%
% See also: main,  main_runnable,  main_plotting

%------------- BEGIN CODE --------------

Nusers = length(flows);
ratioOK = zeros(Nusers,1);
ratioNOK = zeros(Nusers,1);
fprintf('\n***** REPORT *****\n');
for id = 1:Nusers
    totOK = sum(flows(id).success);
    totNOK = sum(flows(id).failed);
    ratioOK(id) = 100 * totOK / (totOK+totNOK);
    ratioNOK(id) = 100 * totNOK / (totOK+totNOK);
    if DEBUG
        fprintf('ID=%d\n\tOK=%.2f(%%)\n\tNOK=%.2f(%%)\n',...
                                              id,ratioOK(id),ratioNOK(id));
    end
end
if DEBUG
    fprintf('Overall Performance: OK=%.2f(%%)\n',mean(ratioOK));
    fprintf('******************\n');
end

%EOF