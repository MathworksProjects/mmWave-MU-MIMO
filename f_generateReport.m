%% function [flows,selectedFlows] = f_generateReport(flows)
% This function blabla
%   INPUTS:
% - InputHere
%   OUTPUTS:
% - OutputHere
function [ratioOK,ratioNOK] = f_generateReport(flows)
    Nusers = length(flows);
    ratioOK = zeros(Nusers,1);
    ratioNOK = zeros(Nusers,1);
    fprintf('\n***** REPORT *****\n');
    for id = 1:Nusers
        totOK = sum(flows(id).success);
        totNOK = sum(flows(id).failed);
        ratioOK(id) = 100 * totOK / (totOK+totNOK);
        ratioNOK(id) = 100 * totNOK / (totOK+totNOK);
        fprintf('ID=%d\n\tOK=%.2f(%%)\n\tNOK=%.2f(%%)\n',...
                                              id,ratioOK(id),ratioNOK(id));
    end
    fprintf('Overall Performance: OK=%.2f(%%)\n',mean(ratioOK));
    fprintf('******************\n');
end