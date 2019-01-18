function val = o_amplitude_Greater_Than_Zero(X,conf,problem)
%AMPLITUDE_GREATER_THAN_ZERO Non-linear constraint for Antenna Subarray matching
%   Non-linear constraint to detect amplitudes equal to 0
    if strcmp(conf.genStructure, 'nchoosek')
        val = sum(X(2:problem.Nmax+1) == 0) > 0;
    elseif strcmp(conf.genStructure, 'allAntennas')
        val = sum(X(1:problem.N_Subarrays) == 0) < problem.Nmax;
    else
        % In this case we include X(end), i.e. Fbb
        val = sum([X(end),X((problem.Nmax+1):(2*problem.Nmax))] == 0) > 0;
    end
    % Zeros found if sum > 0 (val == 1)
end

