function val = o_antenna_combination_constraint(X,conf,problem)
%ANTENNA_COMBINATION_CONSTRAINT Non-linear constraint for Antenna Subarray matching
%   Non-linear constraint to detect antennas placed in
%   the same location
    if strcmp(conf.genStructure, 'allAntennas')
        antenna_selection = find(X(1:round(end/2)));
        num_subarrays_selected = numel(antenna_selection);
        val = (num_subarrays_selected ~= problem.Nmax);
    else
        u = unique(X(1:problem.Nmax)');
        v = problem.Nmax; %v = sum(X(1:conf.Nmax) == sort(X(1:conf.Nmax))); %
        val = size(u,1) < problem.Nmax || v < problem.Nmax; % Duplicates found if not equal (val == 1)
    end
end

