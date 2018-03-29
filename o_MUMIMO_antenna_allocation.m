function [sol_found,W,averageCap,totTime,usersAssigned] = o_MUMIMO_antenna_allocation(conf,problem)
%% Antenna allocation solver for massive arrays in mm-Wave and THz-band
    % This script is the main program of the MATLAB solver for the antenna
    % allocation problem in massive arrays for mm-Wave and THz-band.
    % Version 2, not backward-compatible.  

    % We only need to execute the corresponding method

    if strcmp(conf.solver,'LuckAndChoice')
        [sol_found,W,averageCap,totTime,usersAssigned] = ...
            o_luck_and_choice_solver(problem,conf);
    elseif strcmp(conf.solver,'BestOrRemove')
        [sol_found,W,averageCap,totTime,usersAssigned] = ...
            o_best_or_remove_solver(problem,conf);
    end
end

