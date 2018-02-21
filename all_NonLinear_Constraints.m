function [c, ceq] = all_NonLinear_Constraints(X,conf,problem)
%ALL_NONLINEAR_CONSTRAINTS Wraps all non-linear constraints
%   Wrapping of all non-linear constraints for the MATLAB opt toolbox
    
    c(1) = amplitude_Greater_Than_Zero(X,conf,problem);
    if ~strcmp(conf.genStructure, 'nchoosek')
        c(2) = antenna_combination_constraint(X,conf,problem);
        if ~strcmp(conf.genStructure, 'allAntennas')
            c(3) = position_Constraint(X,problem);
        end
    end
    ceq = [];
end

