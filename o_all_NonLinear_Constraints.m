function [c, ceq] = o_all_NonLinear_Constraints(X,conf,problem)
%ALL_NONLINEAR_CONSTRAINTS Wraps all non-linear constraints
%   Wrapping of all non-linear constraints for the MATLAB opt toolbox
    
    c(1) = o_amplitude_Greater_Than_Zero(X,conf,problem);
    if ~strcmp(conf.genStructure, 'nchoosek')
        c(2) = o_antenna_combination_constraint(X,conf,problem);
        if ~strcmp(conf.genStructure, 'allAntennas')
            c(3) = o_position_Constraint(X,problem);
        end
    end
    ceq = [];
end

