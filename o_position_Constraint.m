function val = o_position_Constraint(X,problem)
%POSITION_CONSTRAINT Non-linear constraint for Antenna Subarray matching
%   	The positions must be integer values
    act_pos = X(1:problem.Nmax);
    val = sum(act_pos == floor(act_pos)) < problem.Nmax;
    % Double values found if sum < N_MaxAnt (val == 1)
end

