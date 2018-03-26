function [x,fval,exitflag,output] = o_CA_Position_Objective_optim_pattern(x0,conf,problem,lb,ub,StepTolerance_Data,FunctionTolerance_Data,MaxIterations_Data)
%% This is an auto generated MATLAB file from Optimization Tool.
% Copyright 2017  The MathWorks, Inc.
%% Start with the default options
options = optimoptions('patternsearch');
%% Modify options setting
options = optimoptions(options,'StepTolerance', StepTolerance_Data);
options = optimoptions(options,'FunctionTolerance', FunctionTolerance_Data);
options = optimoptions(options,'MeshTolerance',1e-1);
options = optimoptions(options,'MaxIterations', MaxIterations_Data);
options = optimoptions(options,'AccelerateMesh', true);
options = optimoptions(options,'UseCompletePoll', true);
options = optimoptions(options,'PollOrderAlgorithm', 'Consecutive');
options = optimoptions(options,'SearchFcn', @GPSPositiveBasis2N);
options = optimoptions(options,'UseCompleteSearch', true);
if conf.verbosity > 1
    options = optimoptions(options,'Display', 'iter');
elseif conf.verbosity <= 1
    options = optimoptions(options,'Display', 'off');
end
options = optimoptions(options,'PlotFcn', {@psplotbestf,@psplotfuncount});
%options = optimoptions(options,'PlotInterval', 1);
options = optimoptions(options,'UseVectorized', false);
options = optimoptions(options,'UseParallel', true);
if conf.multiPath
    [x,fval,exitflag,output] = ...
    patternsearch(@(x)o_Position_Objective_optim_cost_multipath(x,conf,problem),x0,[],[],[],[],lb,ub,@(x)o_all_NonLinear_Constraints(x,conf,problem),options);
else
    [x,fval,exitflag,output] = ...
    patternsearch(@(x)o_Position_Objective_optim_cost_singlepath(x,conf,problem),x0,[],[],[],[],lb,ub,@(x)o_all_NonLinear_Constraints(x,conf,problem),options);
end
end
