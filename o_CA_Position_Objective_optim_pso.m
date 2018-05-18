function [x,fval,exitflag,output] = o_CA_Position_Objective_optim_pso(x0,conf,problem,lb,ub,InitialPopulation_Data,MaxIterations_Data)
%% This is an auto generated MATLAB file from Optimization Tool.
% Copyright 2017  The MathWorks, Inc.
%% Add path
addpath('./constrainedpso');
%% Start with the default options
if conf.verbosity > 1
    options = psooptimset('InitialPopulation',InitialPopulation_Data,...
                      'UseParallel','always',...
                      'Display','diagnose',...
                      'DemoMode', 'pretty',...
                      'PlotFcns', {@psoplotbestf},...
                      'ConstrBoundary', 'soft');
elseif conf.verbosity == 1
        options = psooptimset('InitialPopulation',InitialPopulation_Data,...
                      'UseParallel','always',...
                      'Display','final',...
                      'DemoMode', 'pretty',...
                      'PlotFcns', {@psoplotbestf},...
                      'ConstrBoundary', 'soft');
elseif conf.verbosity == 0
        options = psooptimset('InitialPopulation',InitialPopulation_Data,...
                      'UseParallel','always',...
                      'Display','off',...
                      'DemoMode', 'pretty',...
                      'PlotFcns', {@psoplotbestf},...
                      'ConstrBoundary', 'soft');
end
%     AccelerationFcn: [Function handle | {@psoiterate}]
% CognitiveAttraction: [Positive scalar | {0.5}]
%      ConstrBoundary: ['soft' | 'penalize' | 'reflect' | 'absorb' | {'penalize'}]
%             Display: ['off' | 'final' | 'diagnose' | {'final'}]
%            DemoMode: ['fast' | 'pretty' | 'on' | 'off' | {'off'}]
%        FitnessLimit: [Scalar | {-Inf}]
%         Generations: [Positive integer | {200}]
%           HybridFcn: [@fminsearch | @patternsearch | @fminunc | @fmincon | {[]}]
%   InitialPopulation: [empty matrix | nxnvars matrix | {[]}]
%   InitialVelocities: [empty matrix | nxnvars matrix | {[]}]
%            PlotFcns: [Cell array of fcn handles | {{}}]
%        PlotInterval: [Positive integer | {1}]
%        PopInitRange: [2x1 vector | 2xnvars matrix | {[0;1]}]
%      PopulationSize: [Positive integer | {40}]
%      PopulationType: ['bitstring' | 'doubleVector' | {'doubleVector'}]
%    SocialAttraction: [Positive scalar | {1.25}]
%       StallGenLimit: [Positive integer | {50} ]
%      StallTimeLimit: [Positive scalar (seconds) | {Inf} ]
%           TimeLimit: [Positive scalar (seconds) | {50} ]
%              TolFun: [Positive scalar | {1e-06}]
%              TolCon: [Positive scalar | {1e-06}]
%         UseParallel: ['always' | 'never' | {'never'}]
%          Vectorized: ['on' | 'off' | {'off'}]
%       VelocityLimit: [Positive scalar | {[]}]
nvars = numel(lb); % == numel(ub) == num. variables in the gen
if conf.multiPath
    [x,fval,exitflag,output,population,scores] = ...
    pso(@(x)o_Position_Objective_optim_cost_multipath(x,conf,problem),nvars,[],[],[],[],lb,ub,@(x)o_all_NonLinear_Constraints(x,conf,problem),options);
else
    [x,fval,exitflag,output,population,scores] = ...
    pso(@(x)o_Position_Objective_optim_cost_singlepath(x,conf,problem),nvars,[],[],[],[],lb,ub,@(x)o_all_NonLinear_Constraints(x,conf,problem),options);
end
