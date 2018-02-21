%% Script containing the solver configuration parameters
% Algorithm used for MATLAB optimization
% Genetic -> 'GA'; PatternSearch -> 'PS'; Particle Swarm Optimization -> 'PSO'
conf.algorithm = 'PSO';
% Sets whether we consider multi-path channel or LOS channel. In the first
% case, we need the information of the paths in the input file
conf.multiPath = true;
% Debug boolean
conf.DEBUG = false;
% Sets whether we want to compute the beam pattern for every single
% solution (MUST be true if conf.PlotAssignments == true) (TIME consuming)
conf.IncludeBeamPattern = false;
% If set to true, the Alg. 1 solution is random
conf.randomSolution = false;
% Sets whether we want to refine the solution given for a certain group of
% users when there is room for improvement
conf.RefineSolution = false;
% Sets whether we want to plot the polar plots of the antenna radiation
% patterns for each user (function parse_antenna_assignment). You should check
% if conf.PlotAssignments is TRUE. (TIME consuming)
conf.PlotAssignments = false;
% Sets whether we want to plot all the combinations in get_all_feasible_combs
conf.PlotCombinations = false;
% Sets wether we want to plot the initial and final assignments for every
% <user,Nmax> combination in Alg. 1
conf.plotAssignmentInitialAndFinal = false;
% Sets whether we want to plot the antenna selection displacements for
% every <user,Nmax> combination in Alg. 2 (getComb function) (TIME and
% MEMORY consuming)
conf.PlotDisplacements = false;
% Sets whether we want to compare our internal beamwidth function with
% findLobes computations
conf.findLobes = false;
% Sets the data structure for the Algorithm 1 solver (either genetic or
% pattern) 'onlyAssigned' recommended
conf.genStructure = 'onlyAssigned'; %'nchoosek' / 'onlyAssigned' / 'allAntennas'
% Sets the initial figure Id for combination of assignments' plots (to
% avoid collisions with radiation pattern plots)
conf.figIdx = 100;
% Sets the minimum number of users to be considered when exploring all the
% combinations possible (script get_all_feasible_combs)
conf.nUsersMin = 1;
% Color list for the plots
% http://www.rapidtables.com/web/color/RGB_Color.htm
conf.colorList = {[255 51 51]./255, [128 255 0]./255,...
             [51 255 255]./255, [0 128 255]./255, [0 0 255]./255,...
             [127 0 255]./255, [255 0 255]./255, [255 0 127]./255,...
             [102 0 0]./255, [102 51 0]./255, [102 102 0]./255,...
             [51 102 0]./255, [0 51 102]./255, [0 0 102]./255,...
             [51 0 102]./255, [102 0 102]./255, [102 0 51]./255};
% Color for unused antennas
conf.colorEmpty = [152 152 152]./255;