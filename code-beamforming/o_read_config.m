function [conf] = o_read_config(config_file)
    %% Script containing the solver configuration parameters
    if isempty(config_file)
        % Algorithm used for MATLAB optimization
        % Genetic -> 'GA'; PatternSearch -> 'PS'; Particle Swarm Optimization -> 'PSO'
        conf.algorithm = 'GA';
        % Sets whether f_heuristics will work with SINR or capacities
        conf.MinObjFIsSNR=true;
        % Sets whether we consider multi-path channel or LOS channel. In the first
        % case, we need the information of the paths in the input file
        conf.multiPath = true;
        % Verbosity (0 no output, 1 some informative output, 2 DEBUG mode)
        conf.verbosity = 0;
        % Sets whether we want to compute the beam pattern for every single
        % solution (MUST be true if conf.PlotAssignments == true) (TIME consuming)
        conf.IncludeBeamPattern = false;
        % If set to true, the Alg. 1 solution is random
        conf.randomSolution = false;
        % Weights used in the optimization function for each user alloc. subproblem
        conf.Fweights = [0.4 0.3 0];
        % Solver 'LuckAndChoice' or 'BestOrRemove'
        conf.solver = 'BestOrRemove';
        % Objective function for 'LuckAndChoice'
        conf.ObjFunc = 'compute_averageCap';
        % 
        conf.NoCut = false;
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
        % Delay profile for 5G channel used in simulations. Possible values
        % are ('CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E')
        conf.DelayProfile = 'CDL-A';
        % Sets the initial figure Id for combination of assignments' plots (to
        % avoid collisions with radiation pattern plots)
        % conf.figIdx = 100;
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
    else
        % Arguments that will be read by the readDATInputData function
        % We read first all those that are not related to the number of users
        inputArgList = {
                    % Basic parameters
                    struct('name','detLocation','type','bool');
                    struct('name','useCasesLocation','type','bool');
                    struct('name','useCaseLocation','type','float');
                    % Optimization parameters
                    struct('name','algorithm','type','char');
                    struct('name','Fweights','type','floatArray','size',3);
                    struct('name','solver','type','char');
                    struct('name','ObjFunc','type','char');
                    struct('name','RefineSolution','type','bool');
                    struct('name','findLobes','type','bool');
                    struct('name','randomSolution','type','bool');
                    struct('name','MinObjFIsSNR','type','bool');
                    struct('name','nUsersMin','type','int');
                    struct('name','genStructure','type','char');
                    struct('name','NoCut','type','bool');
                    struct('name','feasibility','type','bool');
                    % Extra Optimization parameters
                    struct('name','NumPhaseShifterBits','type','int');
                    struct('name','NbitsAmplitude','type','int');
                    struct('name','percMutatedGens','type','float');
                    struct('name','mutationImpact','type','float');
                    struct('name','PopulationSize_Data','type','float');
                    struct('name','EliteCount_Data','type','float');
                    struct('name','CrossoverFraction_Data','type','float');
                    struct('name','Maxgenerations_Data','type','float');
                    struct('name','MaxStallgenerations_Data','type','float');
                    struct('name','FunctionTolerance_Data','type','float');
                    struct('name','GArndImpact','type','float');
                    struct('name','GArndmodifyAmpl','type','bool');
                    % Channel parameters
                    struct('name','Use5GChannel','type','bool');
                    struct('name','multiPath','type','bool');
                    struct('name','DelayProfile','type','char');
                    % Debugging parameters
                    struct('name','verbosity','type','int');
                    % Plotting parameters
                    struct('name','PlotAssignments','type','bool');
                    struct('name','IncludeBeamPattern','type','bool');
                    struct('name','PlotCombinations','type','bool');
                    struct('name','PlotDisplacements','type','bool');
                    struct('name','figIdx','type','int');
                    struct('name','plotAssignmentInitialAndFinal','type','bool');
                    };
        conf = o_readDATInputData(config_file,inputArgList);
        conf.colorList = {[255 51 51]./255, [128 255 0]./255,...
                     [51 255 255]./255, [0 128 255]./255, [0 0 255]./255,...
                     [127 0 255]./255, [255 0 255]./255, [255 0 127]./255,...
                     [102 0 0]./255, [102 51 0]./255, [102 102 0]./255,...
                     [51 102 0]./255, [0 51 102]./255, [0 0 102]./255,...
                     [51 0 102]./255, [102 0 102]./255, [102 0 51]./255};
        % Color for unused antennas
        conf.colorEmpty = [152 152 152]./255;
    end
end

    