function [opt_solution,min_value] = o_CA_Position_Objective_optim_exhaustive(conf,...
                                                        problem,lb,ub)
%% Start with the default options
    min_value = Inf;
    opt_solution = zeros(1,problem.Nmax*3+1);
    availableAntennas = 1:problem.N_Subarrays;  %THIS IS YOUR VECTOR
    N = problem.N_Subarrays;            %LENGTH OF V
    K = problem.Nmax;           
    m = 0;
    h = K;
    iteration = 1;         
    A = 1:K;          %FIRST COMBINATION
    while(iteration>0)
        [A,m,h,iteration] = GetNextCombination(N, K, A, m, h, iteration);
        selectedAntennas = availableAntennas(A);
        
        amplitudesIndexes = zeros(1,problem.Nmax+1);
        first_iteration = true;
        while (sum(amplitudesIndexes ~= zeros(1,problem.Nmax+1)) || first_iteration)
            first_iteration = false;
            amplitudes = problem.fixedPointScale*(amplitudesIndexes(1:problem.Nmax) + 1);
            Fbb = problem.fixedPointScale*(amplitudesIndexes(end)+1);
            tmp_solution = [selectedAntennas,amplitudes,lb((2*problem.Nmax+1):(3*problem.Nmax)),Fbb];
            if conf.multiPath
                obj_value = o_Position_Objective_optim_cost_multipath(tmp_solution,conf,problem);
            else
                obj_value = o_Position_Objective_optim_cost_singlepath(tmp_solution,conf,problem);
            end
            if obj_value < min_value
                min_value = obj_value;
                opt_solution = tmp_solution;
                disp(min_value);
                disp(opt_solution);
            end
            amplitudesIndexes = (o_sumOneToCombination(amplitudesIndexes',ones(1,problem.Nmax+1)*(2^conf.NbitsAmplitude-1)))';
        end
    end
end