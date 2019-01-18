function Val = CBG_Position_Objective_optim_cost_singlepath(myGene, conf, problem)
% CBG_Position_Objective_optim_cost_singlepath - It computes the quality of
% the current genes (GA), which maps to an actual antenna allocation per
% user. First, it converts the genes into actual allocations. Second, it
% computes the quality value of such allocation using transfer Scoring
% functions, implemented under o_transferScore.
%
% Syntax:  Val = CBG_Position_Objective_optim_cost_singlepath(myGene, conf, problem)
%
% Inputs:
%    myGene - output configuration from ga in the form of genes
%    conf - struct containint configuration in data/config_test.dat
%    InitialPopulation - Initial antenna configuration
%
% Outputs:
%    Val - Value that reflects the quality of the current antenna
%    allocation configuration, dependent on the transferScore functions,
%    which maps directivities (towards intended and unintended users) to
%    quality.
%
% Other m-files required: CBG_geneToAssignment, o_transferScore
% Subfunctions: p_crossoverArrayGA, p_mutationArrayGA
% MAT-files required: none
%
% See also: CBG_solveit , o_transferScore, CBG_geneToAssignment

%------------- BEGIN CODE --------------
% Convert current gene into generated directivities
[~,~,DirOK,DirNOK,~,~] = CBG_geneToAssignment(myGene,problem,conf);

% Compute generated interference
DirNOK_lin = db2pow(DirNOK);
DirNOK_pcvd_lin = sum(DirNOK_lin,1);
DirNOK_pcvd = pow2db(DirNOK_pcvd_lin);

Val = 0;
for id = 1:problem.nUsers
    % Compute scores per user
    scorePrx = o_transferScore(DirOK(id),1);
    scoreInt = o_transferScore(DirNOK_pcvd(id),0);
    % Compute final score
    Val = Val + (-1)*( conf.Fweights(1)*scorePrx + conf.Fweights(2)*scoreInt );
end

% Display score and performance metrics
if conf.verbosity > 1
    fprintf('Val = %f\n',Val);
    fprintf('scorePrx = %.2f\tscoreInt = %.2f\n',scorePrx,scoreInt);
    fprintf('Dir(dB)= %.2f\tDir_int(dB)= %.2f\n',DirOK,DirNOK);
end



% EOF