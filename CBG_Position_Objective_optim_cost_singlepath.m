function Val = CBG_Position_Objective_optim_cost_singlepath(myGene, conf, problem)
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
end