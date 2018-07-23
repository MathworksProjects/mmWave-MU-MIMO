function [handle_Conf_Array,W,Dir_OK,Dir_NOK,Cap_lin,SINR_PB_lin] = ...
                 CBG_geneToAssignment(myGene,problem,conf)
    handle_Conf_Array = phased.ConformalArray('Element',problem.handle_Ant,...
                        'ElementPosition',[problem.possible_locations(1,:);...
                                           problem.possible_locations(2,:);...
                                           problem.possible_locations(3,:)]);
    
    possible_locations = handle_Conf_Array.getElementPosition;
    elementPos = possible_locations./problem.lambda;  % Normalized
    
    W_LCMV = zeros(length(problem.nUsers),length(problem.N_Antennas));
    W_CBF = zeros(length(problem.nUsers),length(problem.possible_locations));
    nAntennas = length(myGene)/2;  % total number of antennas to distribute
    for id = 1:problem.nUsers
        antennaSelected = find(myGene(nAntennas+1:2*nAntennas)==id);
        elementPos1 = elementPos(:,antennaSelected);
        angles = [-problem.phiUsers ; -problem.thetaUsers];
        sv = steervec(elementPos1,angles);
        Sn = eye(length(antennaSelected));
        resp = zeros(problem.nUsers,1);
        resp(id) = db2pow(33);
        % Call LCMV
        W_LCMV1 = lcmvweights(sv,resp,Sn);
        W_LCMV(id,antennaSelected) = W_LCMV1.';
        % Call Conventional (CBF)
        W2 = cbfweights(elementPos1,angles(:,id));
        W_CBF(id,antennaSelected) = W2.';
        % Normalize weigths
        W_LCMV(id,:) = (1/sqrt(W_LCMV(id,:)*W_LCMV(id,:)'))*W_LCMV(id,:);
        W_CBF(id,:) = (1/sqrt(W_CBF(id,:)*W_CBF(id,:)'))*W_CBF(id,:);
    end

    % Asign final configuration
    W = W_LCMV;

    % Retrieve results
    [Dir_OK,Dir_NOK,Cap_lin,SINR_PB_lin]  = f_BF_results(W_LCMV,handle_Conf_Array,problem.candSet,problem,conf,false);
end