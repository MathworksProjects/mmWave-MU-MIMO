function [W_LCMV,W_CBF,handle_ConformalArray] = f_conventionalBF(problem,candSet)
    % Restrict sub-arrays to Localized for LCMV
    problem.arrayRestriction = 'Localized';
    % Compute number of sub-arrays to assign per user
    problem.NySubarrays = 2;
    problem.NxSubarrays = problem.nUsers/problem.NySubarrays;
    problem = o_compute_antennas_per_user(problem,candSet);
    % Create subarray partition
    problem = o_create_subarray_partition(problem);
    totSubArrays = (problem.NxSubarrays * problem.NySubarrays) / problem.nUsers;
    problem.NzPatch = problem.NxPatch;
    problem.dz = problem.dx;
    problem.handle_Ant = phased.CosineAntennaElement('FrequencyRange',...
                            [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
                            'CosinePower',[1.5 2.5]); % [1.5 2.5] values set porque si
    handle_ConformalArray = phased.URA([problem.NyPatch,problem.NzPatch],...
                              'Lattice','Rectangular','Element',problem.handle_Ant,...
                              'ElementSpacing',[problem.dy,problem.dz]);
    problem.possible_locations = handle_ConformalArray.getElementPosition;
    
    % Antennas assigned to each user (fixed)
    mySubArray = (1:1:(problem.NxSubarrays * problem.NySubarrays));
    relevant_positions = cell(problem.nUsers,1);
    for valID = 1:problem.nUsers
        partAssignation = mySubArray(1:totSubArrays);
        temp = [];
        for ass = partAssignation
            temp = [temp problem.Partition{ass}];  %#ok<AGROW>
        end
        relevant_positions{valID} = temp;
        mySubArray(mySubArray==partAssignation) = [];  % delete assigned antennas
    end
    
    % Compute weights (beamforming)
    W_LCMV = zeros(problem.nUsers, problem.NxPatch*problem.NyPatch);
    W_CBF = zeros(problem.nUsers, problem.NxPatch*problem.NyPatch);
    for id = 1:problem.nUsers
        Nant_user = length(relevant_positions{id});
        % Convert antenna ID's into physical locations
        elementPos = [problem.possible_locations(1,relevant_positions{id});...
                      problem.possible_locations(2,relevant_positions{id});...
                      problem.possible_locations(3,relevant_positions{id})];
        elementPosNorm = elementPos./problem.lambda;
        
        PhiTheta = ([-problem.phiUsers ; -problem.thetaUsers]);
        
        % Apply LCMV Beamformer for selected user
        sv = steervec(elementPosNorm,PhiTheta);
        Sn = eye(Nant_user);
        resp = zeros(problem.nUsers,1) + eps;
        resp(id) = 1;  % Maximum restricted to limit (33dB)
        w_lcmv = lcmvweights(sv,resp,Sn);  % LCMV Beamformer method
        
        % Apply Convencional Beamformer for selected user
        w_cbf = cbfweights(elementPosNorm,PhiTheta(:,id));  % conventional beamformer
        
        % Store results in global W
        W_LCMV(id,relevant_positions{id}) = w_lcmv.';
        W_CBF(id,relevant_positions{id}) = w_cbf.';
    end
end