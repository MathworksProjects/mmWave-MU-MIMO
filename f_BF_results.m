function [DirOK,DirNOK]  = f_BF_results(W,handle_ConformalArray,problem,conf,plotFLAG)
    % Output parameters
	DirOK = -Inf(problem.nUsers,1);  % Directivity target (heuristics)
    DirNOK = -Inf(problem.nUsers,problem.nUsers);  % Directivity others (heuristics)
    % Antenna location in the array
    possible_locations = handle_ConformalArray.getElementPosition;

    for id = 1:problem.nUsers
        relevant_positions = (W(id,:)~=0);
        Taper_user = W(id,relevant_positions);

        handle_Conf_Array_USER = phased.ConformalArray(...
                              'Element',handle_ConformalArray.Element,...
                              'ElementPosition', [possible_locations(2,relevant_positions);...
                                                  possible_locations(2,relevant_positions);...
                                                  possible_locations(3,relevant_positions)],...
                              'Taper',Taper_user);

        % Extract Rx Power (in dB)
        DirOK(id) = patternAzimuth(handle_Conf_Array_USER,problem.freq,problem.thetaUsers(id),'Azimuth',problem.phiUsers(id),'Type','powerdb');
        fprintf('* Directivity IDmax: %.2f (dB)\n',DirOK(id));
        % Extract interference generated to others (in dB)
        for id1 = 1:1:problem.nUsers
            if id1~=id
                DirNOK(id,id1) = patternAzimuth(handle_Conf_Array_USER,problem.freq,problem.thetaUsers(id1),'Azimuth',problem.phiUsers(id1),'Type','powerdb');
                fprintf('  Directivity IDmin(%d): %.2f (dB)\n',id1,DirNOK(id,id1));
            end
        end
        
        
        if plotFLAG
            % Plot beamforming per user
            problem.IDUserAssigned = id;
            o_plotAssignment_mod(problem, handle_Conf_Array_USER);
            % Plot assignation
            px = possible_locations(3,:);  % Antenna allocation on x-axis
            py = possible_locations(2,:);  % Antenna allocation on y-axis
            pz = possible_locations(1,:);  % Antenna allocation on z-axispatch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
            patch = o_getPatch(problem.NxPatch,problem.NyPatch,px,py);
            arrays = o_getArrays(problem.nUsers,W,px,py,pz);
            o_plot_feasible_comb(problem,conf,patch,arrays);
        end
    end
end