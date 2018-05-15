function [Val] = o_Position_Objective_optim_cost_multipath(Position_Taper, conf, problem)
    
    % Initialize optimization variables
    PR   = 0;  % Received power at intended user
    Int  = 0;  % Interference generated to other users
    HPBW = 0;  % Beamdiwth achieved
    
    % Define constants
    if ~strcmp(conf.genStructure, 'allAntennas')
        nck = strcmp(conf.genStructure, 'nchoosek');
        if sum(Position_Taper(1:problem.Nmax-(problem.Nmax*nck)+1) ~= ...
                floor(Position_Taper(1:problem.Nmax-(problem.Nmax*nck)+1))) > 0
            Val = Inf;
            % In the case where Position Taper is formed by any double, we shall discard it
            return
        end
    end

    % Create antenna subarray
    [Position_Taper_elem,problem] = o_subarrayToAntennaElements(Position_Taper,conf,problem);
    % Extracting taper values from the input vector
    Taper_value = Position_Taper_elem(problem.ant_elem+1:problem.ant_elem*2) .* ...
        exp(1i.*Position_Taper_elem(problem.ant_elem*2+1:problem.ant_elem*3));
    if isempty(Taper_value)
        fprintf('warning: Taper is empty.\n');
    end
    % Creating a Conformal Array with cosine elements.
    % Conformal array will be limited to a single plane
    handle_Conf_Array = phased.ConformalArray('Element',problem.handle_Ant,...
                        'ElementPosition',[zeros(1,problem.ant_elem);...
                        problem.possible_locations(2,Position_Taper_elem(1:problem.ant_elem));...
                        problem.possible_locations(3,Position_Taper_elem(1:problem.ant_elem))],...
                        'Taper',Taper_value);
    
    % Consider rx power and interference generated into opt. (IUI metric)
    if conf.Fweights(1)+conf.Fweights(2) > 0
        % Computing the received power in every user
        PRx = zeros(1,problem.nUsers);
        for u1=1:problem.nUsers
            for i=1:problem.maxnChannelPaths
                if problem.phiChannels(u1,i) ~= -Inf
                    % We assume that alphaChannels include the increase in 
                    % distance from dUsers(u1) to the real distance traver-
                    % sed by the signal in on each path
                    
                    % Pathloss into consideration (old)
                    PRx(u1) = PRx(u1) + problem.alphaChannels(u1,i)*...
                              pattern(handle_Conf_Array,problem.freq,...
                                      problem.phiChannels(u1,i),...
                                      problem.thetaChannels(u1,i),...
                                      'Type','Power','Normalize',false)*...
                                      (problem.lambda/(4*pi*problem.dUsers(u1)))^2;
                                  
                    % Consider only antenna directivity (seems more fair)
                    PRx(u1) = PRx(u1) + problem.alphaChannels(u1,i)*...
                              patternAzimuth(handle_Conf_Array, problem.freq, ...
                                             problem.thetaUsers(u), ...
                                             'Azimuth', problem.phiUsers(u),...
                                             'Type', 'power');
                        
                end
            end
        end
        if conf.verbosity > 1
            fprintf('PRx: %f\n',PRx(problem.IDUserAssigned));
            if PRx(problem.IDUserAssigned) == 0 && Position_Taper(end) ~= 0
                fprintf('Chachi!');
            end
        end
        PR = PRx(problem.IDUserAssigned);
        Int = sum(PRx)-PRx(problem.IDUserAssigned);
        Int = Int/(numel(PRx)-1);
    end
    
    % Consider beamwidth into opt. (efficiency metric)
    if conf.Fweights(2) > 0
        % Define variable range
        azimuth = -180:2:180;
        elevation = -90:2:90;
        % Compute bandwidth azymuth
        bw_az = o_beamwidth(handle_Conf_Array,problem.freq,...
            azimuth,problem.thetaChannels(problem.IDUserAssigned),3,...
            conf.PlotAssignments);
        % Compute bandwidth elevation
        bw_el = o_beamwidth(handle_Conf_Array,problem.freq,...
            problem.phiChannels(problem.IDUserAssigned),elevation,3,...
            conf.PlotAssignments);
        % Compute final bandwidth
        HPBW = bw_az * bw_el;
        if conf.verbosity > 1
            fprintf('The HPBW (az) computed by the helper function is %f\n',bw_az);
            fprintf('The HPBW (el) computed by the helper function is %f\n',bw_el);
        end
    end

    % Computing the function value from the SL_Mag
    Val = ( conf.Fweights(1)*(1/PR) * conf.Fweights(2)*(Int) ) + conf.Fweights(3)*(HPBW/(360^2));
        
    if conf.verbosity > 1
        fprintf('Val = %f\n',Val);
    end
end
