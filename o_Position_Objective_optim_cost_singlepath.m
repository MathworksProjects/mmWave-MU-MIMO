function [Val] = o_Position_Objective_optim_cost_singlepath(Position_Taper, conf, problem)
    
    % Initialize optimization variables
    PR   = 0;  % Received power at intended user
    Int  = 0;  % Interference generated to other users
    HPBW = 0;  % Beamdiwth achieved
    
    % Define constants
    if ~strcmp(conf.genStructure, 'allAntennas')
        nck = strcmp(conf.genStructure, 'nchoosek');
        if sum(Position_Taper(1:problem.Nmax-nck*(problem.Nmax-1)) ~= ...
                floor(Position_Taper(1:problem.Nmax-nck*(problem.Nmax-1)))) > 0
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
        for u=1:problem.nUsers
            % Pathloss into consideration (old)
%             PRx(u) = pattern(handle_Conf_Array,problem.freq,...
%                     problem.phiUsers(u),...
%                     problem.thetaUsers(u),...
%                     'Type','Power')*...
%                     (problem.lambda/(4*pi*problem.dUsers(u)))^2;
            % Consider only antenna directivity (seems more fair)
            PRx(u) = patternAzimuth(handle_Conf_Array, problem.freq, ...
                                    problem.thetaUsers(u), ...
                                    'Azimuth', problem.phiUsers(u),...
                                    'Type', 'power');
        end

        PR = PRx(problem.IDUserAssigned);  % Linear
        Int = sum(PRx)-PRx(problem.IDUserAssigned);
        Int = Int/(numel(PRx)-1);    % Linear
    end
    % Compute values in dB
    PRdB  = (10*log10(PR));  % in DB
	IntdB = (10*log10(Int));  % in DB
    
    % Consider beamwidth into opt. (efficiency metric)
    if conf.Fweights(2) > 0
        % Define variable range
        azimuth = -180:2:180;
        elevation = -90:2:90;
        % Compute bandwidth azymuth
        bw_az = o_beamwidth(handle_Conf_Array,problem.freq,...
            azimuth,problem.thetaUsers(problem.IDUserAssigned),3,...
            conf.PlotAssignments);
        % Compute bandwidth elevation
        bw_el = o_beamwidth(handle_Conf_Array,problem.freq,...
            problem.phiUsers(problem.IDUserAssigned),elevation,3,...
            conf.PlotAssignments);
        % Compute final beamwidth
        HPBW = bw_az * bw_el;
        if conf.verbosity > 1
            fprintf('The HPBW (az) is %f\n',bw_az);
            fprintf('The HPBW (el) is %f\n',bw_el);
        end
    end    

    % Computing the function value from the SL_Mag
    
%     % Opt 1
%     Val = ( conf.Fweights(1)*(1/PR) * conf.Fweights(2)*(Int) ) + conf.Fweights(3)*(HPBW/(360^2));
%     % Opt 2
%     corrPR = transferFunction(PRdB);
%     corrInt = transferFunction(IntdB);
%     Val = abs ( ( conf.Fweights(1)*corrPR*(1/PR) * conf.Fweights(2)*corrInt*(Int) ) + conf.Fweights(3)*(HPBW/(360^2)) );
%     % Opt 3
%     WeightIDmax = (1/problem.nUsers);
%     WeightIDmin = (1 - 1/problem.nUsers);
%     Val = ( WeightIDmax*(1/PR) * WeightIDmin*(Int) ) + conf.Fweights(3)*(HPBW/(360^2));
%     % Opt 4
%     WeightIDmax = (-1)*(1/problem.nUsers);
%     WeightIDmin = (+1)*(1 - 1/problem.nUsers);
%     Val = WeightIDmax*abs(PRdB) + WeightIDmin*abs(IntdB) + conf.Fweights(3)*(HPBW/(360^2));
    % Opt 5
    scorePrx = transferScore(PRdB,1);
    scoreInt = transferScore(IntdB,0);
    Val = (-1)*( conf.Fweights(1)*scorePrx + conf.Fweights(2)*scoreInt );

    if conf.verbosity > 1
        fprintf('Val = %f\n',Val);
        fprintf('scorePrx = %.2f\tscoreInt = %.2f\n',scorePrx,scoreInt);
        fprintf('Dir(dB)= %.2f\tDir_int(dB)= %.2f\tHPBW= %.2f\n',PRdB,IntdB,(HPBW/(360^2)));
    end
end

function correction = transferFunction(PowerdB)                        %#ok
    % Interpolate
    x = [-500 -30 -20 -15 -10 -5 0 5 10 15 +20 +30 +500];  % Range
    y = [1 1 0.75 0.5 0.25 0.1 0.1 0.1 0.25 0.5 0.75 1 1];  % Values
    xq = -500:0.1:500;  % Interpolate Range
    vq1 = interp1(x,y,xq);    % Interpolate Values
    [~,idx] = min(abs(PowerdB - xq));
    correction = vq1(idx);
end

function score = transferScore(PowerdB,mode)                           
    x = [-500 -5 0 5 10 15 20 25 500];
    if mode
        % Prx (intended user)
        y = [0 0 0.3 0.5 0.7 0.9 1 1 1];
    else
        % Int (interfereed users)
        y = 1 - [0 0 0.3 0.5 0.7 0.9 1 1 1];
    end
    xx = -500:0.05:500;
    yy = interpn(x,y,xx,'pchip');
    [~,idx] = min(abs(PowerdB - xx));
    score = yy(idx);
end
