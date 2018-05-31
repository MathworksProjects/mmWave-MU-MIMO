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
    % TRANSFERSCORE - Computes a score that reflects how good the
    % directivity value is in the MU-MIMO scenario. Directivities start
    % getting an score when their absolute value is greater than 0.
    % Directivities are capped up at the maximum tolerable antenna gain by
    % the FCC community (33dB if Ptx=10dBm)
    %
    % For more information, see slide 7-4 in:
    % https://www.cse.wustl.edu/~jain/cse574-14/ftp/j_07sgh.pdf
    %
    % Syntax:  score = transferScore(PowerdB,mode)
    %
    % Inputs:
    %    PowerdB - Description
    %    mode - 1 for Directivity to intended user. 0 to compute
    %    directivity towards interfeered users.
    %
    % Outputs:
    %    score - Value between 0 and 1 that reflects how good the current
    %    antenna selection is.
    %
    % Example: 
    %    score = transferScore(10,1)  % Compute score for intended user
    %    score = transferScore(-10,0)  % Compute score to other users
    %
    %------------- BEGIN CODE --------------
    
    % Define the scoring function here. Left hand-side values reflect the
    % power levels in dB of the directivity. The right hand-side values
    % reflect the score obtained
    p = [ -500 0     1;
           -50 0     1;
           -40 0   0.9;
           -33 0  0.85;
           -30 0   0.8;
           -25 0   0.7;
           -20 0   0.6;
           -15 0   0.4;
           -10 0   0.3;
            -5 0   0.2;
             0 0     0;
             5 0.2   0;
            10 0.3   0;
            15 0.4   0;
            20 0.6   0;
            25 0.8   0;
            30 0.9   0;
            33 1     0;
            40 1     0;
            50 1     0;
           500 1     0].';
	x = p(1,:);  % Power in dB
    if mode
        % Prx (intended user)
        y = p(2,:);  % Score
    else
        % Int (interfereed users)
        y = p(3,:);  % Score
    end
    xx = min(x):0.05:max(x);
    yy = interpn(x,y,xx,'linear');
    [~,idx] = min(abs(PowerdB - xx));
    score = yy(idx);

    % Alternative (linear)
%     if mode    
%         mindB = 0;
%         maxdB = 33;
%         minScore = 0;
%         maxScore = 1;
%         m = (maxScore-minScore)/(maxdB-mindB);  % Slope 
%         x = (mindB:0.01:maxdB);
%         y = x.*m + minScore;
%         % Include limits
%         x = [-500 x 500];
%         y = [0 y 1];
%     else
%         % 1st stage
%         mindB1st = 0;
%         maxdB1st = -50;
%         minScore1st = 0;
%         maxScore1st = 0.85;
%         m = (maxScore1st-minScore1st)/(maxdB1st-mindB1st);  % Slope 
%         x1 = (maxdB1st:0.01:mindB1st);
%         n1 = 0;
%         y1 = x1.*m + n1;
%         % 2nd stage
%         mindB2nd = maxdB1st;
%         maxdB2nd = -100;
%         minScore2nd = maxScore1st;
%         maxScore2nd = 1;
%         m = (maxScore2nd-minScore2nd)/(maxdB2nd-mindB2nd);  % Slope 
%         x2 = (maxdB2nd:0.01:mindB2nd);
%         n2 = maxScore2nd - m*maxdB2nd;
%         y2 = x2.*m + n2;
%         % Append stages
%         x = [x2 x1];
%         y = [y2 y1];
%         % Include limits
%         x = [-500 x 500];
%         y = [1 y 0];
%     end
%     [~,idx] = min(abs(PowerdB - x));
%     score = y(idx);
end
