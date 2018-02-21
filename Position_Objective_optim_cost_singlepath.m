function [Val] = Position_Objective_optim_cost_singlepath(Position_Taper, conf, problem)
% Copyright 2017  The MathWorks, Inc.
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

    % Setting figures to invisible
    %set(groot,'defaultFigureVisible','off')

    azimuth = -180:2:180;
    elevation = -90:2:90;

    % Creating a cosine Antenna Element
    handle_Ant =phased.CosineAntennaElement('FrequencyRange',...
        [problem.freq-(problem.Bw/2) problem.freq+(problem.Bw/2)],...
        'CosinePower',[1.5 2.5]);

    [Position_Taper_elem,problem] = subarrayToAntennaElements(Position_Taper,conf,problem);
    % Extracting taper values from the input vector
    Taper_value = Position_Taper_elem(problem.ant_elem+1:problem.ant_elem*2) .* ...
        exp(1i.*Position_Taper_elem(problem.ant_elem*2+1:problem.ant_elem*3));

    % Creating a Conformal Array with cosine elements.
    % Conformal array will be limited to a single plane

%     disp('Positions');
%     [zeros(1,problem.ant_elem);...
%         problem.possible_locations(2,Position_Taper_elem(1:problem.ant_elem));...
%         problem.possible_locations(3,Position_Taper_elem(1:problem.ant_elem))]
%     disp('Taper value');
    if isempty(Taper_value)
        disp('uyuyuy');
    end
    
    handle_Conf_Array = phased.ConformalArray('Element',handle_Ant,...
        'ElementPosition',[zeros(1,problem.ant_elem);...
        problem.possible_locations(2,Position_Taper_elem(1:problem.ant_elem));...
        problem.possible_locations(3,Position_Taper_elem(1:problem.ant_elem))],...
        'Taper',Taper_value);

    % Extracting Azimuth and Elevation Patterns
    [fieldval_az] = patternAzimuth(handle_Conf_Array,problem.freq,...
        problem.thetaUsers(problem.IDUserAssigned),'Azimuth',azimuth,'Type','Power');
    [fieldval_el] = patternElevation(handle_Conf_Array,problem.freq,...
        problem.phiUsers(problem.IDUserAssigned),'Elevation',elevation,'Type','Power');

    % Replacing any inf and zeros with a very low value
    fieldval_az(isinf(fieldval_az)) = -400;
    fieldval_az((fieldval_az == 0)) = -400;

    % Genarating a polar plot with the extracted data
    % the polar plot method "FindLobes" is used to fetch the required 
    % information about the beam.
    if conf.findLobes
        figure(1);
        subplot(1,2,2);
        %f3 = figure('Visible','off');
        handle_polar_az = polarpattern(azimuth,mag2db(fieldval_az));
        handle_polar_az.AntennaMetrics = 1;
        %close(gcf)
        lobe_info = findLobes(handle_polar_az);
    end
    
    [bw_az,SLL_az] = helperFindLobes(handle_Conf_Array,problem.freq,azimuth,problem.thetaUsers(problem.IDUserAssigned),3,conf.PlotAssignments);
    if conf.DEBUG
        fprintf('The HPBW (az) computed by the helper function is %f\n',bw_az);
        fprintf('The SLL (az) computed by the helper function is %f\n',SLL_az);
    end
        
    if conf.findLobes
        HPBW_az = lobe_info.HPBW;%sideLobes.magnitude;
        SLL_az = lobe_info.SLL;
        if conf.DEBUG
            fprintf('The HPBW (az) computed by the findLobes function is %f\n',HPBW_az);
            fprintf('The SLL (az) computed by the findLobes function is %f\n',SLL_az);
        end
    else
        HPBW_az = bw_az;
    end

    if(isempty(HPBW_az))

        HPBW_az = 180;

    end

    % Replacing any inf and zeros with a very low value
    fieldval_el(isinf(fieldval_el)) = -400;
    fieldval_el((fieldval_el == 0)) = -400;

    % Fetching the Sidelobe Magnitude 
    % Magnitude of the elevation beam pattern
    if conf.findLobes
        figure(1);
        subplot(1,2,1);
        %f4 = figure('Visible','off');
        handle_polar_el = polarpattern(elevation,mag2db(fieldval_el));
        handle_polar_el.AntennaMetrics = 1;
        %close(gcf)
        lobe_info = findLobes(handle_polar_el);
    end
    
    [bw_el, SLL_el] = helperFindLobes(handle_Conf_Array,problem.freq,problem.phiUsers(problem.IDUserAssigned),elevation,3,conf.PlotAssignments);
    if conf.DEBUG
        fprintf('The HPBW (el) computed by the helper function is %f\n',bw_el);
        fprintf('The SLL (el) computed by the helper function is %f\n',SLL_el);
    end
    
    if conf.findLobes
        HPBW_el = lobe_info.HPBW;%sideLobes.magnitude;
        SLL_el = lobe_info.SLL;
        if conf.DEBUG
            fprintf('The HPBW (el) computed by the findLobes function is %f\n',HPBW_el);
            fprintf('The SLL (el) computed by the findLobes function is %f\n',SLL_el);
        end
    else
        HPBW_el = bw_el;
    end
    
    %Using a very high value if the value are empty

    if(isempty(HPBW_el))

        HPBW_el = 360;

    end

    % Summing the azimuth and elevation values
    % Adding term to minimize the difference between the Az and El patterns

    HPBW = HPBW_az * HPBW_el;
    SLL = min(SLL_az, SLL_el);
    
    % Computing the received power in every user
    PRx = zeros(1,problem.nUsers);
    for u=1:problem.nUsers
        PRx(u) = pattern(handle_Conf_Array,problem.freq,...
                problem.phiUsers(u),...
                problem.thetaUsers(u),...
                'Type','Power');
    end

    PR = PRx(problem.IDUserAssigned);
    Int = sum(PRx)-PRx(problem.IDUserAssigned);
    Int = Int/(numel(PRx)-1);
    
    % Computing directivity
    %D = directivity(handle_Conf_Array,conf.freq,[0;0]);

    % Computing the function value from the SL_Mag 
    Val = Int/PR + HPBW/(360^2) + 1/SLL; %+ 1/10^(D/10)
    if conf.DEBUG
        fprintf('Val = %f\n',Val);
    end
end
