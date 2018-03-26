function [Val] = o_Position_Objective_optim_cost_multipath(Position_Taper, conf, problem)
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

%     azimuth = -180:2:180;
%     elevation = -90:2:90;

    [Position_Taper_elem,problem] = o_subarrayToAntennaElements(Position_Taper,conf,problem);
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
    
    handle_Conf_Array = phased.ConformalArray('Element',problem.handle_Ant,...
        'ElementPosition',[zeros(1,problem.ant_elem);...
        problem.possible_locations(2,Position_Taper_elem(1:problem.ant_elem));...
        problem.possible_locations(3,Position_Taper_elem(1:problem.ant_elem))],...
        'Taper',Taper_value);
    
%     bw_az = -Inf;
%     for i=1:problem.maxnChannelPaths
%         if problem.thetaChannels(problem.IDUserAssigned,i) ~= -Inf
%             bw_az = max(bw_az,beamwidth(handle_Conf_Array,problem.freq,...
%                 azimuth,problem.thetaChannels(problem.IDUserAssigned,i),3,...
%                 conf.PlotAssignments));
%         end
%     end
%     if conf.DEBUG
%         fprintf('The HPBW (az) computed by the helper function is %f\n',bw_az);
%     end
%     bw_el = -Inf;
%     for i=1:problem.maxnChannelPaths
%         if problem.phiChannels(problem.IDUserAssigned,i) ~= -Inf
%             bw_el = max(bw_el,beamwidth(handle_Conf_Array,problem.freq,...
%                 problem.phiChannels(problem.IDUserAssigned,i),elevation,3,...
%                 conf.PlotAssignments));
%         end
%     end
%     if conf.DEBUG
%         fprintf('The HPBW (el) computed by the helper function is %f\n',bw_el);
%     end
%     
% 
%     % Summing the azimuth and elevation values
%     % Adding term to minimize the difference between the Az and El patterns
% 
%     HPBW = bw_az * bw_el;
    
    % Computing the received power in every user
    PRx = zeros(1,problem.nUsers);
    for u1=1:problem.nUsers
        for i=1:problem.maxnChannelPaths
            if problem.phiChannels(u1,i) ~= -Inf
                PRx(u1) = PRx(u1) + problem.alphaChannels(u1,i)*...
                    pattern(handle_Conf_Array,problem.freq,...
                    problem.phiChannels(u1,i),...
                    problem.thetaChannels(u1,i),...
                    'Type','Power','Normalize',false);
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
    
    % Computing directivity
    %D = directivity(handle_Conf_Array,conf.freq,[0;0]);

    % Computing the function value from the SL_Mag
    if problem.nUsers > 1
        Val = Int/PR;% + HPBW/(360^2); %+ 1/10^(D/10)
    else
        Val = 1/PR;
    end
    if conf.verbosity > 1
        fprintf('Val = %f\n',Val);
    end
end
