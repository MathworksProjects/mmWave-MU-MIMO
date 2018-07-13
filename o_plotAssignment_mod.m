function o_plotAssignment_mod(problem, handle_Conf_Array)
%PLOTASSIGNMENT Plot a given Assignment 
    % Preprocess data
    IDmax = problem.IDUserAssigned;  % Desired user
    IDmin = 1:1:problem.nUsers;
    IDmin(IDmin==problem.IDUserAssigned) = [];  % Users interfereed
    azimuth = -180:1:180;  % Azymuth angle range
    elevation = -90:1:90;  % Elevation angle range
    
    % Array geometry and pattern for optimized values.
%     figure;
%     handle_Conf_Array.viewArray;
%     figure;
%     handle_Conf_Array.pattern(problem.freq,'Type','powerdb');
%     hold on
%     dx = 0.1;   dy = 0.1;   dz = 0.3;
%     [x_u,y_u,z_u] = sph2cart(problem.phiUsers/360*2*pi,...
%         problem.thetaUsers/360*2*pi,...
%         3+50*ones(size(problem.phiUsers)));%3*10^(max(max(p))/10)*
%     % Visualize user location/channels
%     for i = 1:problem.nUsers
%         line([0,x_u(i)],[0,y_u(i)],[0,z_u(i)]);
%         if i == IDmax;   scatter3(x_u(i),y_u(i),z_u(i),'r');
%         else;            scatter3(x_u(i),y_u(i),z_u(i),'g');
%         end
%         text(x_u(i)+dx, y_u(i)+dy, z_u(i)+dz, ['User ',num2str(i)]);
%     end
    fprintf('Plotting selected pattern for user %d\n',IDmax);
    fprintf('* %d [Theta: %.1f - Phi: %.1f]\n',IDmax,problem.thetaUsers(IDmax),problem.phiUsers(IDmax));
    for idx = 1:length(IDmin)
        fprintf('  %d [Theta: %.1f - Phi: %.1f]\n',IDmin(idx),problem.thetaUsers(IDmin(idx)),problem.phiUsers(IDmin(idx)));
    end
    
    figure;
    % Plot Theta cuts - specify Phi
    % Plot designated device
    subplot(2,problem.nUsers,1);
    tt = strcat('User',{' '},num2str(IDmax),{' '},'Phi:',{' '},num2str(problem.phiUsers(IDmax)));
    pattAz = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(IDmax),'Azimuth',azimuth,'Type','powerdb'); 
    polarpattern(azimuth,pattAz,'MagnitudeLim',max(pattAz)+[-50 0],'AntennaMetrics',0,'TitleTop',tt{:},'TitleTopFontSizeMultiplier',2);
    hold on;
    mark = personalMarker(azimuth,pattAz,problem.phiUsers(IDmax));
    polarpattern(azimuth,mark);
    % Plot interfeered devices
    for idx = 1:length(IDmin)
        tt = strcat('User',{' '},num2str(IDmin(idx)),{' '},'Phi:',{' '},num2str(problem.phiUsers(IDmin(idx))));
        subplot(2,problem.nUsers,idx+1);
        theta = problem.thetaUsers(IDmin(idx));
        pattAz = patternAzimuth(handle_Conf_Array,problem.freq,theta,'Azimuth',azimuth,'Type','powerdb');
        polarpattern(azimuth,pattAz,'MagnitudeLim',max(pattAz)+[-50 0],'AntennaMetrics',0,'TitleTop',tt{:},'TitleTopFontSizeMultiplier',2);
        hold on;
        mark = personalMarker(azimuth,pattAz,problem.phiUsers(IDmin(idx)));
        polarpattern(azimuth,mark);
    end
    
    % Plot Phi cuts - specify Theta
    % Plot designated device
    subplot(2,problem.nUsers,problem.nUsers+1)
    tt = strcat('User',{' '},num2str(IDmax),{' '},'Theta:',{' '},num2str(problem.thetaUsers(IDmax)));
    pattEl = patternElevation(handle_Conf_Array,problem.freq,problem.phiUsers(IDmax),'Elevation',elevation,'Type','powerdb');
    polarpattern(elevation,pattEl,'MagnitudeLim',max(pattEl)+[-50 0],'AntennaMetrics',0,'TitleTop',tt{:},'TitleTopFontSizeMultiplier',2);
    hold on;
    mark = personalMarker(elevation,pattEl,problem.thetaUsers(IDmax));
    polarpattern(elevation,mark);
    % Plot interfeered devices
    for idx = 1:length(IDmin)
        subplot(2,problem.nUsers,problem.nUsers+idx+1);
        tt = strcat('User',{' '},num2str(IDmin(idx)),{' '},'Theta:',{' '},num2str(problem.thetaUsers(IDmin(idx))));
        phi = problem.thetaUsers(IDmin(idx));
        pattEl = patternElevation(handle_Conf_Array,problem.freq,phi,'Elevation',elevation,'Type','powerdb');
        polarpattern(elevation,pattEl,'MagnitudeLim',max(pattEl)+[-50 0],'AntennaMetrics',0,'TitleTop',tt{:},'TitleTopFontSizeMultiplier',2);
        hold on;
        mark = personalMarker(elevation,pattEl,problem.thetaUsers(IDmin(idx)));
        polarpattern(elevation,mark);
    end
    
    suptt = strcat('Pattern for user',{' '},num2str(IDmax),{' '},'(Nant =',{' '},num2str(problem.N_Antennas),')');
    tit = suptitle(suptt{:});
    set(tit,'FontSize',14)
end

function mark = personalMarker(angleRange,directivity,userAngle)
    % Creates a marker to visualize the user's angle (LoS) or angles
    % (multipath/multipath-5G)
    mark = (-1).*Inf(size(angleRange));
    [~,tIdx] = min(abs(angleRange-userAngle));
    mark(tIdx) = max(directivity);
end

