function o_plotAssignment_mod(problem, handle_Conf_Array)
%PLOTASSIGNMENT Plot a given Assignment 
    % Array geometry and pattern for optimized values.
    f1 = figure;
    handle_Conf_Array.viewArray;
    f2 = figure;
    handle_Conf_Array.pattern(problem.freq,'Type','powerdb');
    %%
    hold on
    dx = 0.1;
    dy = 0.1;
    dz = 0.3;
    [x_u,y_u,z_u] = sph2cart(problem.phiUsers/360*2*pi,...
        problem.thetaUsers/360*2*pi,...
        3+50*ones(size(problem.phiUsers)));%3*10^(max(max(p))/10)*
    for i = 1:problem.nUsers
        line([0,x_u(i)],[0,y_u(i)],[0,z_u(i)]);
        if i == problem.IDUserAssigned
            scatter3(x_u(i),y_u(i),z_u(i),'r');
        else
            scatter3(x_u(i),y_u(i),z_u(i),'g');
        end
        text(x_u(i)+dx, y_u(i)+dy, z_u(i)+dz, ['User ',num2str(i)]);
    end
    %%
    azimuth = -180:2:180;
    elevation = -90:2:90;
    
    f3 = figure;
    % Plot designated device
    subplot(2,problem.nUsers,1);
    pattAz = patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(problem.IDUserAssigned),'Azimuth',azimuth,'Type','powerdb'); 
    polarpattern(circshift(pattAz,90));
%     tt = strcat('User',{' '},num2str(problem.IDUserAssigned),{' '},'Phi:',{' '},num2str(problem.phiUsers(problem.IDUserAssigned)));
%     sp11.title.set_text(tt);
    % Plot interfeered devices
    intVect = 1:1:problem.nUsers;
    intVect(intVect==problem.IDUserAssigned) = [];
    for idx = 1:length(intVect)
        subplot(2,problem.nUsers,idx+1)
        theta = problem.thetaUsers(intVect(idx));
        pattAz = patternAzimuth(handle_Conf_Array,problem.freq,theta,'Azimuth',azimuth,'Type','powerdb');
        polarpattern(circshift(pattAz,90));
%         tt = strcat('User',{' '},num2str(intVect(idx)),{' '},'Phi:',{' '},num2str(problem.phiUsers(intVect(idx))));
%         title(tt,'FontSize',14);
    end

    % Plot designated device    
    subplot(2,problem.nUsers,problem.nUsers+1)
    pattEl = patternElevation(handle_Conf_Array,problem.freq,problem.phiUsers(problem.IDUserAssigned),'Elevation',elevation,'Type','powerdb');
    polarpattern(circshift(pattEl,90));
%     tt = strcat('User',{' '},num2str(problem.IDUserAssigned),{' '},'Theta:',{' '},num2str(problem.thetaUsers(problem.IDUserAssigned)));
%     title(tt,'FontSize',14);
    % Plot interfeered devices
    intVect = 1:1:problem.nUsers;
    intVect(intVect==problem.IDUserAssigned) = [];
    for idx = 1:length(intVect)
        subplot(2,problem.nUsers,idx+problem.nUsers+1)
        phi = problem.thetaUsers(intVect(idx));
        pattEl = patternElevation(handle_Conf_Array,problem.freq,phi,'Elevation',elevation,'Type','powerdb');
        polarpattern(circshift(pattEl,90));
%         tt = strcat('User',{' '},num2str(intVect(idx)),{' '},'Theta:',{' '},num2str(problem.thetaUsers(intVect(idx))));
%         title(tt,'FontSize',14);
    end
%     suptt = strcat('Radiation pattern for user',{' '},num2str(problem.IDUserAssigned));
%     suptitle(suptt,'FontSize',14);
end

