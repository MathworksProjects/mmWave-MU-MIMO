function o_plotAssignment(problem, handle_Conf_Array)
%PLOTASSIGNMENT Plot a given Assignment 
    % Array geometry and pattern for optimized values.
    f1 = figure;
    handle_Conf_Array.viewArray;
    f2 = figure;
    p = handle_Conf_Array.pattern(problem.freq,'Type','powerdb');
    
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
    f3 = figure;
    azimuth = -180:2:180;
    elevation = -90:2:90;
    patternAzimuth(handle_Conf_Array,problem.freq,problem.thetaUsers(problem.IDUserAssigned),'Azimuth',azimuth,'Type','powerdb');
    
    f4 = figure;
    patternElevation(handle_Conf_Array,problem.freq,problem.phiUsers(problem.IDUserAssigned),'Elevation',elevation,'Type','powerdb');
    %handle_Conf_Array.patternElevation(problem.freq,'Type','powerdb');
    disp('Program paused, press any key to continue')
    pause
    if isvalid(f1)
        close(f1)
    end
    if isvalid(f2)
        close(f2)
    end
    if isvalid(f3)
        close(f3)
    end
    if isvalid(f4)
        close(f4)
    end
end

