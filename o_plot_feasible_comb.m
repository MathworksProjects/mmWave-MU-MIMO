function o_plot_feasible_comb(problem,conf,patch,arrays)
    % Plot Input Arrangements (as computed by Algorithm 1)
    % Get the dimensions for the subplot
    for j = 1:problem.nUsers; if (j^2)>=problem.nUsers; sPtDim = j; break; end; end
    % Get real number of users
    nUsers = size(arrays,3);
    % Plot antenna assignation
    figure;
    for n = 1:problem.nUsers
        subplot(sPtDim,sPtDim,n); hold on;
        o_plotArray(patch,conf.colorEmpty,50,true);
        if n<=nUsers
            o_plotArray(arrays(:,:,n),'r',200,false);
        end
    end
    % Plot antenna assignation (altogether)
    figure; hold on
    o_plotArray(patch,conf.colorEmpty,50,true);
    for n = 1:problem.nUsers
        if n<=nUsers
            o_plotArray(arrays(:,:,n),conf.colorList{n},1600,false);
        end
    end
end