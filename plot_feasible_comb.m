% ----------------------------------------------------------------------- %
% ----------------------------- PLOTTING -------------------------------- %
% ----------------------------------------------------------------------- %

% Plot Input Arrangements (as computed by Algorithm 1)
% Get the dimensions for the subplot
for j = 1:problem.nUsers; if (j^2)>=problem.nUsers; sPtDim = j; break; end; end
for n = 1:problem.nUsers
    figure(conf.figIdx);
    subplot(sPtDim,sPtDim,n); hold on;
    plotPatch(patch,conf.colorEmpty,50);
    plotArray(arrays(:,:,n),'r',70);
end
conf.figIdx = conf.figIdx + 1;

figure(conf.figIdx); hold on
plotPatch(patch,conf.colorEmpty,50);
for n = 1:problem.nUsers
    plotArray(arrays(:,:,n),conf.colorList{n},70);
end
conf.figIdx = conf.figIdx + 1;