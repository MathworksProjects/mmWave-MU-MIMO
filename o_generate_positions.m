function [thetaPos, phiPos, dPos] = o_generate_positions(conf,nUsers,...
            maxdUsers,mindUsers)
%% Randomization and conversion of the users' positions in the space
    %%%% Only in the case where no determinist position is given in the input
    %%%% file defining the problem. The conversion to degrees is performed if
    %%%% angles were given in radians
    pd = makedist('Normal');
    pd.sigma = 45;
    t = truncate(pd,-45,45);
    % if ~isfield(problem,'thetaUsers') || isempty(problem.thetaUsers)
    if conf.verbosity >= 1
        fprintf('===================================================\n');
        fprintf('The users were not assigned positions in the space:\n');
        fprintf('Assigning random values...\n');
    end
    thetaPos = random(t,1,nUsers);
    if conf.verbosity >= 1
        fprintf('New elevations assigned:\n');
        display(thetaPos);
    end
%     elseif conf.anglesInRadians
%         thetaPos = problem.thetaUsers/(2*pi)*360;
%     else % Nothing to do
%         thetaPos = problem.thetaUsers;
%     end
    % if ~isfield(problem,'phiUsers') || isempty(problem.phiUsers)
    phiPos = random(t,1,nUsers);
    if conf.verbosity >= 1
        fprintf('New azimuths assigned:\n');
        display(phiPos);
    end
%     elseif conf.anglesInRadians
%         phiPos = problem.phiUsers/(2*pi)*360;
%     else % Nothing to do
%         phiPos = problem.phiUsers;
%     end
    % if ~isfield(problem,'phiUsers') || isempty(problem.phiUsers)
    dPos = rand(1,10) * (maxdUsers-mindUsers) + mindUsers;
%     else % Nothing to do
%         dPos = problem.dUsers;
%     end
end

