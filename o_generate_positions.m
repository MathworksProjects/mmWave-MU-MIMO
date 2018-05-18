function [thetaPos, phiPos, dPos] = o_generate_positions(conf,problem)
    % Generate Normal distribution
    pd = makedist('Normal');
    pd.sigma = 45;
    t = truncate(pd,-45,45);

    % Check if users were assigned Theta angles (deterministic)
    if ~conf.detLocation || ~isfield(problem,'thetaUsers')
        % Generate Elevation angles
        thetaPos = random(t,1,problem.nUsers);
        if conf.verbosity >= 1
            fprintf('Users were not assigned Theta values (positions) in the space.\n');
        end
    elseif problem.anglesInRadians
        % Convert Theta from radians to degrees
        thetaPos = problem.thetaUsers(1:problem.nUsers)/(2*pi)*360;
    else
        % Deterministic Theta angles in degrees
        thetaPos = problem.thetaUsers(1:problem.nUsers);
    end

    % Check if users were assigned Phi angles (deterministic)
    if ~conf.detLocation || ~isfield(problem,'phiUsers')
        % Generate Azymuth angles
        phiPos = random(t,1,problem.nUsers);
        if conf.verbosity >= 1
            fprintf('Users were not assigned Theta values (positions) in the space.\n');
        end
    elseif problem.anglesInRadians
        % Convert Phi from radians to degrees
        phiPos = problem.phiUsers(1:problem.nUsers)/(2*pi)*360;
    else
        % Deterministic Phi angles in degrees
        phiPos = problem.phiUsers(1:problem.nUsers);
    end

    % Check if users were assigned Distances (deterministic)
    if ~conf.detLocation || ~isfield(problem,'dUsers')
        dPos = rand(1,problem.nUsers) * (problem.maxdUsers-problem.mindUsers) + problem.mindUsers;
    else
        dPos = problem.dUsers(1:problem.nUsers);
    end

    % Display user distribution in space
    if conf.verbosity >= 1
        fprintf('Elevations assigned:\n');
        display(thetaPos);
        fprintf('Azimuths assigned:\n');
        display(phiPos);
        fprintf('Distances to BS assigned:\n');
        display(dPos);
    end 
end