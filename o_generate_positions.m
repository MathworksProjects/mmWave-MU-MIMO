function [thetaPos, phiPos, dPos] = o_generate_positions(conf,nUsers,...
                                                        maxdUsers,mindUsers)
    % Generate Normal distribution
    pd = makedist('Normal');
    pd.sigma = 45;
    t = truncate(pd,-45,45);

    % Check if users were assigned Theta angles (deterministic)
    if ~conf.detLocation || ~isfield(conf,'thetaUsers')
        % Generate Elevation angles
        thetaPos = random(t,1,nUsers);
        if conf.verbosity >= 1
            fprintf('Users were not assigned Theta values (positions) in the space.\n');
        end
    elseif conf.anglesInRadians
        % Convert Theta from radians to degrees
        thetaPos = conf.thetaUsers(1:nUsers)/(2*pi)*360;
    else
        % Deterministic Theta angles in degrees
        thetaPos = conf.thetaUsers(1:nUsers);
    end

    % Check if users were assigned Phi angles (deterministic)
    if ~conf.detLocation || ~isfield(conf,'phiUsers')
        % Generate Azymuth angles
        phiPos = random(t,1,nUsers);
        if conf.verbosity >= 1
            fprintf('Users were not assigned Theta values (positions) in the space.\n');
        end
    elseif conf.anglesInRadians
        % Convert Phi from radians to degrees
        phiPos = conf.phiUsers(1:nUsers)/(2*pi)*360;
    else
        % Deterministic Phi angles in degrees
        phiPos = conf.phiUsers(1:nUsers);
    end

    % Check if users were assigned Distances (deterministic)
    if ~conf.detLocation || ~isfield(conf,'dUsers')
        dPos = rand(1,nUsers) * (maxdUsers-mindUsers) + mindUsers;
    else
        dPos = conf.dUsers(1:nUsers);
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