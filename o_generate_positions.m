function [thetaPos, phiPos, dPos] = o_generate_positions(conf,problem)
    % Generate Normal distribution
    pd = makedist('Normal');
    pd.sigma = 45;
    t = truncate(pd,-45,45);
    % Check for consistency in number of users
    if problem.nUsers~=2 && conf.useCasesLocation
        error('When User Location are set to UserCase, the number of users should be 2.');
    end
    % Generate Angles for usee case (uc) distribution
    % UC 1
    uc_el(1,:) = [0 0];      uc_az(1,:) = [+15 -15];  uc_dist(1,:) = [5 5];
    % UC 2
    uc_el(2,:) = [0 0];      uc_az(2,:) = [+30 -30];  uc_dist(2,:) = [5 5];
    % UC 3
    uc_el(3,:) = [+15 -15];  uc_az(3,:) = [0 0];      uc_dist(3,:) = [5 5];
    % UC 4
    uc_el(4,:) = [+15 -15];  uc_az(4,:) = [15 -15];   uc_dist(4,:) = [5 5];
    % UC 5
    uc_el(5,:) = [+15 -15];  uc_az(5,:) = [15 15];    uc_dist(5,:) = [5 5];
    % UC 6
    uc_el(6,:) = [+15 +15];  uc_az(6,:) = [15 -15];   uc_dist(6,:) = [5 5];

    % Check if users were assigned Theta angles (deterministic)
    if ~conf.detLocation || ~isfield(problem,'thetaUsers')
        % Generate Random Elevation angles
        thetaPos = random(t,1,problem.nUsers);
        if conf.verbosity >= 1
            fprintf('Users were not assigned Theta values (positions) in the space.\n');
        end
    elseif conf.useCasesLocation
        % Generate Elevation angles from Use case
        thetaPos = uc_el(conf.useCaseLoation,:);
    else
        % Retrieve Elevation angles from config file
        if problem.anglesInRadians
            % Convert Theta from radians to degrees
            thetaPos = problem.thetaUsers(1:problem.nUsers)/(2*pi)*360;
        else
            % Deterministic Theta angles in degrees
            thetaPos = problem.thetaUsers(1:problem.nUsers);
        end
    end

    % Check if users were assigned Phi angles (deterministic)
    if ~conf.detLocation || ~isfield(problem,'phiUsers')
        % Generate Azymuth angles
        phiPos = random(t,1,problem.nUsers);
        if conf.verbosity >= 1
            fprintf('Users were not assigned Theta values (positions) in the space.\n');
        end
	elseif conf.useCasesLocation
        % Generate Azymuth angles from Use case
        phiPos = uc_az(conf.useCaseLoation,:);
    else
        % Retrieve Elevation angles from config file
        if problem.anglesInRadians
            % Convert Phi from radians to degrees
            phiPos = problem.phiUsers(1:problem.nUsers)/(2*pi)*360;
        else
            % Deterministic Phi angles in degrees
            phiPos = problem.phiUsers(1:problem.nUsers);
        end
    end

    % Check if users were assigned Distances (deterministic)
    if ~conf.detLocation || ~isfield(problem,'dUsers')
        dPos = rand(1,problem.nUsers) * (problem.maxdUsers-problem.mindUsers) + problem.mindUsers;
    elseif conf.useCasesLocation
        dPos = uc_dist(conf.useCasesLocation,:);
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