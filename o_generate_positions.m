function [thetaPos, phiPos, dPos] = o_generate_positions(conf,problem)
    % Generate Normal distribution
    pd = makedist('Normal');
    pd.sigma = 45;
    t = truncate(pd,-45,45);
    dUser = 5 * ones(1, problem.nUsers);
    
    % Generate Angles for use case (uc) distribution - dynamic
    Delta = 90/problem.nUsers;
    if mod(problem.nUsers,2)==0  % even
        ini = Delta/2;
        vect = [(-problem.nUsers/2:1:-1) (1:1:problem.nUsers/2)];
    else                         % odd
        ini = Delta;
        vect = [(-floor(problem.nUsers/2):1:-1) 0 (1:1:floor(problem.nUsers/2))];
    end
    vect = ini.*vect;
    % UC 1 - Located horizontally (no elevation)
    uc_el(1,:) = zeros(1,problem.nUsers);
    uc_az(1,:) = vect;
    uc_dist(1,:) = dUser.*ones(1,problem.nUsers);
    % UC 2 - Located horizontally (no elevation - a bit more separation)
    uc_el(2,:) = zeros(1,problem.nUsers);
    uc_az(2,:) = 1.5.*vect;
    uc_dist(2,:) = dUser.*ones(1,problem.nUsers);
    % UC 3 - Located vertically (no azymuth)
    uc_el(3,:) = vect;
    uc_az(3,:) = zeros(1,problem.nUsers);
    uc_dist(3,:) = dUser.*ones(1,problem.nUsers);
    % UC 4 - Located diagonaly
    uc_el(4,:) = vect;
    uc_az(4,:) = vect;
    uc_dist(4,:) = dUser.*ones(1,problem.nUsers);
    % UC 5 - Located vertically (15 deg azymuth)
    uc_el(5,:) = vect;
    uc_az(5,:) = 15.*ones(1,problem.nUsers);
    uc_dist(5,:) = dUser.*ones(1,problem.nUsers);
    % UC 6 - Located horizontally (15 deg elevation)
    uc_el(6,:) = 15.*ones(1,problem.nUsers);
    uc_az(6,:) = vect;
    uc_dist(6,:) = dUser.*ones(1,problem.nUsers);

    % Check if users were assigned Theta angles (deterministic)
    if ~conf.detLocation || ~isfield(problem,'thetaUsers')
        % Generate Random Elevation angles
        thetaPos = random(t,1,problem.nUsers);
        if conf.verbosity >= 1
            fprintf('Users were not assigned Theta values (positions) in the space.\n');
        end
    elseif conf.useCasesLocation
        % Generate Elevation angles from Use case
        thetaPos = uc_el(conf.useCaseLocation,:);
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
        phiPos = uc_az(conf.useCaseLocation,:);
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