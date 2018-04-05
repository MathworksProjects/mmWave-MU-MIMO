function problem = f_configuration(problem)
    % Define classes of traffic (Map to 5G applications)
    problem.Nclasses = 1;  % Number of classes (applications running)
	problem.FLAGagg = true;  % Aggregate traffic
    problem.DEBUG = true;    % Enable Debugger and logs
    % Define traffic types per users - For now, simplistic asumption of
    % same traffic per user
    for u = 1:problem.nUsers
        problem.class(u).iat = 100;  % Deterministic (constant) inter-arrival times (iat)
        problem.class(u).deadline = 50;  % Deadline in ms
        problem.class(u).numPkts = 100;  % number of packets for the class
        problem.class(u).Payload = 1500*8;  % Default payload in bits
    end
end