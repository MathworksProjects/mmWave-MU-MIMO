function conf = f_configuration
    % Define classes of traffic (Map to 5G applications)
    conf.Nclasses = 1;  % Number of classes (applications running)
    conf.Nusers = 2;
    % Define traffic types per users - For now, simplistic asumption of
    % same traffic per user
    for u = 1:conf.Nusers
        conf.class(u).iat = 100;  % Deterministic (constant) inter-arrival times (iat)
        conf.class(u).deadline = 100;  % Deadline in ms
        conf.class(u).numPkts = 100;  % number of packets for the class
        conf.class(u).Payload = 1500*8;  % Default payload in bits
    end
    conf.FLAGagg = true;  % Aggregate traffic
end