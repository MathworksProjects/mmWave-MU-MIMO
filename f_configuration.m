function [problem,traffic,flows] = f_configuration(conf,problem)
    % Geographic distribution of users
    [problem.thetaUsers, problem.phiUsers, problem.dUsers] = ...
        o_generate_positions(conf, problem.nUsers, problem.maxdUsers,problem.mindUsers);
    % Generate channels per user
    if conf.Use5GChannel
        % 5G 3GPP ETSI TR 38.901 compliant channel model
        [~, problem.thetaChannels, problem.phiChannels, problem.alphaChannels] = ...
                            o_generate_5Gchannels(conf,problem.nUsers,...
                            problem.thetaUsers, problem.phiUsers, ...
                            problem.NxPatch, problem.NyPatch, problem.freq,...
                            problem.maxnChannelPaths);
    else
        % Basic channel model (for the moment, including channel gains (alphas) !!!
        [problem.thetaChannels, problem.phiChannels, problem.alphaChannels] = ...
                            o_generate_channels(conf,problem.nUsers,...
                            problem.maxnChannelPaths);
    end
    % Configure Traffic
    problem = f_configureTraffic(problem);  % Struct with configuration parameters
    % Generate upper layer traffic
    [traffic,maxTime] = f_genDetTraffic(problem.class,problem.trafficType,problem.DEBUG);
    % Adequate Simulation time to the last packet arrival
    if maxTime~=0; problem.Tsym = maxTime; end
    % Generate PHY-flows: Convert traffic (arrivals) into individual Flow for
    % each user. Flows may overlap in time as the inter-arrival time may be
    % less than Tslot
    [flows] = f_arrivalToFlow(problem.Tslot,traffic);
end