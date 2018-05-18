function [thetaChannels, phiChannels, alphaChannels] = o_generate_channels(conf,nUsers,...
                    maxnChannelPaths)
% Boolean to check whether we have already selected the channels to remove
    % for each user
    pd = makedist('Normal');
    pd.sigma = 45;
    t = truncate(pd,-45,45);
    if conf.verbosity >= 1
        fprintf('===================================\n');
        fprintf('The channels were not assigned AoA:\n');
        fprintf('Assigning random values...\n');
    end
    thetaChannels = random(t,nUsers,maxnChannelPaths);
    if conf.verbosity >= 1
        fprintf('New elevations assigned:\n');
        display(thetaChannels);
    end
    
    phiChannels = random(t,nUsers,maxnChannelPaths);
    if conf.verbosity >= 1
        fprintf('New azimuths assigned:\n');
        display(phiChannels);
    end
    % Finally, alpha gains of each path
    pd = makedist('Normal');
    pd.mu = 0.5;
    pd.sigma = 0.25;
    t = truncate(pd,0,1);
    alphaChannels = random(t,nUsers,maxnChannelPaths);
    userWithMaxNChannels = randi(nUsers);
    indexToRemove = ones(1,nUsers)*(maxnChannelPaths+1);
    for i=setdiff(1:nUsers,userWithMaxNChannels)
        % Now we select some channels (from indexToRemove to end) and remove 
        % them (not every user should have maxnChannelPaths
        indexToRemove(i) = max(1,randi(maxnChannelPaths)+1);
        thetaChannels(i,indexToRemove(i):end) = -Inf;
        phiChannels(i,indexToRemove(i):end) = -Inf;
        alphaChannels(i,indexToRemove(i):end) = -Inf;
    end
    if conf.verbosity >= 1
        fprintf('New gains assigned:\n');
        display(alphaChannels);
    end
end

