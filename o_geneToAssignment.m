function [handle_Conf_Array,W,PRx,I] = o_geneToAssignment(gene,problem,conf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Extracting taper values from the input vector
    Taper_value = gene(problem.ant_elem+1:problem.ant_elem*2) .* ...
        exp(1i.*gene(problem.ant_elem*2+1:problem.ant_elem*3));

    % Creating a Conformal Array with cosine elements.
    % Conformal array will be limited to a single plane
    handle_Conf_Array = phased.ConformalArray('Element',problem.handle_Ant,...
        'ElementPosition',[zeros(1,problem.ant_elem);...
        problem.possible_locations(2,gene(1:problem.ant_elem));...
        problem.possible_locations(3,gene(1:problem.ant_elem))],...
        'Taper',Taper_value);

    %%
    PotRx = zeros(1,problem.nUsers);
    for u1=1:problem.nUsers
        if conf.multiPath
            for i=1:problem.maxnChannelPaths
                if problem.phiChannels(u1,i) ~= -Inf
                    PotRx(u1) = PotRx(u1) + problem.alphaChannels(u1,i)*...
                        pattern(handle_Conf_Array,problem.freq,...
                        problem.phiChannels(u1,i),...
                        problem.thetaChannels(u1,i),...
                        'Type','Power','Normalize',false)*...
                        (problem.lambda/(4*pi*problem.dUsers(u1)))^2;
                        % We assume that alphaChannels include the increase in 
                        % distance from dUsers(u1) to the real distance traver-
                        % sed by the signal in on each path
                end
            end
        else
            PotRx(u1) = pattern(handle_Conf_Array,problem.freq,...
                    problem.phiUsers(u1),...
                    problem.thetaUsers(u1),...
                    'Type','Powerdb','Normalize',false)*...
                    (problem.lambda/(4*pi*problem.dUsers(u1)))^2;
                    % We assume that alphaChannels include the increase in 
                    % distance from dUsers(u1) to the real distance traver-
                    % sed by the signal in on each path
        end
    end
    
    PotRx = 10*log10(PotRx);
    
    PRx = PotRx(problem.IDUserAssigned);
    PotRx(problem.IDUserAssigned) = [];
    % The interferences vector in the solution files does not contain the
    % interference inflicted to the user being analyzed [U] (it's nonsense)
    % Therefore, we need to shift the read index once we have read the
    % users with ID lower than U, and assign 0 to the interference
    % inflicted to himself.
    shift = 0;
    I = zeros(1,problem.nUsers);
    for m=1:problem.nUsers
        if m == problem.IDUserAssigned
            I(m) = 0;
            shift = -1;
        else
            I(m) = PotRx(m+shift);
        end
    end
    n_selected = length(gene)/3;
    amplitude = gene(n_selected+1:2*n_selected);
    phase = gene(2*n_selected+1:end);
    W = zeros(1,problem.NxPatch*problem.NyPatch);
    W(gene(1:n_selected)) = amplitude.*cos(phase) + 1i*amplitude.*sin(phase);
end

