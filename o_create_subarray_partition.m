function [problem] = o_create_subarray_partition(problem)
    % Now we fill in the partition cell, depending on the restriction
    % imposed on it:
    % - 'Localized': the subarrays are rectangular contiguous divisions of the
    %    array
    % - 'DiagInterleaved': the subarrays are interleaved in a diagonal fashion
    % - 'Interleaved': the subarrays are interleaved forming rectangular arrays
    %    with larger distance between elements.
    if strcmp(problem.arrayRestriction,'Localized')
        % partition contains the antenna arrangement into the different
        % subarrays, i.e. for each subarrays, you will find the set of antenna
        % indexes associated to it.
        problem.Partition = cell(1,problem.N_Subarrays);
        for nx=0:(problem.NxSubarrays-1)
            for ny=0:(problem.NySubarrays-1)
                x0 = floor(nx/problem.NxSubarrays*problem.NxPatch);
                xend = floor((nx+1)/problem.NxSubarrays*problem.NxPatch);
                for x=x0:(xend-1)
                    y0 = floor(ny/problem.NySubarrays*problem.NyPatch);
                    yend = floor((ny+1)/problem.NySubarrays*problem.NyPatch);
                    for y=y0:(yend-1)
                        problem.Partition{nx*problem.NySubarrays+ny+1} = ...
                            [problem.Partition{nx*problem.NySubarrays+ny+1},...
                            (x*problem.NyPatch+y+1)];
                    end
                end
            end
        end
    elseif strcmp(problem.arrayRestriction,'DiagInterleaved')
        % partition contains the antenna arrangement into the different
        % subarrays, i.e. for each subarrays, you will find the set of antenna
        % indexes associated to it.
        if problem.NxSubarrays ~= problem.NySubarrays
            if  (problem.NxSubarrays <= problem.NySubarrays && ...
                    problem.NxPatch ~= problem.NxSubarrays) || ...
                    (problem.NySubarrays <= problem.NxSubarrays && ...
                    problem.NyPatch ~= problem.NySubarrays)
                fprintf(['ERROR: problem.NxSubarrays (%d) and',...
                    ' problem.NySubarrays',...
                    ' (%d) do not match with a correct diagonally interleaved',...
                    ' configuration\n'], ...
                problem.NxSubarrays, problem.NySubarrays);
                return
            end
        end
        problem.N_Subarrays = max(problem.NxSubarrays,problem.NySubarrays);
        problem.Partition = cell(1,problem.N_Subarrays);
        for x=0:(problem.NxPatch-1)
            for y=0:(problem.NyPatch-1)
                problem.Partition{mod(x+y,problem.NySubarrays)+1} = ...
                    [problem.Partition{mod(x+y,problem.NySubarrays)+1},...
                    x*problem.NyPatch+y+1];
            end
        end
    elseif strcmp(problem.arrayRestriction,'Interleaved')
        % partition contains the antenna arrangement into the different
        % subarrays, i.e. for each subarrays, you will find the set of antenna
        % indexes associated to it.
        problem.Partition = cell(1,problem.N_Subarrays);
        for x=0:(problem.NxPatch-1)
            for y=0:(problem.NyPatch-1)
                problem.Partition{mod(y,problem.NySubarrays)+1+...
                    mod(x,problem.NxSubarrays)*problem.NySubarrays} = ...
                    [problem.Partition{mod(y,problem.NySubarrays)+1+...
                    mod(x,problem.NxSubarrays)*problem.NySubarrays},...
                    x*problem.NyPatch+y+1];
            end
        end    
    elseif strcmp(problem.arrayRestriction,'None')
        problem.Partition = cell(1,problem.N_Antennas);
        for a=1:problem.N_Antennas
            problem.Partition{a} = a;
        end
        problem.N_Subarrays = problem.N_Antennas;
    end
end

