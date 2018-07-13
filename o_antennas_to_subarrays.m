function relevant_partition = o_antennas_to_subarrays(antenna_list,partition)
% O_ANTENNAS_TO_SUBARRAYS - Converts antenna elements into a list of
% partition indices that contain them
%
% Syntax:  subarrays = o_antennas_to_subarrays(antenna_list,partition)
%
% Inputs:
%    antenna_list - List containing the antenna elements
%    partition - Cell containing the elements in each subarray partition
%
% Outputs:
%    relevant_partition - Relevant indices in the partition list where
%                         elements in antenna_list belong
%
% Example: 
%    antenna_list = [1 2 3 4];
%    partition{1} = [1 2];
%    partition{2} = [3 4];
%    partition{3} = [5 6];
%    relevant_partition = o_antennas_to_subarrays(antenna_list,partition)
%    display(relevant_partition);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: f_heuristics.m

%------------- BEGIN CODE --------------

    relevant_partition = [];
    for partID=1:length(partition)
        if all(ismember(partition{partID}, antenna_list))
            relevant_partition = [relevant_partition partID];  %#ok<AGROW>
        end
    end
end