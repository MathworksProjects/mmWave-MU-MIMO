function [new_partition] = o_delete_subarrays_from_partition(partition,subarrays)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    new_partition = partition;
    positions_to_delete = zeros(1,length(subarrays));
    positions_iterator = 1;
    for n=1:length(partition)
        if(sum(ismember(subarrays,partition{n})) > 0)
            positions_to_delete(positions_iterator) = n;
            positions_iterator = positions_iterator + 1;
        end
    end
    new_partition(:,positions_to_delete) = [];
end

