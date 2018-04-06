function subarrays = antennas_to_subarrays(antenna_list,partition)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
subarrays = [];
for n=1:length(partition)
    if(sum(ismember(antenna_list,partition{n})) > 0)
        subarrays = [subarrays, n];
    end
end

