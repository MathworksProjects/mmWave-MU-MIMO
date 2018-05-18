function [InitialValue] = o_generateInitialValue(NAntToBeAssigned,problem,conf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if strcmp(conf.genStructure, 'nchoosek') %CHECK
        InitialValue = [randi([1 nchoosek(problem.N_Subarrays,...
            NAntToBeAssigned)],1,1),...
            ones(1,NAntToBeAssigned),...
            zeros(1,NAntToBeAssigned)];
    elseif strcmp(conf.genStructure, 'allAntennas') % CHECK
        positions = randperm(problem.N_Subarrays,NAntToBeAssigned);
        weights_amp = zeros(1,problem.N_Subarrays);
        weights_amp(positions) = 1;
        InitialValue = [weights_amp,...
            zeros(1,problem.N_Subarrays)];
    else
         InitialValue = [randperm(problem.N_Subarrays,NAntToBeAssigned),...
            ones(1,NAntToBeAssigned),...
            zeros(1,NAntToBeAssigned),rand];
    end
end

