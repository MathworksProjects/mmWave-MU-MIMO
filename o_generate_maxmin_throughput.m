function [MaxThr,MinThr] = o_generate_maxmin_throughput(nUsers,...
            maxMinThr,minMinThr,cteMaxThr)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    MaxThr = repmat(cteMaxThr,1,nUsers);
    pd = makedist('Normal');
    pd.sigma = (maxMinThr-minMinThr)/4;
    pd.mu = (maxMinThr-minMinThr)/2 + minMinThr;
    t = truncate(pd,minMinThr,maxMinThr);
    MinThr = random(t,1,nUsers);
end

