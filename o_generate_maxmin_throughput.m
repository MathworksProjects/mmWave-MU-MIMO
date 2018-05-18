function [MaxObjF,MinObjF] = o_generate_maxmin_throughput(nUsers,...
            maxMinObjF,minMinObjF,cteMaxObjF)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    MaxObjF = repmat(cteMaxObjF,1,nUsers);
    pd = makedist('Normal');
    pd.sigma = (maxMinObjF-minMinObjF)/4;
    pd.mu = (maxMinObjF-minMinObjF)/2 + minMinObjF;
    t = truncate(pd,minMinObjF,maxMinObjF);
    MinObjF = random(t,1,nUsers);
end

