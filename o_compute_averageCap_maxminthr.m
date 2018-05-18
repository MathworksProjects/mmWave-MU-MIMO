function [aveCap, cap] = o_compute_averageCap_maxminthr(PRx,I,Noise,MaxObjF,MinObjF)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    nUsers = length(MinObjF);
    interf = zeros(1,nUsers); % linear!
    PRx_sel = zeros(1,nUsers); % linear!
    for u2=1:nUsers
        for u1=1:nUsers
            if I(u1,u2) ~= 0
                interf(u2) = interf(u2) + 10^(I(u1,u2)/10); % Interference inflicted in u2 by the selection of u1
            end
        end
        PRx_sel(u2) = 10^(PRx(u2)/10);
    end
    SNR = PRx_sel./(interf+10^(Noise/10)); % In watts for the capacity computation!!
    cap = log2(1+SNR);
    if any((cap-MinObjF)<0) || any((MaxObjF-cap)<0); aveCap = -Inf;
    else;                                            aveCap = sum(cap)/nUsers;
    end
end

