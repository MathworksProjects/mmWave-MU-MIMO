function [Cap_tot] = o_compute_avCap_maxminthr_LAC(selection,PRx,I,Noise,MaxThr,MinThr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    interf = zeros(1,length(selection));
    PRx_sel = zeros(1,length(selection));
    for u2=1:length(selection)
        for u1=1:length(selection)
            if selection(u1) ~= 0 && I(u1,selection(u1),u2) ~= 0
                interf(u2) = interf(u2) + 10^(I(u1,selection(u1),u2)/10); % Interference inflicted in u2 by the selection of u1
            end
        end
        if selection(u2) ~= 0
            PRx_sel(u2) = 10^(PRx(u2,selection(u2))/10);
        end
    end
    SNR = PRx_sel./(interf+10^(Noise/10)); % In watts for the capacity computation!!
    Cap = log2(1+SNR);
    % If there is any Cap below the MinThr, or above the MaxThr (we
    % consider MaxThr and MinThr to be in bits/Hz...), the total Cap is
    % capped to -Inf (infeasible)
    if any((Cap-MinThr)<0) || any((MaxThr-Cap)<0)
        Cap_tot = -Inf;
    else
        Cap_tot = sum(Cap)/length(selection);
    end
end

