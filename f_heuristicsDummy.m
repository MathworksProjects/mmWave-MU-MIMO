%% [tentativeIds,tentativeTH] = f_heuristicsDummy(combIDs,combTH)
% Dummy function that emulates the real heuristics (developped by Santi
% Rodrigo)
% function [tentativeIds,tentativeTH] = f_heuristicsDummy(combIDs,combTH)
function [estObj] = f_heuristicsDummy(MinObjF,MinObjFIsSNR,snrRange)
    if MinObjFIsSNR    
        MinObjFdB = max(10*log10(MinObjF),min(snrRange));  % in dB
        MaxObjFdB = 30;  % in dB
        estObjdB  = rand(1,length(MinObjFdB))*range([MinObjFdB MaxObjFdB]) + MinObjFdB;  % SNR (dB)
        estObj    = 10.^(estObjdB/10);  % SNR (linear)
        estObj    = MinObjF;  % SNR (linear)
    else
        estObj = MinObjF;  % Capacity in bps/s/Hz
    end
end