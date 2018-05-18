% This function selects the MCS upon the estimated SNR with the optimal
% Beamformin configuration. The PER-MCS mapping is determined previously
% and stored in problem.MCSPER. A primary approach is taken using the AWGN
% channel, but we are still developping the same tables using the 5GNGR
% channel model available in Matlab.
% INPUTS
% - candSet: User ID being considered in the current slot
% - SNRList: SNR estimated by the heuristics (dB)
% - PERtarget: The maximum tolerable PER.
% - DEBUG: If True, it plots the results
% OUTPUTS:
% - MCS: Selected MCS
% - PER: Statistical PER estimated for the given configuration.
function [MCS,PER] = f_selectMCS(candSet,SNRList,PERtarget,MCSPER,DEBUG)
    snrRange = MCSPER.snrRange;
    mcsRange = MCSPER.mcsRange;
    table = MCSPER.MCSPERTable;
    strategy = 'aggressive';  % Possibilities are 'aggressive', 'conservative' and 'middle'
    risky = 'false';
    
    Nusers = length(SNRList);
    MCS = zeros(1,Nusers);
    PER = zeros(1,Nusers);
    for k =1:Nusers
        [~,snrIdx] = min(abs(SNRList(k)-snrRange));
        possibilities = table(snrIdx,:);
        tempPoss = possibilities<=PERtarget;
        if strcmp(strategy,'aggressive') && any(tempPoss~=0)
            % Aggressive - Risk it with the highest
            mcsIdx = find(possibilities<=PERtarget,1,'last');
        elseif strcmp(strategy,'conservative') && any(tempPoss~=0)
            % Conservative - Take always the lowest
            mcsIdx = find(possibilities<=PERtarget,1,'first');
        elseif strcmp(strategy,'middle') && any(tempPoss~=0)
            % Middle - amongst the possibilities, we select the one
            % standing right in the middle
            mcsList = find(possibilities<=PERtarget,length(mcsRange));
            mcsIdx = ceil(mean(mcsList));
        else
            mcsIdx = 1;
            risky = true;
        end
        MCS(k) = mcsRange(mcsIdx);
        PER(k) = table(snrIdx,mcsIdx);
        
        if DEBUG==1 * risky==0
            fprintf('\t\t\tID=%d\tMCS=%d\tPER=%.1f(%%)\tSNR=%.2f\n',candSet(k),MCS(k),PER(k).*1e2,SNRList(k));
        elseif DEBUG==1 && risky==1
            fprintf('\t\t\tID=%d\tMCS=%d\tPER=%.1f(%%)\tSNR=%.2f  (Failed to meet PER required)\n',candSet(k),MCS(k),PER(k).*1e2,SNRList(k));
        end
    end
end