function [MCS,PER, RATE] = f_selectMCS(candSet,SNRList_lin,PERtarget,MCSPER,mcsPolicy,DEBUG)
% f_selectMCS - This function selects the MCS upon the estimated SNR with
% the optimal Beamformin configuration. The PER-MCS mapping is determined
% previously and stored in problem.MCSPER. A primary approach is taken
% using the AWGN channel, but we are still developping the same tables
% using the 5GNGR channel model available in Matlab.
%
% Syntax:  [MCS,PER] = f_selectMCS(candSet,SNRList,PERtarget,MCSPER,DEBUG)
%
% Inputs:
%    candSet -  User ID being considered in the current slot.
%    SNRList_lin - SNR estimated by the heuristics (dB).
%    PERtarget - The maximum tolerable PER.
%    mcsPolicy - Policy as to how the MCS is selected.
%    DEBUG -  If True, it plots the results.
%
% Outputs:
%    MCS - Selected MCS.
%    PER - Statistical PER estimated for the given configuration.
%    RATE - Data rate specified by 802.11ad in bps.
%
% Example: 
%    candSet = [1 2 3 4 5];
%    threshold = [0.1 0.2 0.3 0.7 0.4];
%    finalSet = f_PERtentative(candSet,threshold)
%    disp('finalSet')
%
% Other m-files required: f_configuration
% Subfunctions: none
% MAT-files required: none
%
% See also: f_configuration, main
%
%------------- BEGIN CODE --------------

snrRange = MCSPER.snrRange;
mcsRange = MCSPER.mcsRange;
table = MCSPER.MCSPERTable;
risky = 'false';

SNRList = pow2db(SNRList_lin);  % Convert linear to dB
Nusers = length(SNRList);  % Number of selected users to be served
MCS = zeros(1,Nusers);  % Final MCS selection
PER = zeros(1,Nusers);  % Statistical PER after selection
RATE = zeros(1,Nusers);  % Data rate specified by 802.11ad in bps
for k =1:Nusers
    [~,snrIdx] = min(abs(SNRList(k)-snrRange));
    possibilities = table(snrIdx,:);
    tempPoss = possibilities<=PERtarget;
    if strcmp(mcsPolicy,'aggressive') && any(tempPoss~=0)
        % Aggressive - Risk it with the highest
        mcsIdx = find(possibilities<=PERtarget,1,'last');
    elseif strcmp(mcsPolicy,'conservative') && any(tempPoss~=0)
        % Conservative - Take always the lowest
        mcsIdx = find(possibilities<=PERtarget,1,'first');
    elseif strcmp(mcsPolicy,'middle') && any(tempPoss~=0)
        % Middle - amongst the possibilities, we select the one
        % standing right in the middle
        mcsList = find(possibilities<=PERtarget,length(mcsRange));
        mcsIdx = ceil(mean(mcsList));
    elseif strcmp(mcsPolicy,'random') && any(tempPoss~=0)
        mcsIdx = ceil(rand*length(mcsRange));
    else
        mcsIdx = 1;
        risky = true;
    end
    % Output
    MCS(k)  = mcsRange(mcsIdx);
    PER(k)  = table(snrIdx,mcsIdx);
    RATE(k) = achievableRate(mcsIdx);

    if DEBUG==1 * risky==0
        fprintf('\t\t\tID=%d\tMCS=%d\tPER=%.1f(%%)\tSNR=%.2f(dB)\n',candSet(k),MCS(k),PER(k).*1e2,pow2db(SNRList_lin(k)));
    elseif DEBUG==1 && risky==1
        fprintf('\t\t\tID=%d\tMCS=%d\tPER=%.1f(%%)\tSNR=%.2f(dB)  (Failed to meet PER required)\n',candSet(k),MCS(k),PER(k).*1e2,pow2db(SNRList_lin(k)));
    end
end
end

function rate = achievableRate(MCSIdx)
    % Achievable Data Rates in bps
    % http://www.rfwireless-world.com/Tutorials/WLAN-802-11ad-tutorial.html
    table = [385e6, ...
             770e6, ...
             962.5e6, ...
             1155e6, ...
             1251.25e6, ...
             1540e6, ...
             1925e6, ...
             2310e6, ...
             2502e6, ...
             3080e6, ...
             3850e6, ...
             4620e6];
    rate = table(MCSIdx);
end


% EOF