% We load the raw and unparsed results from the modified Matlab function
% "test_DMGOFDMPHYPacketErrorRateExample.m", which evaluates the PER for a
% given range of MCSs using a AWGN channel and the standard 802.11ad. This
% is taken as a baseline for the project.
% 
% The aim of this script is to parse that data and create the PERMCSTable 
% that will be used in the main. The final result is stored in data/.
%
% Usage example of the generated MCSPERTable:
%   snr = 1;   % Evaluating 1dB SNR
%   mcs = 13;  % Frame is configured using mcs 13
%   % The output of this script returns the statistical PER
%   mcsTable.MCSPERTable(snr==mcsTable.snrRange,mcs==mcsTable.mcsRange)

load('data/MCSPERTable_unparsed.mat');
snrRange = (-1:0.5:20);
mcsRange = (13:1:24);
MCSPERTable = zeros(length(snrRange),length(mcsRange));
for k=1:length(mcsRange)
    len = length(snrRanges{k});
    idVector = zeros(1,len);
    for i = 1:len
        idVector(i) = find(snrRange==snrRanges{k}(i));
    end
    MCSPERTable(idVector,k) = packetErrorRate(k,:).';
    if idVector(1)>1;                   MCSPERTable(1:idVector-1,k) = 1;  end
    if idVector(end)<length(snrRange);  MCSPERTable(idVector(end)+1:end,k) = 0;  end
end
mcsTable.MCSPERTable = MCSPERTable;
mcsTable.snrRange = snrRange;
mcsTable.mcsRange = mcsRange;
save('data/MCSPERTable','mcsTable');