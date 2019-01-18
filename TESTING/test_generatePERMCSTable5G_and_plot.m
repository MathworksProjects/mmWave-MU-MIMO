% Usage example of the generated MCSPERTable:
%   snr = 6;   % Evaluating 1dB SNR
%   mcs = 8;  % Frame is configured using mcs 13
%   % The output of this script returns the statistical PER, in percent
%   mcsTable.MCSPERTable(snr==mcsTable.snrRange,mcs==mcsTable.mcsRange)

addpath('../UTILITIES','-end');  % Add utilities folder at the end of search path
addpath('../code-systems','-end');  % Add system's folder at the end of search path
addpath('../code-beamforming','-end');  % Add beamforming folder at the end of search path
addpath('../code-wirelessEmulation','-end');  % Add channel folder at the end of search path

T = readtable('PER_TABLES/CDL_ALL_DelaySpread16ns.txt');
profile = 'CDL-C';

mcsTable.snrRange = unique(T.Var3);
mcsTable.mcsRange = unique(T.Var2);
index = find(contains(T.Var1, profile));
tempTable = sortrows(T{index, 2:end}, 1);
tempTable = sortrows(tempTable, 2);

mcsTable.MCSPERTable = reshape(tempTable(:, 3), length(mcsTable.mcsRange), length(mcsTable.snrRange))';

markers = 'ox*sd^v><ph+';
color = 'bmcrgbrkymcr';
figure;

for mcs = 1 : length(mcsTable.mcsRange)
    semilogy(mcsTable.snrRange, mcsTable.MCSPERTable(:, mcsTable.mcsRange(mcs))/100, ['-' markers(mcs) color(mcs)], 'LineWidth', 2);
    hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('PER');
dataStr = arrayfun(@(x)sprintf('MCS %d', x), mcsTable.mcsRange,'UniformOutput',false);
legend(dataStr, 'Location', 'best');
title(['PER for DMG SC-PHY with 3GPP TR 38.901 Channel, ', profile]);