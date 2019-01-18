% The purpose of this script is to simulate 802.11ad performance of MCS 1 - 12

clear; close all; clc;
addpath('../UTILITIES','-end');  % Add utilities folder at the end of search path
addpath('../code-systems','-end');  % Add system's folder at the end of search path
addpath('../code-beamforming','-end');  % Add beamforming folder at the end of search path
addpath('../code-wirelessEmulation','-end');  % Add channel folder at the end of search path
addpath('../data','-end');  % Add data folder at the end of search path

diary('SNR_TEMP.log');

lengthPSDU = 1000;

nTx_row = 1;
nTx_col = 1;
nRx = 1;

tic
for profile = {'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E'}
    for SNR = 0 : 20
        parfor MCS = 1 : 12
            tx_phy = s_phy_tx( ...
                'MCS', MCS, ...
                'PSDULength', lengthPSDU);
            channel_pha = s_phased_channel( ...
                'numInputElements_row',     1, ...
                'numInputElements_col',     1, ...
                'numOutputElements_row',    1, ...
                'numOutputElements_col',    1, ...
                'SNR',                      SNR, ...
                'applyPathLoss',            false, ...
                'profile',                  cell2mat(profile));
            rx_phy = s_phy_rx();
            
            totPkt = 2000;
            numError = 0;
            
            for i = 1 : totPkt
                psdu = randi([0 1], lengthPSDU * 8, 1);
                [txSymbols, cfgDMG] = tx_phy(psdu);
                txWaveforms_afterChannel = channel_pha(txSymbols, 0);
                [psdu_rx, rxflag] = rx_phy(txWaveforms_afterChannel, cfgDMG);
                
                if ~isempty(psdu_rx)
                    bitErrorFlag = any(biterr(psdu, psdu_rx));
                else
                    bitErrorFlag = true;
                end
                % if we receive nothing, or rxFlag is false, 
                % or we have a bit error, we claim packet error
                if isempty(psdu_rx) || ~rxflag || bitErrorFlag 
                    numError = numError + 1;
                end
            end
            fprintf('%s\t%d\t%d\t%.2f\n', cell2mat(profile), MCS, SNR, numError * 100 / totPkt);
        end
    end
end
toc
diary off