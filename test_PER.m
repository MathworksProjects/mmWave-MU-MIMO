% clear; close all;

lengthPSDU = 1000;
nTx_row = 2;
nTx_col = 2;

tic
for SNR = 5
    for MCS = 1
        tx_phy = s_phy_tx( ...
            'MCS', MCS, ...
            'PSDULength', lengthPSDU);
        rx_phy = s_phy_rx();
        tx_pha = s_phased_tx( ...
            'numTxElements_row', nTx_row, ...
            'numTxElements_col', nTx_col, ...
            'visualization', false);
        channel_pha = s_phased_channel( ...
            'numInputElements_row',     nTx_row, ...
            'numInputElements_col',     nTx_col, ...
            'numOutputElements_row',    1, ...
            'numOutputElements_col',    1, ...
            'SNR',                      SNR);
        
        totPkt = 1;
        numErr = 0;
        
        for i = 1 : totPkt
            psdu = randi([0 1], lengthPSDU*8, 1);
            [pkt, cfgDMG] = tx_phy(psdu);
            waveform = tx_pha(pkt, randi([-30 30]), []);
            waveforms = channel_pha(waveform);
            [psdu_rx, rxflag] = rx_phy(waveforms, [], cfgDMG);
            if ~isempty(psdu_rx)
                numErr = any(biterr(psdu,psdu_rx)) + numErr;
            end
        end
        fprintf('MCS = %d, SNR = %.2fdB, PER = %.2f%%\n', MCS, SNR, numErr * 100 / totPkt);
        
    end
    
end
toc