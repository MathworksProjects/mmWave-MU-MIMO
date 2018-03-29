function finalSet = f_PER(candSet,problem,TXbits,W)
%f_PER tries for every candidate set, will the packet go through given
%inputs?
%   Supposely, 
%   *candSet* is a set with candidates (Carlos: please tell me what's inside); 
%   *problem* defines simulation scenarios -- the info I may need are,
%       1. nTx_row -- how many elements in URA row-wise
%       2. nTx_col -- how many elements in URA col-wise
%   *TxBits* defines the packet length in PSDU;
%   *W* the weights, a complex-valued double matrix with [nTx_row *
%   nTx_col, 1] dimension
%   Special notes -- this function support codegen and par-for, i.e., can
%   be run parallely to evaluate all candidates in candSet.

%#codegen

%% Inputs conversions
lengthPSDU = TXbits;

%% Unknown input-structures conversions
nTx_row = 8;
nTx_col = 8;

%% Validations -- very important
% validateattributes(W, {'double'}, {'2d','finite', 'complex', 'ncols', nTx_row, 'nrows', nTx_col});

%% Simulation constants
SNR = 5;
MCS = 1;
totPkt = 100;
numErr = 0;

%% Constructors
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

goThrough = false;
for i = 1 : totPkt %% Par-for can be costly -- overhead is huge
    psdu = randi([0 1], lengthPSDU*8, 1);
    [pkt, cfgDMG] = tx_phy(psdu);
    waveform = tx_pha(pkt, randi([-30 30]));
    waveforms = channel_pha(waveform);
    [psdu_rx, rxflag] = rx_phy(waveforms, [], cfgDMG);
    goThrough = any(biterr(psdu,psdu_rx)) && rxFlag;
end
% fprintf('MCS = %d, SNR = %.2fdB, PER = %.2f%%\n', MCS, SNR, numErr * 100 / totPkt);

end

