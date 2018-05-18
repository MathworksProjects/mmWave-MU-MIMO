function finalSet = f_PER(candSet, problem, W, TXbits, MCS, channel_handles)
%This function conduct communication layer experiments 
%   This function takes in several parameters,
%   1. candSet -- potential candidate set; numerical vector, index start with 1;
%   2. problem -- struct
%   3. W -- the weights matrix, shall be nUsers x (nRow*nCol elements) dim
%   4. TXbits -- length of payload; integer number
%   5. MCS -- modulation and coding scheme, ranging from 1 to 12
%   6. channel_handles -- CDLChannel handles stored in a cell, indexed by
%   specific user
%   Generates one logical array, contains true or false, indicating for
%   that user, the packet go through or not

%% All parameters needed from upper layer
% PHY
lengthPSDU = TXbits;

% PHASED
nTx_row = problem.NxPatch;
nTx_col = problem.NyPatch;

% Simulation scenarios
nUsers = length(candSet);
angleToRx = phitheta2azel([problem.phiUsers; problem.thetaUsers]);
distance = problem.dUsers(candSet); %% !!!ASSUME CANDSET START WITH 1 WHEN COUNTING!!!

%% Parameters privite to this function
SNR = 10;
nRx = 1;

%% Configure system objects
tx_phy = s_phy_tx( ...
    'MCS', MCS, ...
    'PSDULength', lengthPSDU);
tx_pha = s_phased_tx( ...
    'numTxElements_row', nTx_row, ...
    'numTxElements_col', nTx_col, ...
    'txGain',            20, ...
    'visualization', false);
channel = s_phased_channel_handle_version( ...
    'SNR',               SNR, ...
    'applyPathLoss',     true);
rx_pha = s_phased_rx( ...
    'numRxElements', nRx, ...
    'rxGain',        20);
rx_phy = s_phy_rx();

%% Simulations
psdu = cell(nUsers, 1);
txWaveforms = cell(nUsers, 1);
response = cell(nUsers, 1);
finalSet = false(nUsers, 1);

for user_iter = 1 : nUsers
    psdu{user_iter} = randi([0 1], lengthPSDU * 8, 1);
    [txSymbols, cfgDMG] = tx_phy(psdu{user_iter});
    [txWaveforms{user_iter}, response{user_iter}] = tx_pha(txSymbols, angleToRx(:, user_iter), W(:, user_iter));
end

for outer_iter = 1 : nUsers
    combined_tx_waveforms = txWaveforms{outer_iter};
    for inner_iter = 1 : nUsers
        if inner_iter ~= outer_iter
            combined_tx_waveforms = combined_tx_waveforms + txWaveforms{inner_iter} * response{outer_iter};
        end
    end

    txWaveforms_afterChannel = channel(combined_tx_waveforms, distance(outer_iter), channel_handles{outer_iter});
    rxSymbols = rx_pha(txWaveforms_afterChannel, angleToRx(outer_iter));
    [psdu_rx, rxflag] = rx_phy(rxSymbols, cfgDMG);

    if ~isempty(psdu_rx)
        finalSet(outer_iter) = ~(any(biterr(psdu{outer_iter}, psdu_rx)) && rxflag);
    end
end
end