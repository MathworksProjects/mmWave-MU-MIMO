function finalSet = f_PER(candSet, problem, W, TXbits, MCS, channel_handles, array_handle)
% f_PER - This function conduct communication layer experiments
%
% Syntax:  finalSet = f_PER(candSet, problem, W, TXbits, MCS, ...
%                           channel_handles, array_handle)
%
% Inputs:
%	candSet - potential candidate set; numerical vector, index start with 1;
%   problem - struct
%   W - the weights matrix, shall be nUsers x (nRow*nCol elements) dim
%   TXbits - length of payload; integer number
%   MCS - modulation and coding scheme, ranging from 1 to 12
%   channel_handles - CDLChannel handles stored in a cell, indexed by users
%   array_handle - handle of phased array definition (conformal array)
%
% Outputs:
%   finalSet - Generates one logical array, contains true or false, indicating 
%              for that user, the packet go through or not
%
% Other m-files required: s_* functions
% Subfunctions: none
% MAT-files required: none
%
%------------- BEGIN CODE --------------

%% All parameters needed from upper layer
% PHY
lengthPSDU = ceil(TXbits / 8); % in bytes
noiseFigure = problem.Noise;

% PHASED
nTx_row = problem.NxPatch;
nTx_col = problem.NyPatch;
centerfreq = problem.freq;

% Simulation scenarios
nUsers = length(candSet);
% angleToRx = phitheta2azel([problem.phiUsers; problem.thetaUsers]);
angleToRx = [problem.phiUsers; problem.thetaUsers];
distance_3d = problem.dUsers(candSet); %% !!!ASSUME CANDSET START WITH 1 WHEN COUNTING!!!
distance_2d = problem.dUsers(candSet); %% not so precise

%% Parameters privite to this function
nRx = 1;

%% Configure system objects
for i = 1 : nUsers
    tx_phy{i} = s_phy_tx( 'PSDULength', lengthPSDU(i), 'MCS', MCS(i));
end
tx_pha = s_phased_tx( ...
    'numTxElements_row',    nTx_row, ...
    'numTxElements_col',    nTx_col, ...
    'antenna_array',        array_handle, ...
    'txGain',               0, ...
    'center_frequency',     centerfreq, ...
    'visualization',        false);
channel = s_phased_channel_handle_version( ...
    'noiseFigure',          noiseFigure, ...
    'center_frequency',     centerfreq, ...
    'applyPathLoss',        false);
rx_pha = s_phased_rx( ...
    'numRxElements',        nRx, ...
    'rxGain',               0);
rx_phy = s_phy_rx();

% resp = cell(nUsers, 1);
possible_locations = array_handle.getElementPosition;
for id = 1:problem.nUsers
    relevant_positions = (W(id,:)~=0);
    Taper_user = W(id,relevant_positions);
    
    handle_Conf_Array_USER{id} = phased.ConformalArray(...
        'Element',array_handle.Element,...
        'ElementPosition', [possible_locations(1,relevant_positions);...
        possible_locations(2,relevant_positions);...
        possible_locations(3,relevant_positions)],...
        'Taper',Taper_user);
    %     resp{id} = phased.ArrayResponse( ...
    %         'SensorArray', handle_Conf_Array_USER{id}, ...
    %         'WeightsInputPort', true);
end

%% Simulations
psdu = cell(nUsers, 1);
txWaveforms = cell(nUsers, 1);
finalSet = [];
waveform_size = zeros(nUsers, 1);

for user_iter = 1 : nUsers
    psdu{user_iter} = randi([0 1], lengthPSDU(user_iter) * 8, 1, 'int8');
    [txSymbols, cfgDMG] = tx_phy{user_iter}(psdu{user_iter});
    txWaveforms{user_iter} = tx_pha(txSymbols, angleToRx(:, user_iter), W(user_iter, :).');
    waveform_size(user_iter) = size(txSymbols, 1);
end
maximumSize = max(waveform_size);

% Get the response -- response is a nUsers x nUsers, e.g., response(1, 3)
% means, response at angle @ user 1, with W specified by user 3, i.e., the
% response at interference angle @ user 1
response = zeros(nUsers, nUsers);
for outer_iter = 1 : nUsers
    for inner_iter = 1 : nUsers
        response(outer_iter, inner_iter) =  patternAzimuth(handle_Conf_Array_USER{outer_iter}, ...
            problem.freq,problem.thetaUsers(inner_iter),'Azimuth',problem.phiUsers(inner_iter),'Type','power');
        %          relevant_positions = (W(inner_iter,:)~=0);
        %          Taper_user = W(inner_iter,relevant_positions);
        %          resp_from_arrayResponse = abs(resp{inner_iter}(centerfreq, angleToRx(:, outer_iter), ones(8, 1))) .^ 2
    end
end

for outer_iter = 1 : nUsers
    combined_tx_waveforms = sqrt(response(outer_iter, outer_iter)) * [txWaveforms{outer_iter}; zeros(maximumSize - length(txWaveforms{outer_iter}), nTx_row * nTx_col)];
    for inner_iter = 1 : nUsers
        if inner_iter ~= outer_iter
            combined_tx_waveforms = combined_tx_waveforms + ...
                [txWaveforms{inner_iter}; zeros(maximumSize - length(txWaveforms{inner_iter}), nTx_row * nTx_col)] ...
                * sqrt(response(inner_iter, outer_iter));
        end
    end
    
    txWaveforms_afterChannel = channel(combined_tx_waveforms, channel_handles{outer_iter}, distance_3d(outer_iter), distance_2d(outer_iter));
    rxSymbols = rx_pha(txWaveforms_afterChannel, angleToRx(outer_iter));
    [psdu_rx, rxflag] = rx_phy(rxSymbols, cfgDMG);
    
    if ~isempty(psdu_rx)
        gothrough = ~(any(biterr(psdu{outer_iter}, psdu_rx)) && rxflag);
        if gothrough
            finalSet = [finalSet outer_iter];  %#ok<AGROW>
        end
    end
end
%------------- END CODE --------------
end