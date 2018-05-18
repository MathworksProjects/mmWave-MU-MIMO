%% 802.11ad Packet Error Rate Simulation for OFDM PHY
%
% This example shows how to measure the packet error rate of an IEEE(R)
% 802.11ad(TM) DMG OFDM PHY link using an end-to-end simulation with an
% AWGN channel.

% Copyright 2017 The MathWorks, Inc.

%% Introduction
% In this example an end-to-end simulation is used to determine the packet
% error rate for an 802.11ad DMG [ <#10 1> ] OFDM link with an AWGN channel
% at a selection of SNR points for a defined modulation and coding scheme
% (MCS). For each SNR point, multiple packets are transmitted through a
% channel, demodulated and the PSDUs recovered. The PSDUs are compared to
% those transmitted to determine the number of packet errors and hence the
% packet error rate. Perfect time and frequency synchronization is assumed
% in this example. The following diagram summarizes the processing for each
% packet.
%
% <<DMGOFDMPERDiagram.png>>
%
% This example also demonstrates how a <matlab:doc('parfor') parfor> loop
% can be used instead of the <matlab:doc('for') for> loop when simulating
% each SNR point to speed up a simulation. <matlab:doc('parfor') parfor>,
% as part of the Parallel Computing Toolbox(TM), executes processing for
% each SNR in parallel to reduce the total simulation time.

%% Waveform Configuration
% An 802.11ad DMG OFDM transmission is simulated in this example. The DMG
% format configuration object contains the format specific configuration of
% the transmission. The object is created using the
% <matlab:doc('wlanDMGConfig') wlanDMGConfig> function. The properties of
% the object contain the configuration. In this example the object is
% configured for an OFDM transmission with MCS 21 and an 8192 byte PSDU. If
% |mcs| is specified as a vector, the simulation is performed for each MCS
% element. The MCS determines the PHY type used, therefore the MCS must be
% within the range 13-24 to use the OFDM PHY.

% Create a format configuration object for a DMG OFDM transmission
cfgDMG = wlanDMGConfig;
cfgDMG.PSDULength = 8192; % bytes

% For DMG OFDM PHY, the valid range of MCS is 13-24(inclusive) 
mcs = 13:1:24; % OFDM PHY, 16QAM, rate 13/16
% mcs = 13;

%% Simulation Parameters
% For each SNR point a number of packets are generated, passed through a
% channel and demodulated to determine the packet error rate. The SNR
% points to simulate are selected from |snrRanges| based on the MCS to
% simulate. The SNR range for each MCS is selected in order to simulate the
% transition from all packets being decoded in error to all packets being
% decoded successfully as the SNR increases.

% SNR ranges to use for AWGN
snrRanges = {...
    -1:0.5:1.5, ...  % MCS 13
    0:0.5:2.5, ...   % MCS 14
    1.5:0.5:4, ...   % MCS 15
    3:0.5:5.5, ...   % MCS 16
    4.5:0.5:7, ...   % MCS 17
    7.5:0.5:10, ...  % MCS 18
    9:0.5:11.5, ...  % MCS 19
    10.5:0.5:13, ... % MCS 20
    12:0.5:14.5, ... % MCS 21
    14.5:0.5:17, ... % MCS 22
    16.5:0.5:19, ... % MCS 23
    17.5:0.5:20, ... % MCS 24
    };

%%
% The number of packets tested at each SNR point is controlled by two
% parameters:
%
% # |maxNumErrors| is the maximum number of packet errors simulated at each
% SNR point. When the number of packet errors reaches this limit, the
% simulation at this SNR point is complete.
% # |maxNumPackets| is the maximum number of packets simulated at each SNR
% point and limits the length of the simulation if the packet error limit
% is not reached. 
%
% The numbers chosen in this example will lead to a very short simulation.
% For meaningful results we recommend increasing these numbers.

maxNumErrors = 1e3;  % The maximum number of packet errors at an SNR point
maxNumPackets = 1e4; % Maximum number of packets at an SNR point

%%
% Set the remaining variables for the simulation.

% Indices of data and pilot occupied subcarriers
cfgDMG.MCS = mcs(1); % Set OFDM MCS to get subcarrier indices
scIndices = helperOccupiedSubcarrierIndices('DMG-Data',cfgDMG);
Nsd = numel(scIndices.Data); % Number of data carrying subcarriers

% OFDM information
ofdmInfo = helperOFDMInfo('DMG-Data',cfgDMG);

%% Processing SNR Points
% For each SNR point a number of packets are tested and the packet error
% rate calculated.
%
% For each packet the following processing steps occur:
%
% # A PSDU is created and encoded to create a single packet waveform.
% # AWGN is added to the waveform to create the desired average SNR per
% subcarrier after OFDM demodulation. <matlab:doc('comm.AWGNChannel')
% comm.AWGNChannel> is configured to provide the correct SNR. The
% configuration accounts for the noise energy in unused subcarriers which
% are removed during OFDM demodulation.
% # The DMG-Data field is extracted from the perfectly synchronized
% received waveform and OFDM demodulated.
% # The pilots are discarded and the remaining OFDM demodulated symbols are
% equalized using the known channel response. As an AWGN link is used in
% this example, the complex channel gain is assumed to be one for each
% subcarrier.
% # The PSDU is recovered from the equalized data symbols.
%
% A <matlab:doc('parfor') parfor> loop can be used to parallelize
% processing of the SNR points, therefore for each SNR point an AWGN
% channel is created and configured with <matlab:doc('comm.AWGNChannel')
% comm.AWGNChannel>. To enable the use of parallel computing for increased
% speed comment out the 'for' statement and uncomment the 'parfor'
% statement below.

numSNR = numel(snrRanges{1}); % Number of SNR points   
numMCS = numel(mcs);          % Number of MCS
packetErrorRate = zeros(numMCS,numSNR);

for imcs = 1:numMCS
    
    cfgDMG.MCS = mcs(imcs);
    if ~strcmp(phyType(cfgDMG),'OFDM')
        error('This example only supports DMG OFDM PHY simulation');
    end

    % Indices of fields within the packet
    fieldIndices = wlanFieldIndices(cfgDMG);
    
    % SNR points to simulate from MCS
    snr = snrRanges{cfgDMG.MCS-12};

    %parfor isnr = 1:numSNR % Use 'parfor' to speed up the simulation
    for isnr = 1:numSNR % Use 'for' to debug the simulation
        % Set random substream index per iteration to ensure that each
        % iteration uses a repeatable set of random numbers
        stream = RandStream('combRecursive','Seed',0);
        stream.Substream = isnr;
        RandStream.setGlobalStream(stream);

        % Create an instance of the AWGN channel per SNR point simulated
        awgnChannel = comm.AWGNChannel;
        awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
        awgnChannel.SignalPower = 1; 
        % Account for noise energy in nulls so the SNR is defined per
        % active subcarrier
        awgnChannel.SNR = snr(isnr)-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones); 

        % Loop to simulate multiple packets
        numPacketErrors = 0;
        numPkt = 1; % Index of packet transmitted
        while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
            % Generate a packet waveform
            txPSDU = randi([0 1],cfgDMG.PSDULength*8,1); % PSDULength in bytes
            tx = wlanWaveformGenerator(txPSDU,cfgDMG);

            % Pass the waveform through AWGN channel model
            rx = awgnChannel(tx);

            % Extract data field
            rxData = rx(fieldIndices.DMGData(1):fieldIndices.DMGData(2));

            % OFDM demodulate
            demodSym = helperOFDMDemodulate(rxData,'DMG-Data',cfgDMG);
            dataSym = demodSym(scIndices.Data,:); % Discard pilots

            % Equalize
            chanSym = ones(Nsd,1); % Set channel gains to 1 as AWGN channel
            nVar = 10^(-snr(isnr)/10); % Noise variance
            [eqSym,csi] = helperSymbolEqualize(dataSym,chanSym,nVar);

            % Recover data
            rxPSDU = wlanDMGDataBitRecover(eqSym,nVar,csi,cfgDMG);

            % Determine if any bits are in error, i.e. a packet error
            packetError = any(biterr(txPSDU,rxPSDU));
            numPacketErrors = numPacketErrors+packetError;
            numPkt = numPkt+1;
        end

        % Calculate packet error rate (PER) at SNR point
        packetErrorRate(imcs,isnr) = numPacketErrors/(numPkt-1);
        disp(['MCS ' num2str(mcs(imcs)) ','...
              ' SNR ' num2str(snr(isnr)) ...
              ' completed after ' num2str(numPkt-1) ' packets,'...
              ' PER:' num2str(packetErrorRate(imcs,isnr))]);
    end
end

save('TotExecMatlab80211ad','packetErrorRate','mcs','snrRanges');

%% Plot Packet Error Rate vs SNR Results

markers = 'ox*sd^v><ph+';
color = 'bmcrgbrkymcr';
figure;
for imcs = 1:numMCS
    semilogy(snrRanges{mcs(imcs)-12},packetErrorRate(imcs,:).',['-' markers(imcs) color(imcs)]);
    hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('PER');
dataStr = arrayfun(@(x)sprintf('MCS %d',x),mcs,'UniformOutput',false);
legend(dataStr);
title('PER for DMG OFDM-PHY with AWGN channel');

%% Further Exploration
% The number of packets tested at each SNR point is controlled by two
% parameters: |maxNumErrors| and |maxNumPackets|. For meaningful results
% these values should be larger than those presented in this example.
% Increasing the number of packets simulated allows the PER under different
% scenarios to be compared. As an example, the figure below was created by
% running the example for longer with |maxNumErrors = 1e3| and
% |maxNumPackets = 1e4|, for |mcs = 13:24|.
%
% <<DMGOFDMPERExample.png>>

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperOccupiedSubcarrierIndices.m') helperOccupiedSubcarrierIndices.m>
% * <matlab:edit('helperOFDMDemodulate.m') helperOFDMDemodulate.m>
% * <matlab:edit('helperOFDMInfo.m') helperOFDMInfo.m>
% * <matlab:edit('helperSymbolEqualize.m') helperSymbolEqualize.m>

%% Selected Bibliography
% # IEEE Std 802.11ad(TM)-2012 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications.
% Amendment 3: Enhancements for Very High Throughput in the 60 GHz Band.

displayEndOfDemoMessage(mfilename)