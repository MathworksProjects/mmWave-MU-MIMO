classdef s_phy_rx < matlab.System
    % phy_rx Add summary here
    %
    % NOTE: When renaming the class name untitled9, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.

    properties

    end

    properties(Nontunable)
        PacketDetectThreshold = 0.5; % Good for low SNRs
    end

    properties(DiscreteState)

    end

    properties(Access = private)
        numPacketErrors
        numPkt
    end

    methods
        % Constructor
        function obj = s_phy_rx(varargin)
            setProperties(obj,nargin,varargin{:})
        end
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            obj.numPacketErrors = 0;
            obj.numPkt = 0;
        end

        function [rxPSDU, rxFlag] = stepImpl(obj, rxWaveform, channelGain, cfgDMG)
            ind = wlanFieldIndices(cfgDMG);
            fs = wlanSampleRate(cfgDMG);
            Ngi = 64; % Fixed GI length defined in the standard (20.6.3.2.5)
            
            % Packet detection
            pktStartOffset = dmgPacketDetect(rxWaveform,0,obj.PacketDetectThreshold);
            if isempty(pktStartOffset) % If empty no STF detected; packet error
                obj.numPacketErrors = obj.numPacketErrors+1;
                obj.numPkt = obj.numPkt+1;
                rxPSDU = [];
                rxFlag = false;
                return % Go to next loop iteration
            end

            % Frequency offset estimation and correction
            stf = rxWaveform(pktStartOffset+(ind.DMGSTF(1):ind.DMGSTF(2)));
            fOffsetEst = dmgCFOEstimate(stf);
            rxWaveform = helperFrequencyOffset(rxWaveform,fs,-fOffsetEst);

            % Symbol timing and channel estimate
            preamblefield = rxWaveform(pktStartOffset+1:pktStartOffset+ind.DMGHeader(2),:);
            [symbolTimingOffset,chanEst] = dmgTimingAndChannelEstimate(preamblefield);
            startOffset = pktStartOffset+symbolTimingOffset;

            % If not enough samples to decode detected data field start,
            % then assume synchronization error and packet error
            if (startOffset+ind.DMGData(2))>size(rxWaveform,1)
               obj.numPacketErrors = obj.numPacketErrors+1;
                obj.numPkt = obj.numPkt+1;
                rxPSDU = [];
                rxFlag = false;
                return
            end

            % Noise estimation using the STF as repeating sequence
            stf = rxWaveform(pktStartOffset+(ind.DMGSTF(1):ind.DMGSTF(2)));
            nVarEst = dmgSTFNoiseEstimate(stf);

            % Extract data field (ignore first GI)
            rxData = rxWaveform(startOffset+((ind.DMGData(1)+Ngi):ind.DMGData(2)));

            % Linear frequency domain equalization
            rxEqDataBlks = dmgSingleCarrierFDE(rxData,chanEst,nVarEst);

            % Unique word phase tracking
            rxEqDataBlks = dmgUniqueWordPhaseTracking(rxEqDataBlks);

            % Discard GI from all blocks
            rxDataSym = rxEqDataBlks(1:end-Ngi,:);

            % Recover the transmitted PSDU from DMG Data field
            rxPSDU = wlanDMGDataBitRecover(rxDataSym,nVarEst,cfgDMG);
            rxFlag = true;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end
