classdef s_phased_channel < matlab.System
    % phased_channel Add summary here
    %
    % NOTE: When renaming the class name phased_channel, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.

    properties

    end

    properties(Nontunable)
        numInputElements_row = 8;
        numInputElements_col = 8;
        numOutputElements_col = 1;
        numOutputElements_row = 1;
        SNR = 5;
        center_frequency = 60e9;
    end

    properties(DiscreteState)

    end

    properties(Access = private)
        MIMOChannel;
        AWGNChannel;
        CDLChannel;
    end

    methods
        % Constructor
        function obj = s_phased_channel(varargin)
            setProperties(obj,nargin,varargin{:})
        end
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            
%             gc.Doppler = 5;
%             gc.chanSRate = 1760e6 / 10; 
%             DelaySpread = 16;
%             numPaths = 5;
%             gc.PathDelays = floor(linspace(0,DelaySpread,numPaths))*(1/gc.chanSRate);
%             gc.PathGains  = zeros(size(gc.PathDelays));
%             for n=2:numPaths
%                 gc.PathGains(n) = gc.PathGains(n-1)-abs(randn);
%             end
%             
%             obj.MIMOChannel = comm.MIMOChannel( ...
%                 'SampleRate',                   gc.chanSRate, ... 
%                 'MaximumDopplerShift',          gc.Doppler, ...
%                 'PathDelays',                   floor(linspace(0, DelaySpread, numPaths))*(1/gc.chanSRate), ...
%                 'AveragePathGains',             gc.PathGains,...
%                 'NumTransmitAntennas',          obj.numTxElements, ...
%                 'NumReceiveAntennas',           obj.numRxElements, ...
%                 'PathGainsOutputPort',          true,...
%                 'NormalizePathGains',           true,...
%                 'NormalizeChannelOutputs',      true, ...
%                 'SpatialCorrelationSpecification', 'None');
            obj.AWGNChannel = comm.AWGNChannel( ...
                'NoiseMethod',                  'Signal to noise ratio (SNR)', ...
                'SNR',                          obj.SNR);
            
            cdlChan = nr5gCDLChannel;
            cdlChan.SampleRate = 1760e6 / 10; %% SC-1760e6 sa/s, OFDM 2640e6 sa/s
            cdlChan.TransmitAntennaArray.Size = [obj.numInputElements_row, obj.numInputElements_col, 1, 1, 1];
            cdlChan.ReceiveAntennaArray.Size = [obj.numOutputElements_col, obj.numOutputElements_row, 1, 1, 1];
            cdlChan.DelayProfile = 'CDL-C';
            cdlChan.DelaySpread = 100e-9;
            cdlChan.CarrierFrequency = obj.center_frequency;
            cdlChan.MaximumDopplerShift = (0 / 3.6) / physconst('lightspeed') * obj.center_frequency;
            
            obj.CDLChannel = cdlChan;
            
        end

        function rxWaveforms = stepImpl(obj, txWaveforms)
%             rxWaveforms = obj.MIMOChannel(txWaveforms);
            rxWaveforms = obj.CDLChannel(txWaveforms);
            rxWaveforms = obj.AWGNChannel(rxWaveforms);
        end

        function resetImpl(obj)

        end
    end
end
