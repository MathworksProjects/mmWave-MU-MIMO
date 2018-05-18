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
        center_frequency = 60.48e9;
        isLoS = true;
    end
    
    properties(DiscreteState)
        
    end
    
    properties(Access = private)
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
            
            obj.AWGNChannel = comm.AWGNChannel( ...
                'NoiseMethod',                  'Signal to noise ratio (SNR)', ...
                'SNR',                          obj.SNR);
            
            cdlChan = nr5gCDLChannel;
            cdlChan.SampleRate = 1760e6 / 1; %% SC-1760e6 sa/s, OFDM 2640e6 sa/s
            cdlChan.TransmitAntennaArray.Size = [obj.numInputElements_row, obj.numInputElements_col, 1, 1, 1];
            cdlChan.ReceiveAntennaArray.Size = [obj.numOutputElements_col, obj.numOutputElements_row, 1, 1, 1];
            cdlChan.DelayProfile = 'CDL-C';
            cdlChan.DelaySpread = 100e-9;
            cdlChan.CarrierFrequency = obj.center_frequency;
            cdlChan.MaximumDopplerShift = (5 / 3.6) / physconst('lightspeed') * obj.center_frequency; % 5km/h pedestrian
            
            obj.CDLChannel = cdlChan;
            
        end
        
        function [rxWaveforms] = stepImpl(obj, txWaveforms, distance_3d)
            [rxWaveforms]       = obj.CDLChannel(txWaveforms);
            [~, rxWaveforms]    = apply_pathloss(obj, distance_3d, rxWaveforms);
            rxWaveforms         = obj.AWGNChannel(rxWaveforms);
        end
        
        function resetImpl(obj)
            
        end
        
        
    end
    
    methods(Access = private)
        function [pl_db, attenuated_waveforms] = apply_pathloss(obj, distance_3d, txWaveforms)
            pl_los = 32.4 + 17.3 * log10(distance_3d) + 20 * log10(obj.center_frequency);
            if obj.isLoS
                pl_db = pl_los;
            else %NLoS
                pl_nlos = 38.3 * log10(distance_3d) + 17.3 + 24.9 * log10(obj.center_frequency);
                pl_db = max(pl_los, pl_nlos);
            end
            
            pl_lin = db2pow(pl_db);
            attenuated_waveforms = txWaveforms ./ pl_lin;
            
        end
    end
end
