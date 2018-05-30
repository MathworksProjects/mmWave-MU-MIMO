classdef s_phased_channel < matlab.System
    % s_phased_channel Simulate 3GPP TR 38.901 channel + pathloss + noise
    
    properties(Nontunable)
        numInputElements_row = 8;
        numInputElements_col = 8;
        numOutputElements_col = 1;
        numOutputElements_row = 1;
        SNR = 5;
        center_frequency = 60.48e9;
        applyPathLoss = false;
        profile = 'CDL-C';
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
            cdlChan.DelayProfile = obj.profile;
            cdlChan.DelaySpread = 100e-9;
            cdlChan.CarrierFrequency = obj.center_frequency;
            cdlChan.MaximumDopplerShift = (0 / 3.6) / physconst('lightspeed') * obj.center_frequency; % 0km/h pedestrian
            
            obj.CDLChannel = cdlChan;
            
        end
        
        function [rxWaveforms] = stepImpl(obj, txWaveforms, distance_3d, distance_2d)
            [rxWaveforms]       = obj.CDLChannel(txWaveforms);
            if obj.applyPathLoss
                [~, rxWaveforms]= apply_pathloss(obj, rxWaveforms, distance_3d, distance_2d);
            end
            rxWaveforms         = obj.AWGNChannel(rxWaveforms);
        end
        
        function resetImpl(obj)
            release(obj.CDLChannel);
            release(obj.AWGNChannel);
        end
    end
    
    methods(Access = private)
        function [pl_db, attenuated_waveforms] = apply_pathloss(obj, txWaveforms, distance_3d, distance_2d)
            % Calculate LoS probability threshold (indoor profile)
            prob_los_threshold = los_probability(distance_2d);
            rand_number = rand();
            % If the random number is less than the threshold, LoS
            pl_los = 32.4 + 17.3 * log10(distance_3d) + 20 * log10(obj.center_frequency) + lognrnd(0, 3);
            if rand_number < prob_los_threshold
                pl_db = pl_los;
            else % NLoS
                pl_nlos = 38.3 * log10(distance_3d) + 17.3 + 24.9 * log10(obj.center_frequency) + lognrnd(0, 8.03);
                pl_db = max(pl_los, pl_nlos);
            end
            pl_lin = db2pow(pl_db);
            attenuated_waveforms = txWaveforms ./ pl_lin;
        end
        
        function p_los = los_probability(distance_2d) % Refer to section 7.6.3.3 in TR 38.901
            if distance_2d <= 1.2
                p_los = 1;
            elseif distance_2d > 1.2 && (distance_2d < 6.5)
                p_los = exp(-(distance_2d - 1.2)/4.7);
            else
                p_los = exp(-(distance_2d - 6.5)/32.6) * 0.32;
            end
        end
    end
end
