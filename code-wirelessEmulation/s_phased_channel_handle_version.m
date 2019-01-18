classdef s_phased_channel_handle_version < matlab.System
    % s_phased_channel_handle_version Simulate 3GPP TR 38.901 channel + 
    % pathloss + noise, the difference is this system object takes in
    % handles for nr5gCDLChannel
    
    properties(Nontunable)
        SNR = 5;
        isLoS = true;
        applyPathLoss = false;
        noiseFigure = -104.5;
        center_frequency = 60.48e9;
    end
    
    properties(Access = private)

    end
    
    methods
        % Constructor
        function obj = s_phased_channel_handle_version(varargin)
            setProperties(obj,nargin,varargin{:})
        end
    end
    
    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            
        end
        
        function [rxWaveforms] = stepImpl(obj, txWaveforms, cdlChanHandle, distance_3d, distance_2d)
            rxWaveforms = cdlChanHandle(txWaveforms);
            if obj.applyPathLoss
                [~, rxWaveforms] = apply_pathloss_deterministic(obj, rxWaveforms, distance_3d);
            end
            rxWaveforms         = rxWaveforms + sqrt(db2pow(obj.noiseFigure)) * randn(size(rxWaveforms));
        end
        
        function resetImpl(obj)
            
        end

    end
    
    methods(Access = private)
        function [pl_db, attenuated_waveforms] = apply_pathloss_random(obj, txWaveforms, distance_3d, distance_2d)
            % Calculate LoS probability threshold (indoor profile)
            prob_los_threshold = los_probability(distance_2d);
            rand_number = rand();
            % If the random number is less than the threshold => LoS
            pl_los = 32.4 + 17.3 * log10(distance_3d) + 20 * log10(obj.center_frequency / 1e9) + lognrnd(0, 3);
            if rand_number < prob_los_threshold
                pl_db = pl_los;
            else % NLoS
                pl_nlos = 38.3 * log10(distance_3d) + 17.3 + 24.9 * log10(obj.center_frequency / 1e9) + lognrnd(0, 8.03);
                pl_db = max(pl_los, pl_nlos);
            end
            pl_lin = db2pow(pl_db);
            attenuated_waveforms = txWaveforms ./ pl_lin;
        end
        
        function [pl_db, attenuated_waveforms] = apply_pathloss_deterministic(obj, txWaveforms, distance_3d)
            pl_los = 32.4 + 17.3 * log10(distance_3d) + 20 * log10(obj.center_frequency / 1e9) + lognrnd(0, 3);
            if obj.isLoS
                pl_db = pl_los;
            else % NLoS
                pl_nlos = 38.3 * log10(distance_3d) + 17.3 + 24.9 * log10(obj.center_frequency / 1e9) + lognrnd(0, 8.03);
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
