classdef s_phased_channel_handle_version < matlab.System
    % s_phased_channel_handle_version Simulate 3GPP TR 38.901 channel + 
    % pathloss + noise, the difference is this system object takes in
    % handles for nr5gCDLChannel
    
    properties(Nontunable)
        SNR = 5;
        isLoS = true;
        applyPathLoss = false;
        NoiseFigure = -104.5;
    end
    
    properties(Access = private)

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
            
        end
        
        function [rxWaveforms] = stepImpl(obj, txWaveforms, distance_3d, cdlChanHandle)
            rxWaveforms = cdlChanHandle(txWaveforms);
            if obj.applyPathLoss
                [~, rxWaveforms]= apply_pathloss(obj, distance_3d, rxWaveforms);
            end
            rxWaveforms         = rxWaveforms + sqrt(db2pow(obj.NoiseFigure)) * randn(size(rxWaveforms));
            rxWaveforms         = obj.AWGNChannel(rxWaveforms);
        end
        
        function resetImpl(obj)
            release(obj.AWGNChannel);
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
