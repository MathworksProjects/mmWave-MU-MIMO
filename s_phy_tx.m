classdef s_phy_tx < matlab.System
    % s_phy_tx Implements a 802.11ad phy transmitter

    properties(Nontunable)
        %% Waveform specification
        MCS = 1;
        TrainingLength = 4;
        PacketType = 'TRN-R';
        PSDULength = 1000; % in Bytes
    end

    properties(Access = private)
        dmg;
    end

    methods
        % Constructor
        function obj = s_phy_tx(varargin)
            setProperties(obj,nargin,varargin{:})
        end
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            obj.dmg = wlanDMGConfig( ...
                'MCS',              obj.MCS, ...
                'TrainingLength',   obj.TrainingLength, ...
                'PacketType',       obj.PacketType, ...
                'PSDULength',       obj.PSDULength);
        end

        function [txWaveform, cfgDMG] = stepImpl(obj, psdu)
            txWaveform = wlanWaveformGenerator(psdu, obj.dmg);
            cfgDMG = obj.dmg;
        end

        function resetImpl(obj)

        end
    end
end
