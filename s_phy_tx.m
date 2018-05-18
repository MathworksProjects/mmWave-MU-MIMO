classdef s_phy_tx < matlab.System
    % phy_tx Add summary here
    %
    % NOTE: When renaming the class name untitled, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.

    properties(Nontunable)
        %% Waveform specification
        MCS = 1;
        TrainingLength = 4;
        PacketType = 'TRN-R';
        PSDULength = 1000; % in Bytes
    end

    properties(DiscreteState)

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
