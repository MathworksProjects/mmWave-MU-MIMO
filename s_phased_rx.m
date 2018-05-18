classdef s_phased_rx < matlab.System
    % untitled4 Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.
    
    % Public, tunable properties
    properties
        
    end
    
    properties(Nontunable)
        %% Center frequency
        center_frequency = 60.48e9;
        %% Rx gain
        rxGain = 1;           % in dB
        % Noise temperature
        %         nT = 290;               % deg K
        %% number of receive antenna
        numRxElements = 1;
        
        %% Vis
        visualization = false
    end
    
    % Pre-computed constants
    properties(Access = private)
        antenna_array
        collector
        receiver
    end
    
    methods
        function obj = s_phased_rx(varargin)
            setProperties(obj,nargin,varargin{:})
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            light_speed = physconst('LightSpeed');
            lambda = light_speed / obj.center_frequency;
            antenna_spacing = lambda / 2;
            
            % Antenna
            if obj.numRxElements == 1
                obj.antenna_array = phased.IsotropicAntennaElement;
            else
                obj.antenna_array =  phased.ULA( obj.numRxElements, ...
                    'ElementSpacing',       antenna_spacing, ...
                    'Element',              phased.IsotropicAntennaElement('BackBaffled', true));
            end
            
            % Establish the collector
            obj.collector = phased.Collector(...
                'Sensor',               obj.antenna_array, ...
                'PropagationSpeed',     light_speed, ...
                'OperatingFrequency',   obj.center_frequency);
            
            obj.receiver = phased.ReceiverPreamp( ...
                'Gain',                 obj.rxGain); %, ...
            %                 'ReferenceTemperature', obj.nT, ...
            %                 'SampleRate',           gc.symRate*gc.upSample);
        end
        
        function rxSymbols = stepImpl(obj, txWaveforms, angle_from_tx)
            rxSymbols = obj.collector(txWaveforms, angle_from_tx);
            rxSymbols = obj.receiver(rxSymbols);
            
            %% Visualize
            if obj.visualization
                array_resp = phased.ArrayResponse( ...
                    'SensorArray', obj.antenna_array);
                scanAz = -180 : 180;
                max_rho = 5e3;
                arrayResponse = array_resp(obj.center_frequency, scanAz);
                figure(2)
                polar(deg2rad(scanAz(:)),abs(arrayResponse(:))*max_rho/4, '-m');
            end
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end
