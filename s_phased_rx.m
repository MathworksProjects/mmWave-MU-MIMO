classdef s_phased_rx < matlab.System
    % s_phased_rx Implements a relative comprehensive phased array receiver
    
    properties(Nontunable)
        %% Center frequency
        center_frequency = 60.48e9;
        
        %% Rx gain
        rxGain = 1;
        
        %% number of receive antenna
        numRxElements = 1;
        
        %% Vis
        visualization = false
        
    end

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
            release(obj.collector);
            release(obj.receiver);
        end
    end
end
