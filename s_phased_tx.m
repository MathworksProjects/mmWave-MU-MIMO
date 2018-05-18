classdef s_phased_tx < matlab.System
    % phased_tx Add summary here
    %
    % NOTE: When renaming the class name untitled5, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.
    
    properties
        
    end
    
    properties(Nontunable)
        %% Antenna configurations
        numTxElements_row = 8;
        numTxElements_col = 8;
        txPower = 8;    % Watts
        txGain = 0;     % dBW
        
        %% Waveform specifications
        center_frequency = 60.48e9;
        
        %% Vis
        visualization = false;
    end
    
    properties(DiscreteState)
        
    end
    
    % Pre-computed constants
    properties(Access = private)
        antenna_array;
        steeringvec;
        transmitter;
        radiator;
    end
    
    methods
        function obj = s_phased_tx(varargin)
            setProperties(obj,nargin,varargin{:})
        end
    end
    
    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Physical constants
            light_speed = physconst('light');
            wavelength = light_speed / obj.center_frequency;
            
            % Antenna definition
            %             obj.antenna_array = phased.ULA( ...
            %                 obj.numTxElements, ...
            %                 'ElementSpacing',           0.5 * wavelength, ...
            %                 'Element',                  phased.IsotropicAntennaElement('BackBaffled', true));
            
            obj.antenna_array = phased.URA( ...
                'Size',                     [obj.numTxElements_row, obj.numTxElements_col], ...
                'ElementSpacing',           0.5 * wavelength, ...
                'Element',                  phased.IsotropicAntennaElement('BackBaffled', true));
            
            obj.steeringvec = phased.SteeringVector( ...
                'SensorArray',              obj.antenna_array, ...
                'PropagationSpeed',         light_speed);
            
            %Gain per antenna element
            obj.transmitter = phased.Transmitter( ...
                'PeakPower',                obj.txPower / obj.numTxElements, ...
                'Gain',                     obj.txGain);
            
            %Transmit array
            obj.radiator = phased.Radiator( ...
                'Sensor',                   obj.antenna_array, ...
                'WeightsInputPort',         true, ...
                'PropagationSpeed',         light_speed, ...
                'OperatingFrequency',       obj.center_frequency,...
                'CombineRadiatedSignals',   false);
        end
        
        function txWaveforms = stepImpl(obj, txBits, toRxAngle, W)
            
            %% Temp-debug part -- 2 sub array beamformings
            if isempty(W) %% W unspecified
                % Hybrid beamforming -- 2 sub arrays
                numSubarray = 2;
                numSubElem = obj.numTxElements_col * obj.numTxElements_row / numSubarray;
                
                subpos = (-(numSubarray-1)/2 : (numSubarray-1)/2)*(0.5*numSubElem);
                subelempos = (-(numSubElem-1)/2 : (numSubElem-1)/2)*0.5;
                
                % complex weights used as part of digital baseband precoding
                wT_digital = steervec(subpos, [toRxAngle(1); 0]);
                
                % analog phase shift values used as part of RF precoding
                wT_analog = exp(1i*angle(steervec(subelempos, [toRxAngle(1); 0])));
                
                % From the system perspective, the effect of the hybrid beamforming can be
                % represented by hybrid weights as shown below.
                wT_hybrid = kron(wT_digital, wT_analog);
                
                % % Model signal traveling to mobile by applying a phaseshift on the
                % % elements since phased.Radiator does not do it when radiated signals are
                % % uncombined.
                % wR = obj.steeringvec(obj.center_frequency, -toRxAngle);
                % weight = wT_hybrid .* wR;
                
                
                % then radiate it -- since we uses weights to "mimic" the
                % boresight, we do not need to steer them physically
                W = wT_hybrid;
            end
            %% Radiate signals out of weights
            txWaveforms = obj.radiator(txBits, ...
                repmat([0; 0], 1, obj.numTxElements_row * obj.numTxElements_col), conj(W));
            
            %% Amplify those signals through transmitter
            txWaveforms = obj.transmitter(txWaveforms)

            if obj.visualization
                array_resp = phased.ArrayResponse( ...
                    'SensorArray', obj.antenna_array, ...
                    'WeightsInputPort', true);
                scanAz = -180 : 180;
                max_rho = 5e3;
                arrayResponse = array_resp(obj.center_frequency, scanAz, wT_hybrid);
                figure(1)
                polar(deg2rad(scanAz(:)),abs(arrayResponse(:))*max_rho/4, '-m');
                title('Transmitter Antenna Pattern');
            end
        end
        
        function resetImpl(obj)
            
        end
    end
end
