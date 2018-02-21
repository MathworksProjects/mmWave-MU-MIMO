function [ HPBW, SLL ] = helperFindLobes(array_handle, frequency, azimuth, elevation, dB_threshold, plot)
%3dB BEAMWIDTH Computation of the Half-Power Beam Width
%   We consider N elements, placed in (p_x, p_y, p_z), with weights w
%   px, py, px and w are row vectors
    if numel(azimuth) > 1 && numel(elevation) > 1
        error('SRM:beamwidth:BadArgs', ...
            'either azimuth or elevation must be a scalar value');
    end
    if numel(azimuth) == 1 && numel(elevation) == 1
        error('SRM:beamwidth:BadArgs', ...
            'either azimuth or elevation must be an array');
    end    
    
    if numel(azimuth) == 1
        [fieldval_el] = patternElevation(array_handle, frequency, azimuth, ...
        'Elevation', elevation, 'Type', 'Power');

        % Replacing any inf and zeros with a very low value
        fieldval_el(isinf(fieldval_el)) = -400;
        fieldval_el(fieldval_el == 0) = -400;
        fieldval_el = mag2db(fieldval_el);
        if plot
            figure;
            subplot(1,2,1);
            handle_polar_el = polarpattern(elevation,fieldval_el);
            handle_polar_el.AntennaMetrics = 1;
        end
        % Maximum value of the pattern
        [max_beam, max_index] = max(fieldval_el);
        curr_beam = max_beam;
        theta_ini = elevation(max_index);
        if max_index < numel(elevation)
            curr_index = max_index+1;
        else
            curr_index = 1;
        end
        theta = elevation(curr_index);
        omnidirectional_theta = false;
        increasing = false;
        decreasing = false;
        SLL = Inf;
        while (curr_beam > (max_beam-dB_threshold)) && ~(omnidirectional_theta)
            if curr_index < numel(elevation)
                curr_index = curr_index + 1;
            else
                curr_index = 1;
            end
            theta = elevation(curr_index);
            
            if curr_index == max_index
                omnidirectional_theta = true;
            end
            
            decreasing = (curr_beam-fieldval_el(curr_index)) > 0;
            if (increasing && decreasing)
                SLL = min(max_beam-curr_beam,SLL);
            end
            increasing = (curr_beam-fieldval_el(curr_index)) < 0;
            curr_beam = fieldval_el(curr_index);
        end
        while ~(omnidirectional_theta)
            if curr_index < numel(elevation)
                curr_index = curr_index + 1;
            else
                curr_index = 1;
            end
            
            if curr_index == max_index
                omnidirectional_theta = true;
            end
            
            decreasing = (curr_beam-fieldval_el(curr_index)) > 0;
            if (increasing && decreasing)
                SLL = min(SLL,max_beam-curr_beam);
            end
            increasing = (curr_beam-fieldval_el(curr_index)) < 0;
            curr_beam = fieldval_el(curr_index);
        end
        HPBW = 2*abs(theta_ini - theta);
    elseif numel(elevation) == 1
        [fieldval_az] = patternAzimuth(array_handle, frequency, elevation, ...
        'Azimuth', azimuth, 'Type', 'Power');

        % Replacing any inf and zeros with a very low value
        fieldval_az(isinf(fieldval_az)) = -400;
        fieldval_az(fieldval_az == 0) = -400;
        fieldval_az = mag2db(fieldval_az);
        if plot
            figure;
            subplot(1,2,2);
            handle_polar_az = polarpattern(azimuth,fieldval_az);
            handle_polar_az.AntennaMetrics = 1;
        end
        
        % Maximum value of the pattern
        [max_beam, max_index] = max(fieldval_az);
        curr_beam = max_beam;
        phi_ini = azimuth(max_index);
        if max_index < numel(azimuth)
            curr_index = max_index+1;
        else
            curr_index = 1;
        end
        phi = azimuth(curr_index);
        omnidirectional_phi = false;
        increasing = false;
        decreasing = false;
        SLL = Inf;
        while (curr_beam > (max_beam-dB_threshold)) && ~(omnidirectional_phi)
            if curr_index < numel(azimuth)
                curr_index = curr_index + 1;
            else
                curr_index = 1;
            end
            phi = azimuth(curr_index);
            
            if curr_index == max_index
                omnidirectional_phi = true;
            end
            
            decreasing = (curr_beam-fieldval_az(curr_index)) > 0;
            if (increasing && decreasing)
                SLL = min(max_beam-curr_beam,SLL);
            end
            increasing = (curr_beam-fieldval_az(curr_index)) < 0;
            curr_beam = fieldval_az(curr_index);
        end
        while ~(omnidirectional_phi)
            if curr_index < numel(azimuth)
                curr_index = curr_index + 1;
            else
                curr_index = 1;
            end
            
            if curr_index == max_index
                omnidirectional_phi = true;
            end
            decreasing = (curr_beam-fieldval_az(curr_index)) > 0;
            if (increasing && decreasing)
                SLL = min(max_beam-curr_beam,SLL);
            end
            increasing = (curr_beam-fieldval_az(curr_index)) < 0;
            curr_beam = fieldval_az(curr_index);
        end
        
        HPBW = 2*abs(phi_ini - phi);
    end
end