function x_rnd = apply_random_modification_assignment(x,impact,structure,nAntennas,modifyAmpl,maxAmpl)
%APPLY_RANDOM_MODIFICATION_ASSIGNMENT - Apply random shuffling to antenna
%assignment
%This function applies a random shuffling affecting a number of antennas 
% assigned proportional to 'impact' and optionally modifying the amplitudes
% assigned, but always keeping the phase.
%
% Syntax:  x_rnd = apply_random_modification_assignment(x,impact,structure,nAntennas,modifyAmpl,maxAmpl)
%
% Inputs:
%    x           - 1-D array containing the antennas assignation, following
%                  the structure defined in structure
%    impact      - proportion of matrix cells that will be affected
%    structure   - String containing the atructure of the x variable (which
%                  position means what). For the moment, only
%                  'onlyAssigned' structure is supported:
%                  'onlyAssigned' - the last element is a real value
%                  indicating the baseband weight applied to all antennas.
%                  The rest of the array positions are divided into three
%                  equally sized parts: the first third contains the IDs of
%                  the antennas assigned, the second part contains its
%                  weights (intensity amplitude) and the third contains the
%                  phases applied.
%    nAntennas   - Integer representing the number of antennas available at
%                  the antenna array.
%    modifyAmpl  - boolean variable determining whether we can modify the
%                  value of the complex weights apart from applying random 
%                  shuffling
%    maxAmpl     - Maximum value for the weights' amplitudes
%
% Outputs:
%    x_rnd - Result of the modification
%
% Example: 
%    x_rnd = apply_random_modification(x,0.2,'onlyAssigned',24,true,1)
%    x_rnd = apply_random_modification(x,0.8,'onlyAssigned',12,false)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

%------------- BEGIN CODE --------------
    if strcmp(structure,'onlyAssigned')
        x_rnd = x;
        nAntennasAssigned = (size(x,2)-1)/3;
        nAntennasToModify = ceil(impact*nAntennasAssigned);
        antennasToModify = randperm(nAntennasAssigned,nAntennasToModify);
        antennasNotAssigned = setdiff(1:nAntennas,x(1:nAntennasAssigned));
        if ~isempty(antennasNotAssigned)
            antennasToAssign = datasample(antennasNotAssigned,nAntennasToModify);
            for i = 1:nAntennasToModify
                new_antenna = antennasToAssign(i);
                x_rnd(antennasToModify(i)) = new_antenna;
            end
        end

        if modifyAmpl
            amplitudesToModify = randperm(nAntennasAssigned,nAntennasToModify);
            mod = 2*rand(1,nAntennasToModify) - 1;
            x_rnd(amplitudesToModify) = max(0,min(maxAmpl,x_rnd(amplitudesToModify)+mod));
        end
    else
        fprintf('''%s'' structure is not implemented\n', structure);
        x_rnd = x;
    end
end