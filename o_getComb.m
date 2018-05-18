function [nComb,comb,combMatrix,conf] = o_getComb(patch,dx,dy,array,conf)
    % Get array basic parameters
    array(:,array(1, :)== -Inf) = [];
    Nxy = size(array,2); % Number of antennas selected in this layer
    if Nxy ~= 0
        uPatchX = unique(patch(1,:)); % x coordinates in patch
        uPatchY = unique(patch(2,:)); % y coordinates in patch
        % Find Maximums - array and patch
        aMinX = min(array(1,:)); % min x coordinate among the antennas selected
        aMinY = min(array(2,:)); % min y coordinate among the antennas selected
        pMinX = min(patch(1,:)); % min x coordinate among all the antennas
        pMinY = min(patch(2,:)); % min x coordinate among all the antennas
        % Correction to the Origins
        corrX = aMinX - pMinX;
        corrY = aMinY - pMinY;
        % Corrected array to the Origins
        array(1,:) = array(1,:) - corrX; % displacement of the array pattern to the origin
        array(2,:) = array(2,:) - corrY;
        % Combinations
        amaxX = max(array(1,:));
        amaxY = max(array(2,:));
        movX = round((abs(amaxX - max(patch(1,:))))/dx);  
        movY = round((abs(amaxY - max(patch(2,:))))/dy);  
        combX = round(movX + 1); % # possible Combinations in x-axis
        combY = round(movY + 1); % # possible Combinations in y-axis
        nComb = combX*combY;  % Total combinations
        if conf.PlotDisplacements
            % Get the dimensions for the subplot
            for j = 1:nComb; if (j^2)>=nComb; sPtDim = j; break; end; end
        end
        comb = zeros(3,Nxy,nComb);
        combMatrix = repmat(patch,1,1,nComb);
        listX = (0:1:movX);%(0:dx:movX*dx);
        listY = (0:1:movY);%(0:dy:movY*dy);
        for idy = listY
            for idx = listX
                ID = idy*length(listX) + idx + 1;%round( (idy/dy)*length(listX) + idx/dx ) + 1;
                for a = 1:Nxy
                    comb(1,a,ID) = array(1,a) + idx*dx;
                    comb(2,a,ID) = array(2,a) + idy*dy;
                    [~,x] = min(abs(comb(1,a,ID)-uPatchX));
                    [~,y] = min(abs(comb(2,a,ID)-uPatchY));
                    colXOK = find(combMatrix(1,:,ID)==comb(1,a,ID));
                    colOK = find(combMatrix(2,colXOK,ID)==comb(2,a,ID));
                    %id = (y-1)*Nx + x; %Nxy -> 6
                    combMatrix(3,colXOK(combMatrix(2,colXOK,ID)==comb(2,a,ID)),ID) = 1;  % Assign 1 to variable z (assigned?)
                end
                if conf.PlotDisplacements
                    figure(conf.figIdx);
                    subplot(sPtDim,sPtDim,ID); hold on;
                    o_plotPatch(patch,conf.colorEmpty,50);
                    o_plotArray(comb(:,:,ID),'r',70);
                end
            end
        end
        conf.figIdx = conf.figIdx + 1;
    else
        nComb = 0;
        comb = zeros(3,Nxy,nComb);
        combMatrix = repmat(patch,1,1,nComb);
    end
end