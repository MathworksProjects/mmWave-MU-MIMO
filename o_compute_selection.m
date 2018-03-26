function [sel] = o_compute_selection(combination,orderedIndices)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    nRows = length(combination);
    sel_str = '[';
    for n=1:nRows-1
        if combination(n) ~= 0
            sel_str = [sel_str,'orderedIndices(',sprintf('%d',n),',',sprintf('%d',combination(n)),');'];
        else
            sel_str = [sel_str,'orderedIndices(',sprintf('%d',n),',',sprintf('%d',1),');'];
        end
    end
    if combination(nRows) ~= 0
        sel_str = [sel_str,'orderedIndices(',sprintf('%d',nRows),',',sprintf('%d',combination(nRows)),')];'];
    else
        sel_str = [sel_str,'orderedIndices(',sprintf('%d',nRows),',',sprintf('%d',1),')];'];
    end
    sel = eval(sel_str);
    for n=1:nRows
        if combination(n) == 0
            sel(n) = 0;
        end
    end
end

