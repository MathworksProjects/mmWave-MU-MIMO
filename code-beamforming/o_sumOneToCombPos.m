function [newcomb] = o_sumOneToCombPos(combination, maxN, pos)
    if length(maxN) ~= length(combination)
        fprintf('maxN must have the same length as combination\n');
        return
    end
    newcomb = combination;
    if combination(pos) == maxN(pos)
    	newcomb = [];
    else
        newcomb(pos) = combination(pos) + 1;
    end
end

