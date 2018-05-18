function [newcomb] = o_sumOneToCombination(combination, maxN)
    s = size(combination,1); % Suponemos vector columna
    if length(maxN) ~= s
        fprintf('maxN must have the same length as combination\n');
        return
    end
    newcomb = combination;
    for i = s:-1:1
        if combination(i) == maxN(i)
            newcomb(i) = 0;
            if i == 1
                return
            end
        else
            newcomb(i) = combination(i) + 1;
            return
        end
    end
end

