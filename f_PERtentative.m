function finalSet = f_PERtentative(candSet)
    % Create random PER values
    PER = rand(size(candSet));
    % Random threshold to determine whether or not a packets has been
    % received correctly
    threshold = 0.5;
    % Create Set of users to which we have transmitted AND they have
    % received it successfully
    finalSet = candSet(PER>threshold);
end