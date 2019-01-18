function finalSet = f_PERtentative(candSet,threshold)
% F_PERTENTATIVE - The function evaluates the chances the packets transmited
% to the users in 'candSet' have to be succesfully received. The quality of
% the link is controlled by 'threshold', which represents the PER. The
% function return the subset of 'candSet' that statistically see their
% packets get through the channel in the current transmission.
%
% Syntax:  finalSet = f_PERtentative(candSet,threshold)
%
% Inputs:
%    candSet - Vector containing a list of user-IDs that are being sent a
%              packet in the current slot.
%    threshold - Vector containing the PER of each user. Its size is same
%                as candSet
%
% Outputs:
%    finalSet - Vectpr of the list of user-IDs whose packets are received
%               succesfully at the receiver
%
% Example: 
%    candSet = [1 2 3 4 5];
%    threshold = [0.1 0.2 0.3 0.7 0.4];
%    finalSet = f_PERtentative(candSet,threshold)
%    disp('finalSet')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: f_PER,  main,  main_runnable
%
%------------- BEGIN CODE --------------

% Create random value uniformly distributed between 0 and 1
PER = rand(size(candSet));
% Create Set of users to which we have transmitted AND they have
% received it successfully
finalSet = candSet(PER>threshold);

% [EOF]