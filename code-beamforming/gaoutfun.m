function [state,options,optchanged] = gaoutfun(options,state,flag)
global bestS
optchanged = false;
switch flag
    case 'done'
        bestS = state.Best;
end
