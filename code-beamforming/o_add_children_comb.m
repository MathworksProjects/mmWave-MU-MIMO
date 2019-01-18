function [Set_ret,obj_vals_ret] = o_add_children_comb(comb,maxN,PRx,I,Noise,MaxObjF,MinObjF,Set,objs_val,orderedIndices,ObjFunc,avoid_rep)  %#ok
% o_add_children_comb - 
%
% Syntax:  [Set_ret,obj_vals_ret] = o_add_children_comb(comb,maxN,PRx,I,Noise,MaxObjF,MinObjF,Set,objs_val,orderedIndices,ObjFunc,avoid_rep)  %#ok
%
% Inputs:
%    comb - 
%    maxN - 
%    PRx - 
%    I - 
%    Noise - 
%    MaxObjF - 
%    MinObjF - 
%    Set - 
%    objs_val - 
%    orderedIndices - 
%    ObjFunc - 
%    avoid_rep - 
%
% Outputs:
%    Set_ret -Updated flow structure
%    obj_vals_ret -Updated flow structure
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: main
%
%------------- BEGIN CODE --------------
Set_ret = Set;
obj_vals_ret = objs_val;
for n=1:length(comb)
    if avoid_rep && n ~= 1 && comb(n-1) ~= 0 % We avoid including repeated children
        return
    end
    maxNarray = ones(length(comb),1)*maxN;
    tmp = o_sumOneToCombPos(comb,maxNarray,n);
    if ~isempty(tmp)
        selection = o_compute_selection(tmp,orderedIndices);
        obj_vals_ret = [obj_vals_ret;eval([ObjFunc,'(selection,PRx,I,Noise,MaxObjF,MinObjF)'])];
        Set_ret = [Set_ret,tmp]; % Column vector!!
    end
end



% EOF