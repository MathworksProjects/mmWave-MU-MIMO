function [max_obj, max_sel] = o_find_best_combination(problem,conf,PRx,I,patch,comb,combMat,arraysCell)
%% Algorithm 2, max sum(SNR)
% Order arrays into a 3D-matrices cell: size:3, Nmax, N
%[~,orderedIndices] = sort(PRx,2,'descend');
[~,orderedIndices] = sort(PRx-10*log10(sum(10.^(I/10),3)),2,'descend');
%PRx(user,orderedIndices(user,x)) is the PRx of the xth solution ordered by PRx
%PRx(user,orderedIndicesSNR(user,x)) is the PRx of the xth solution ordered by PRx'SNR'

sol_found = false;
init_comb = zeros(problem.nUsers,1);
selection = o_compute_selection(init_comb,orderedIndices);
max_obj = -Inf;
first_obj = 0;
max_comb = [];
max_sel = [];
first_comb = [];
count = 1;
count_feas = 0;
count_newmax = 0;
Set = [];
obj_vals = [];
[Set,obj_vals] = o_add_children_comb(init_comb,length(problem.NmaxArray),PRx,I,...
    problem.Noise,problem.MaxThr,problem.MinThr,Set,obj_vals,orderedIndices,...
    conf.ObjFunc,true);
[~,SetOrderingIndObj] = sort(obj_vals,'descend');

cut_condition = false;
if conf.NoCut
    cutConditionEval = 'false';
elseif strcmp(conf.ObjFunc,'compute_sumSNR')
    cutConditionEval = '((10*log10(sum(PRx_lin))-problem.Noise) <= max_obj)';
elseif strcmp(conf.ObjFunc,'compute_weightedSNR')
    cutConditionEval = '(problem.nUsers*(10*log10(sum(PRx_lin))-problem.Noise) <= max_obj)';
elseif strcmp(conf.ObjFunc,'compute_averageSNR')
    cutConditionEval = '((10*log10(sum(PRx_lin))-problem.Noise) <= problem.nUsers*max_obj)';
elseif strcmp(conf.ObjFunc,'compute_sumCap')
    cutConditionEval = '(sum(log2(1+(PRx_lin/(10^(problem.Noise/10))))) <= max_obj)';
elseif strcmp(conf.ObjFunc,'compute_averageCap')
    cutConditionEval = '(sum(log2(1+(PRx_lin/(10^(problem.Noise/10))))) <= problem.nUsers*max_obj)';
elseif strcmp(conf.ObjFunc,'o_compute_avCap_maxminthr_LAC')
    cutConditionEval = '(sum(log2(1+(PRx_lin/(10^(problem.Noise/10))))) <= problem.nUsers*max_obj)';
elseif strcmp(conf.ObjFunc,'compute_averageSumCap')
    cutConditionEval = '(log2(1+sum(PRx_lin/(10^(problem.Noise/10)))) < max_obj)'; % No claro que funcione bien
else
    disp('ERROR: Unknown Object Maximization Function for Algorithm 2');
    return
end
if conf.verbosity >= 1
    fprintf('#Iter.\t#New max.\t#Feas. sols.\tObj.Val.\tAnt. used\n');
    fprintf('=========================================================\n');
end
while((~sol_found && ~isempty(Set)) || (sol_found && ~isempty(Set) && ~cut_condition))
    count = count + 1;
    if conf.verbosity >= 1
        fprintf('%d\t\t%d\t\t%d\t\t%f\t\t',count,count_newmax,count_feas,max_obj);
        for user = 1:problem.nUsers
            % The following computation obtains the number of subarrays
            % / antennas actually used by the user
            if selection(user) ~= 0
                fprintf('%d ', ...
                    sum(arraysCell{selection(user)}(1,:,user) ~= -Inf));
            else
                fprintf('0 ');
            end
        end
        fprintf('\n');
    end
    firstEl = SetOrderingIndObj(1);
    curr_comb = Set(:,firstEl);
    selection = o_compute_selection(curr_comb,orderedIndices);
    
    if ~conf.feasibility || (obj_vals(firstEl) ~= -Inf && ...
        o_isFeasComb(patch,comb,combMat,selection,find(selection ~= 0)))
        count_feas = count_feas + 1;
        if conf.verbosity > 1
            disp('New feasible combination found')
            fprintf('Antennas used per user: ');
            for user = 1:problem.nUsers
                % The following computation obtains the number of subarrays
                % / antennas actually used by the user
                if selection(user) ~= 0
                    fprintf('%d\t', ...
                        sum(arraysCell{selection(user)}(1,:,user) ~= -Inf));
                else
                    fprintf('0\t');
                end
            end
            fprintf('\n');
            fprintf('Obj. value: %f\n\n',obj_vals(firstEl));
        end
        if ~sol_found
%             disp('First feasible combination found!')
%             fprintf('Antennas used per user: ');
%             for user = 1:problem.nUsers
%                 % The following computation obtains the number of subarrays
%                 % / antennas actually used by the user
%                 if selection(user) ~= 0
%                     fprintf('%d\t', ...
%                         sum(arraysCell{selection(user)}(1,:,user) ~= -Inf));
%                 else
%                     fprintf('0\t');
%                 end
%             end
%             fprintf('\n');
%             %display(selection');
%             fprintf('Obj. value: %f\n\n',obj_vals(firstEl));
            first_obj = obj_vals(firstEl);
            first_sel = selection;
            first_comb = curr_comb;
            sol_found = true;
        end
        if max_obj < obj_vals(firstEl)
            count_newmax = count_newmax + 1;
            max_obj = obj_vals(firstEl);
            max_sel = selection;
            max_comb = curr_comb;            
%             disp('New max found!');fprintf('Antennas used per user: ');
%             for user = 1:problem.nUsers
%                 % The following computation obtains the number of subarrays
%                 % / antennas actually used by the user
%                 if selection(user) ~= 0
%                     fprintf('%d\t', ...
%                         sum(arraysCell{selection(user)}(1,:,user) ~= -Inf));
%                 else
%                     fprintf('0\t');
%                 end
%             end
%             fprintf('\n');
%             fprintf('Obj. value: %f\n\n',obj_vals(firstEl));
        end
    end
    PRx_lin = zeros(1,problem.nUsers);
    for u=1:problem.nUsers
        if (selection(u) ~= 0)
            PRx_lin(u) = (10^(PRx(u,selection(u))/10));
        end
    end
    cut_condition = eval(cutConditionEval);
    obj_vals(firstEl) = [];
    Set(:,firstEl) = [];
    [Set,obj_vals] = o_add_children_comb(curr_comb,length(problem.NmaxArray),PRx,I,...
        problem.Noise,problem.MaxThr,problem.MinThr,Set,obj_vals,orderedIndices,...
        conf.ObjFunc, true);
    [~,SetOrderingIndObj] = sort(obj_vals,'descend');
end

if conf.verbosity >= 1
    fprintf('%d\t\t%d\t\t%d\t\t%f\t\t',count,count_newmax,count_feas,max_obj);
    for user = 1:problem.nUsers
        % The following computation obtains the number of subarrays
        % / antennas actually used by the user
        if selection(user) ~= 0
            fprintf('%d ', ...
                sum(arraysCell{selection(user)}(1,:,user) ~= -Inf));
        else
            fprintf('0 ');
        end
    end
    fprintf('\n\n');
end

if conf.verbosity >= 1 && ~sol_found
    fprintf('No solution found!\n\n');
end

if conf.verbosity >= 1
    fprintf('# Iterations: %d\t', count);
    fprintf('# New maxima: %d\t', count_newmax);
    fprintf('# Feasible solutions: %d\n\n', count_feas);ç
end

end