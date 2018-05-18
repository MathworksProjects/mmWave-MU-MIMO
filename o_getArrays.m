function array = o_getArrays(N,W,px,py,pz)
% Matrix W is a 2D matrix with dimensions [N,Nmax]
% If a given antenna is selected (W(.,.) ~= 0), we
% save its x,y,z location in the array vector
    Nmax = -Inf;
    for n = 1:N
        Nmax = max(Nmax,length(find(W(n,:))));
    end
    array = ones(3,Nmax,N)*(-Inf);
    for n = 1:N
        ind_array = 1;
        for p = 1:size(W,2)
            if not(W(n,p) == 0)
                rx = px(p);
                ry = py(p);
                rz = pz(p);
                array(:,ind_array,n) = [rx ry rz].';
                ind_array = ind_array + 1;
            end
        end
    end
end