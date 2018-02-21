function patch = getPatch(Nx,Ny,px,py)
    patch = zeros(3,Nx*Ny);
    for idx = 1:1:Nx
        for idy = 1:1:Ny
            id = (idy-1)*Nx + idx;
            id_input = (idx-1)*Ny + idy;
            patch(:,id) = [px(id_input) ; py(id_input) ; 0];
        end
    end
end
