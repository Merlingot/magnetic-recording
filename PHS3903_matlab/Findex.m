function m = Findex(nx, ny, nz)
    global Nx Ny 
% index d'un noeud (ijk) dans une matrice 1D
    m = (nz-1)*Nx*Ny + (ny-1)*Ny + nx;
end

