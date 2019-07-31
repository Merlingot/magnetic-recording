function m = index(nx, ny, nz, Nx, Ny)
% index d'un noeud (ijk) dans une matrice 1D
m = (nz-1)*Nx*Ny + (ny-1)*Nx + nx;
end

