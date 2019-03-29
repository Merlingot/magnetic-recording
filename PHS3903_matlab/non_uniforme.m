%%% Paramètres
global Nx Ny Nz
Nx = 3;
Ny = 3;
Nz = 3;
N = Nx*Ny*Nz;
Ta = 300; % Température ambiante (300 K)

H = zeros(1,N) + 1; %matrice des h !!! À définir

% Construction de la matrice du Laplacien
% Laplacien*T(r) -> L*T où L:matrice 2D, T:vecteur 1D

L=zeros(N,N);

%% Function handles

index = @ fIndex;
h = @ fh;

%% Définition de l'Opérateur L (LAPLACIEN) et de la matrice des h

for k = 1:Nz
    for j=1:Ny
        for i=1:Nx
            % remplir la ligne l de la matrice M
            l = (k-1)*Nx*Ny + (j-1)*Nx + i;
            
            % remplir la colone c : 
            if ((k>1)&&(k<Nz))&&((j>1)&&(j<Ny))&&((i>1)&&(i<Nx))
                % noeud qui est strictement à l'intérieur de la cellule de simulation
                
                % hi = x_(i+1,j,k) - x_(i,j,k) 
                % hi_ = x_(i,j,k) - x_(i-1,j,k)
                
                % À faire %%%%%
                hi = 1; hi_ = 1;
                hj = 1; hj_ = 1;
                hk = 1; hk_ = 1;
                %%%%%%%%%%%%%%%
                
                % Point i,j,k
                c = l;
                L(l,c) = -2*(hi*hi_*hj*hj_*hk_)^(-1);
                
                % Point i+1,j,k
                c = index(i+1,j,k);
                L(l,c) = 2*(hi*(hi+hi_))^(-1);
                % Point i,j+1,k
                c = index(i,j+1,k);
                L(l,c) = 2*(hj*(hj+hj_))^(-1);
                % Point i,j,k+1
                c = index(i,j,k+1);
                L(l,c) = 2*(hk*(hk+hk_))^(-1);
                
                % Point i-1,j,k
                c = index(i-1,j,k);
                L(l,c) = 2*(hi_*(hi+hi_))^(-1);
                % Point i,j-1,k
                c = index(i,j-1,k);
                L(l,c) = 2*(hj_*(hj+hj_))^(-1);
                % Point i,j,k-1
                c = index(i,j,k-1);
                L(l,c) = 2*(hk_*(hk+hk_))^(-1);
                
                
            elseif (i==1)
                % noeud sur le plan de symétrie x=0
                
                
                
                
            elseif (i==Nx)
                % noeud sur r x=Lx
                
            elseif (j==1)
                % noeud sur le plan de symétrie y=0
                
            elseif (j==Ny)
                % noeud à la surface interne  y=Ly
                
            elseif (k==1)
                % noeud à la surface externe  z=0 
                
            elseif (k==Nz)
                % noeud à la surface interne  z=Lz

            else
                disp('Erreur dans la définition de la matrice de coefficients');
            end
        end
    end
end




function m = fIndex(nx, ny, nz)
    global Nx Ny 
% index d'un point d'un matrice 3D dans une matrice 1D
    m = (nz-1)*Nx*Ny + (ny-1) + nx;
end

function h = fh(nx, ny, nz)
% valeur du pas de discrétisation
    m = (nz-1)*Nx*Ny + (ny-1) + nx;
    h = H(m) ; %indexation
end


