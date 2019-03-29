clear all
%% Constantes
rho = 8.9*1.0000e3;     %densité du cobalt
cp = 420;               %capacité thermique du cobalt

k = 100;        %Conductivité thermique du cobalt
xi=10;          % coefficient de transfert thermique Air-Cobalt (sans rotation du disque)
c = cp*rho;     % constante c=cp*rho
a = k/c;        % constante a = k/cp*rho

Tc = 300;       % Température Cobalt (300 K)
Tair = 300;     % Température Air

%% Paramètres 
global Nx Ny
Nx = 5; Lx = 10e-9; hx = 10e-9/(Nx-1);
Ny = 5; Ly = 10e-9; hy = 10e-9/(Ny-1);
Nz = 5; Lz = 10e-9; hz = 10e-9/(Nz-1);
N = Nx*Ny*Nz;
Nt = 10; Lt = 10e-8; dt = Lt/Nt;   %Pas de temps
%%%

T0 = ones(N,1)*Tc; %Distribution de température initiale

%% Initialisation des matrices M et A

M=eye(N,N);    %sparse NxN identity matrix - Ne change pas dans le temps
A=sparse(N,N); %sparse NxN zeros matrix - Ne change pas dans le temps
b0=zeros(N,1); %vecteur Nx1 - b(t=0) - Ne change pas dans le temps (b initial - utilisé une seule fois)

for k = 1:Nz
    for j=1:Ny
        for i=1:Nx
            % remplir la ligne l :
            l = (k-1)*Nx*Ny + (j-1)*Nx + i;
            
            % remplir la colone c : 
            if ((k>1)&&(k<Nz))&&((j>1)&&(j<Ny))&&((i>1)&&(i<Nx))
                % noeud qui est strictement à l'intérieur de la cellule de simulation
                
                % Contribution du point i,j,k
                c = l;
                A(l,c) = -2*( 1/hx^2 + 1/hy^2 + 1/hz^2);
%                 b0(l) =  fSource([i*hx,j*hy,k*hz],0); %% Ce terme va changer dans le temps
                
                % Contribution du point i+1,j,k
                c = index(i+1,j,k);
                A(l,c) = 1/hx^2;
                % Contribution du point i,j+1,k
                c = index(i,j+1,k);
                A(l,c) = 1/hy^2;
                % Contribution du point i,j,k+1
                c = index(i,j,k+1);
                A(l,c) = 1/hz^2;
                
                % Contribution du point i-1,j,k
                c = index(i-1,j,k);
                A(l,c) = 1/hx^2;
                % Contribution du point i,j-1,k
                c = index(i,j-1,k);
                A(l,c) = 1/hy^2;
                % Contribution du point i,j,k-1
                c = index(i,j,k-1);
                A(l,c) = 1/hz^2;
                
            %  Conditions de Neumann --------------  
            elseif (i==1)
                % noeud sur le plan de symétrie x=0
                
                M(l,l) = 0;
                c = index(i,j,k);   A(l,c) = -3; %noeud 1,j,k
                c = index(i+1,j,k); A(l,c) = 4;  %noeud 2,j,k
                c = index(i+2,j,k); A(l,c) = -1; %noeud 3,j,k
                
               
            elseif (j==1)
                % noeud sur le plan de symétrie y=0
                
                M(l,l) = 0;
                c = index(i,j,k);   A(l,c) = -3; %noeud i,1,k
                c = index(i,j+1,k); A(l,c) = 4;  %noeud i,2,k
                c = index(i,j+2,k); A(l,c) = -1; %noeud i,3,k
            
            
            % Conditions de Diriclet --------------
            elseif (i==Nx)
                % noeud à la surface interne  x=Lx
                M(l,l) = 0;
                A(l,l) = 2*hx; b0(l) = -2*Tc*hx; %noeud i,j,k
                
                
            elseif (j==Ny)
                % noeud à la surface interne  y=Ly
                M(l,l) = 0;
                A(l,l) = 2*hy; b0(l) = -2*Tc*hy; %noeud i,j,k

            elseif (k==Nz)
                % noeud à la surface interne  z=Lz
                M(l,l) = 0;
                A(l,l) = 2*hz; b0(l) = -2*Tc*hz; %noeud i,j,k
            
                
            %Condition de Convection  --------------  
            elseif (k==1)
                % noeud à la surface externe  z=0 
                M(l,l) = 0;
                c = index(i,j,k);   A(l,c) = -2*xi*hz - 3*k; %noeud i,j,1
                c = index(i,j,k+1); A(l,c) = 4;  %noeud i,j,2
                c = index(i,j,k+2); A(l,c) = -1; %noeud i,j,3
                b0(l) = 2*Tair*hz;

            else
                disp('Erreur dans la définition de la matrice de coefficients');
            end
        end
    end
end

if (det(M) == 0)
    disp('Matrice M est singulière')
end

%% Résolution dans le temps
bn = b0; Tn = T0;

for n = [0:Nt]  
    % définir bn (terme source) à chaque temps tn:
    % Change seulement pour les noeuds qui sont strictement à l'intérieur de la cellule de simulation  
    for k = 2:Nz-1
        for j=2:Ny-1
            for i=1:Nx
                % remplir la ligne l :
                l = (k-1)*Nx*Ny + (j-1)*Nx + i;
                bn(l) =  fSource([i*hx,j*hy,k*hz],dt*n); 
            end
        end
    end
   
    % bbn:
    bbn = (M + a*dt*A)*Tn + bn*dt/c ;
    
    % Tn+1 = bbn/M
    [L,U]=lu(M);Tnplus1=U\(L\bbn);  
end    

Tr=reshape(Tnplus1,Nx,Ny,Nz);



