clear all
%% Constantes &  Paramètres
rho = 8.9e3;    % densité du cobalt
cp = 420;       % capacité thermique du cobalt

K = 100;        % conductivité thermique du cobalt
E = 10;         % coefficient de transfert thermique Air-Cobalt (sans rotation du disque)
c = cp*rho;     % constante c=cp*rho
a = c/K;        % constante alpha = cp*rho/k

Tc = 300;       % Température Cobalt (300 K)
Tair = 300;     % Température Air

%%% paramètres :
global Nx Ny
Nx = 20; Lx = 10e-9; hx = Lx/(Nx-1);
Ny = 20; Ly = 10e-9; hy = Ly/(Ny-1);
Nz = 20; Lz = 10e-9; hz = Lz/(Nz-1);
N = Nx*Ny*Nz;
Nt = 100;  dt = 100*a*hx^2; Lt = Nt*dt;  %Pas de temps

%%% valeurs initiales:
T0 = ones(N,1)*Tc;  % Distribution de température initiale
t0 = 0;          % Temps initial

%%% function handle :
index = @findex;

%% Initialisation des matrices M et A

H = zeros(N,N);      %sparse NxN zeros matrix - Termes de Maillage non uniforme
M = zeros(N,N);      %sparse NxN identity matrix - Ne change pas dans le temps
A = zeros(N,N);      %sparse NxN zeros matrix - Ne change pas dans le temps
b0 = zeros(N,1);
for k = 1:Nz
    for j=1:Ny
        for i=1:Nx
            % remplir la ligne l :
            l = (k-1)*Nx*Ny + (j-1)*Nx + i;
            
            % remplir la colone c :
            if ((k>1)&&(k<Nz))&&((j>1)&&(j<Ny))&&((i>1)&&(i<Nx))
                % noeud qui est strictement à l'intérieur de la cellule de simulation
                
                % Contribution du point i,j,k
                c = l; A(l,c) = 1 + 2*dt*a^-1*( hx^-2 + hy^-2 + hz^-2);
                M(l,c) = 1;
                
                % Contribution du point i+1,j,k
                c = index(i+1,j,k); A(l,c) = -1*dt*a^-1*hx^-2;
                % Contribution du point i,j+1,k
                c = index(i,j+1,k); A(l,c) = -1*dt*a^-1*hy^-2;
                % Contribution du point i,j,k+1
                c = index(i,j,k+1); A(l,c) = -1*dt*a^-1*hz^-2;
                
                % Contribution du point i-1,j,k
                c = index(i-1,j,k); A(l,c) = -1*dt*a^-1*hx^-2;
                % Contribution du point i,j-1,k
                c = index(i,j-1,k); A(l,c) = -1*dt*a^-1*hx^-2;
                % Contribution du point i,j,k-1
                c = index(i,j,k-1); A(l,c) = -1*dt*a^-1*hx^-2;
                            %Condition de Convection  --------------
                %
                
            %  Conditions de Neumann --------------
                %  0 = -3u1 + 4u2 - u3  @n+1
            elseif (i==1)
                % noeud sur le plan de symétrie x=0
                b0(l) = 0; M(l,l)=0;
                c = index(i,j,k); A(l,c) = -3;  %1,j,k
                c = index(i+1,j,k); A(l,c) = 4; %2,j,k
                c = index(i+2,j,k); A(l,c) = -1;%3,j,k
                
                
            elseif (j==1)
                % noeud sur le plan de symétrie y=0
                % Contribution du point i,j+1,k
                M(l,l) = 0; b0(l) = 0;
                c = index(i,j,k); A(l,c) = -3;   %i,1,k
                c = index(i,j+1,k); A(l,c) = 4;  %i,2,k
                c = index(i,j+2,k); A(l,c) = -1; %i,3,k
                
            % Conditions de Diriclet --------------
                % 0 = u_N@n+1 - Tc
            elseif (i==Nx)
                % noeud à la surface interne  x=Lx
                A(l,l)=1; b0(l) = Tc;
                
            elseif (j==Ny)
                % noeud à la surface interne  y=Ly
                A(l,l)=1; b0(l) = Tc;
                
            elseif (k==Nz)
                % noeud à la surface interne  z=Lz
                A(l,l)=1; b0(l) = Tc;
                       elseif (k==1)
                
                %noeud à la surface externe  z=0
%                 b0(l) = 2*E*Tair*hz;
%                 c = index(i,j,k); A(l,c) = -(2*E*hz + 3*K);   %i,1,k
%                 c = index(i,j,k+1); A(l,c) = 4*K;  %i,2,k
%                 c = index(i,j,k+2); A(l,c) = -1*K; %i,3,k

                % Contribution du point i,j,k
                c = l; A(l,c) = 1 + 2*dt*a^-1*( hx^-2 + hy^-2 + hz^-2);
                M(l,c) = 1;
                
                % Contribution du point i+1,j,k
                c = index(i+1,j,k); A(l,c) = -1*dt*a^-1*hx^-2;
                % Contribution du point i,j+1,k
                c = index(i,j+1,k); A(l,c) = -1*dt*a^-1*hy^-2;
                % Contribution du point i,j,k+1
                c = index(i,j,k+1); A(l,c) = -(1)*dt*a^-1*hz^-2;
                
                
                % Contribution du point i-1,j,k
                c = index(i-1,j,k); 
                A(l,c) = -1*dt*a^-1*hx^-2;
                % Contribution du point i,j-1,k
                c = index(i,j-1,k); A(l,c) = -1*dt*a^-1*hy^-2;
                % Contribution du point i,j,k-1 (condition de convection)
                c = index(i,j,k+1); A(l,c) = A(l,c) - dt*a^-1*hz^-2*4/(3 + E*hz/K) ;
                c = index(i,j,k+2); A(l,c) = A(l,c) - dt*a^-1*hz^-2*(-1)/(3 + E*hz/K) ;
                b0(l) =  dt*a^-1*hz^-2*(E*hz + Tair)/(3*K + E*hz); 
               
                
            end
        end
    end
end


%% Résolution dans le temps
bn = b0; Tn = T0;
Ainv = inv(A);

for n = [0:Nt]
    % Définir bn (terme source) à chaque temps tn:
    % Change seulement pour les noeuds qui sont strictement à l'intérieur de la cellule de simulation
    for k = 1:Nz-1
        for j=2:Ny-1
            for i=2:Nx-1
                % remplir la ligne l :
                l = (k-1)*Nx*Ny + (j-1)*Nx + i;
                bn(l) = b0(l) +  a^-1*K^-1*fSource([i*hx,j*hy,k*hz], t0 + dt*n);
            end
        end
    end
    
    b = M*Tn + bn; 
    Tnplus1 = Ainv*b;
    Tn = Tnplus1;
    PPP(:,n+1) = Tnplus1;
end




%% Visualisation



for  i = 1:Nt
    Tr = reshape(PPP(:,i),[Nx,Ny,Nz]);
    x = linspace(0,Lx,Nx);
    y = linspace(0,Ly,Ny);
    z = linspace(0,Lz,Nz);
    [X,Y,Z] = meshgrid(x,y,z);
    xslice = [0,Lx];
    yslice = [0,Ly];
    zslice = [0, Lz];
    slice(X,Y,Z,Tr,xslice,yslice,zslice)
    colorbar
    xlabel('x')
    ylabel('y')
    zlabel('z')
    dim = [0.08 0.08 0.025 0.025];
    t = t0 + (i-1)*Nt;
    str = ['t = ' num2str(t) ' s'];
    pause
%     annotation('textbox',dim,'String',str,'FitBoxToText','on');
end

%%
function m = findex(nx, ny, nz)
global Nx Ny
% index d'un noeud (ijk) dans une matrice 1D
m = (nz-1)*Nx*Ny + (ny-1)*Nx + nx;
end
