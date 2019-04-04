clear all
%% Constantes &  Paramètres
rho = 8.9e3;    % densité du cobalt
cp = 420;       % capacité thermique du cobalt

K = 100;        % conductivité thermique du cobalt
E = 10;         % coefficient de transfert thermique air-cobalt (sans rotation du disque)
c = cp*rho;     % constante c=cp*rho
a = c/K;        % constante alpha = cp*rho/K

Tc = 300;                   % Température Cobalt (300 K)
Tair = 300;                 % Température Air
TCurie = 320 + 375.15;      % Température de Curie du Cobalt (320 C)

%%% paramètres :
global Nx Ny
Nx = 20; Lx = 20e-9; hx = Lx/(Nx-1);
Ny = 20; Ly = 20e-9; hy = Ly/(Ny-1);
Nz = 20; Lz = 20e-9; hz = Lz/(Nz-1);
N = Nx*Ny*Nz;
Lt = 10e-9; Nt = 100;  dt = Lt/Nt;   %Pas de temps

%%% valeurs initiales:
T0 = ones(N,1)*Tc;  % Distribution de température initiale
t0 = 0e-9;          % Temps initial

%%% function handle :
index = @findex;

%% Initialisation des matrices M et A

% H = sparse(N,N);      %sparse NxN zeros matrix - Termes de Maillage non uniforme (À faire)
M = sparse(N,N);      %sparse NxN identity matrix - Ne change pas dans le temps
A = sparse(N,N);      %sparse NxN zeros matrix - Ne change pas dans le temps
b0 = zeros(N,1);      %vecteur 1xN 
vect = zeros(0,Nt);   %vecteur des temps tn
PPP = zeros(N, Nt);   %vecteur avec les solutions à tous les temps tn

tic
for k = 1:Nz
    for j=1:Ny
        for i=1:Nx
            % remplir la ligne l :
            l = (k-1)*Nx*Ny + (j-1)*Nx + i;
            
            % ---------------------------------------
            % 1. Conditions de Diriclet -------------
            if (i==Nx)
                % noeud à la surface interne  x=Lx
                A(l,l)=1; b0(l) = Tc;
                
            elseif (j==Ny)
                % noeud à la surface interne  y=Ly
                A(l,l)=1; b0(l) = Tc;
                
            elseif (k==Nz)
                % noeud à la surface interne  z=Lz
                A(l,l)=1; b0(l) = Tc;
                % ----------------------------------------
                
                % ----------------------------------------
            else
                    % Contribution du point i,j,k
                    A(l,l) = 1 + 2*dt*a^-1*( hx^-2 + hy^-2 + hz^-2);
                    M(l,l) = 1;
                    % La contribution du terme source sera initialisé dans la résolution
                    % temporelle

                    % ----------------------------------------
                    if (i==1)  %réflexion en x=0
                        A(l,l) = A(l,l) - dt*a^-1*hx^-2*(4/3);
                        % Contribution suplémentaire du point i+1,j,k
                        c = index(i+1,j,k); A(l,c) = -1*dt*a^-1*hx^-2*(1 -1/3);
                    else
                        % Contribution du point i-1,j,k
                        c = index(i-1,j,k); A(l,c) = -1*dt*a^-1*hx^-2;
                        % Contribution du point i+1,j,k
                        c = index(i+1,j,k); A(l,c) = -1*dt*a^-1*hx^-2;
                    end
                    
                    if (j==1) %réflexion en y=0
                        A(l,l) = A(l,l) - dt*a^-1*hy^-2*(4/3);
                        % Contribution suplémentaire du point i,j+1,k
                        c = index(i,j+1,k); A(l,c) = -1*dt*a^-1*hy^-2*(1-1/3);
                    else
                        % Contribution du point i,j-1,k
                        c = index(i,j-1,k); A(l,c) = -1*dt*a^-1*hy^-2;
                        % Contribution du point i,j+1,k
                        c = index(i,j+1,k); A(l,c) = -1*dt*a^-1*hy^-2;
                    end

                    if(k==1) %convection en z=0
                        A(l,l) = A(l,l) - dt*a^-1*hz^-2*(4/3 - 2*E*hz/(3*K));
                        b0(l) = dt*a^-1*hz^-2*( 2*E*hz*Tair/(3*K));
                        % Contribution du point i,j,k+1
                        c = index(i,j,k+1); A(l,c) = -1*dt*a^-1*hz^-2*(1-1/3);
                    else
                        % Contribution du point i,j,k-1
                        c = index(i,j,k-1); A(l,c) = -1*dt*a^-1*hx^-2;

                        % Contribution du point i,j,k+1
                        c = index(i,j,k+1); A(l,c) = -1*dt*a^-1*hz^-2;
                    end
                    
                
            end
        end
    end
end
t1 = toc;


%% Résolution dans le temps
bn = b0; Tn = T0;
Ainv = inv(A);
SSS = zeros(Nx,Ny,Nz, Nt);

tic
for n =0:Nt 
    
    % Définir bn (terme source) à chaque temps tn:
    % Change seulement pour certains noeuds
    
    tn = t0 + dt*n;
    for k = 1:Nz-1
        for j=1:Ny-1
            for i=1:Nx-1
                % remplir la ligne l :
                l = (k-1)*Nx*Ny + (j-1)*Nx + i;
                s = K^-1*fSource([i*hx,j*hy,k*hz], tn);
                bn(l) = b0(l) +  dt*a^-1*s;
                SSS(i,j,k, n+1) = s;
            end
        end
    end
    
    vect(n+1) = tn;
    b = M*Tn + bn; 
    Tnplus1 = Ainv*b;
    PPP(:, n+1) = Tnplus1; 
    Tn = Tnplus1;
end
t2=toc;



%% Visualisation
an = [];
dim = [0.08 0.08 0.025 0.025];

for  i = 1:Nt 
    Tr = reshape(PPP(:,i),[Nx,Ny,Nz]);
    Sr = SSS(:,:,:,i);
    
    x = linspace(0,Lx,Nx);
    y = linspace(0,Ly,Ny);
    z = linspace(0,Lz,Nz);
    [X,Y,Z] = meshgrid(x,y,z);
    xslice = [0,Lx];
    yslice = [0,Ly];
    zslice = [0, Lz];
    ti = vect(i);
    str = ['t = ' num2str(ti) ' s'];
    delete(an);
    
    
    subplot(1,3,1)
    slice(X,Y,Z,Tr,xslice,yslice,zslice)
    colorbar
    colormap('jet')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Température')
    an = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    subplot(1,3,2)
    Tbool = double(Tr>TCurie);
    slice(X,Y,Z,Tbool,xslice,yslice,zslice)
    colorbar
    colormap('jet')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Régions de température supérieure à la température de Curie')
    caxis([0 1])
    
    subplot(1,3,3)
    slice(X,Y,Z,Sr,[Lx],[Ly],[0])
    colorbar
    colormap('jet')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Terme source')
    
    pause %%Pour voir frame par frame
end

fprintf('fin')

%% Choses à améliorer
%1. Vecteur du terme source -> changer les boucles for pour des opérations
%   matricielles


%%
function m = findex(nx, ny, nz)
global Nx Ny
% index d'un noeud (ijk) dans une matrice 1D
m = (nz-1)*Nx*Ny + (ny-1)*Nx + nx;
end
