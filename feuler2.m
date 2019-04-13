function [PPP, t_tot]=feuler2(x, y, z, Hx, Hy, Hz,Nt, pui, ee)

% Constantes &  Paramètres
t4 = tic;

rho = 8.9e3;    % densité du cobalt
cp = 420;       % capacité thermique du cobalt
K = 100;        % conductivité thermique du cobalt
E = 10;         % coefficient de transfert thermique air-cobalt (sans rotation du disque)
c = cp*rho;     % constante c=cp*rho
a = c/K;        % constante alpha = cp*rho/K

Tc = 300;                   % Température Cobalt (300 K)
Tair = 300;                 % Température Air

%%% PARAMÈTRES %%%%%%%%%%
% Temps
t0 = 0e-9; % Temps initial
Lt = 20e-9;  dt = Lt/Nt;   %Pas de temps

Nx = length(x); Ny=length(y); Nz=length(z);
N = Nx*Ny*Nz;
[X,Y,Z] = meshgrid(y,x,z); %Attention! il faut inverser x et y à cause de meshgrid

%% Initialisation des matrices et vecteurs

%statiques:
T0 = ones(N,1)*Tc;    %vecteur Nx1 : Distribution de température initiale
M = sparse(N,N);      %sparse NxN identity matrix
A = sparse(N,N);      %sparse NxN zeros matrix
b0 = zeros(N,1);      %vecteur Nx1 : partie statique du vecteur bn
I = speye(N,N);       %matrice qui multiplie le vecteur bn (éléments nuls pour x,y,z=Nx,Ny,Nz)
%dynamiques:
vect = zeros(0,Nt);   %vecteur des temps tn
PPP = zeros(N, Nt);   %vecteur avec la solution (température) à tous les temps tn


% Définition de la matrice des coefficients
t2=tic;
for k = 1:Nz
    for j=1:Ny
        for i=1:Nx
            
            % remplir la ligne l :
            l = (k-1)*Nx*Ny + (j-1)*Nx + i;
            
            % 1. Conditions de Diriclet -------------
            if (i==Nx)
                % noeud à la surface interne  x=Lx
                A(l,l)=1; b0(l) = Tc; I(l,l)=0;
                
            elseif (j==Ny)
                % noeud à la surface interne  y=Ly
                A(l,l)=1; b0(l) = Tc; I(l,l)=0;
                
            elseif (k==Nz)
                % noeud à la surface interne  z=Lz
                A(l,l)=1; b0(l) = Tc; I(l,l)=0;
                % ----------------------------------------
                
                % ----------------------------------------
            else
                % Contribution du point i,j,k
                
                hx = 1/(Hx(i+1)*Hx(i)); hy = 1/(Hy(j+1)*Hy(j)); hz = 1/(Hz(k+1)*Hz(k));
                A(l,l) = 1 + 2*dt*a^-1*( hx + hy + hz);
                M(l,l) = 1;
                % La contribution du terme source sera initialisée dans la résolution
                % temporelle
                
                % ----------------------------------------
                if (i==1)  %réflexion en x=0
                    
                    hx1 = Hx(i); A(l,l) = A(l,l) - dt*a^-1*hx1^-2*(4/3);
                    % Contribution suplémentaire du point i+1,j,k
                    c = index(i+1,j,k,Nx,Ny); A(l,c) = -1*dt*a^-1*hx1^-2*(1 -1/3);
                else
                    % Contribution du point i-1,j,k
                    c = index(i-1,j,k,Nx,Ny); hx=2/(Hx(i)*(Hx(i+1)+Hx(i)));
                    A(l,c) = -1*dt*a^-1*hx;
                    % Contribution du point i+1,j,k
                    c = index(i+1,j,k,Nx,Ny); hx=2/(Hx(i+1)*(Hx(i+1)+Hx(i)));
                    A(l,c) = -1*dt*a^-1*hx;
                end
                
                if (j==1) %réflexion en y=0
                    
                    hy1 = Hy(j); A(l,l) = A(l,l) - dt*a^-1*hy1^-2*(4/3);
                    % Contribution suplémentaire du point i,j+1,k
                    c = index(i,j+1,k,Nx,Ny); hy=2/(Hy(j+1)*(Hy(j+1)+Hy(j)));
                    A(l,c) = -1*dt*a^-1*hy*(1-1/3);
                else
                    % Contribution du point i,j-1,k
                    c = index(i,j-1,k,Nx,Ny); hy=2/(Hy(j)*(Hy(j+1)+Hy(j)));
                    A(l,c) = -1*dt*a^-1*hy;
                    % Contribution du point i,j+1,k
                    c = index(i,j+1,k,Nx,Ny); hy=2/(Hy(j+1)*(Hy(j+1)+Hy(j)));
                    A(l,c) = -1*dt*a^-1*hy;
                end
                
                if(k==1) %convection en z=0
                    
                    hz1 = Hz(k); A(l,l) = A(l,l) - dt*a^-1*hz1^-2*(4/3 - 2*E*hz1/(3*K));
                    b0(l) = dt*a^-1*( 2*E*Tair/(3*K*hz1));
                    % Contribution du point i,j,k+1
                    c = index(i,j,k+1,Nx,Ny); hz=2/(Hz(k+1)*(Hz(k+1)+Hz(k)));
                    A(l,c) = -1*dt*a^-1*hz*(1-1/3);
                else
                    % Contribution du point i,j,k-1
                    c = index(i,j,k-1,Nx,Ny); hz=2/(Hz(k)*(Hz(k+1)+Hz(k)));
                    A(l,c) = -1*dt*a^-1*hz;
                    
                    % Contribution du point i,j,k+1
                    c = index(i,j,k+1,Nx,Ny); hz=2/(Hz(k+1)*(Hz(k+1)+Hz(k)));
                    A(l,c) = -1*dt*a^-1*hz;
                end
            end
        end
    end
end
tmat = toc(t2);


%% Résolution dans le temps
bn = zeros([N,1]); Tn = T0;
Ainv = inv(A);
SSS = zeros(Nx,Ny,Nz,Nt); %Terme source aux temps tn

t3 =tic;
for n =0:Nt-1
    % Temps
    tn = t0 + dt*n;
    vect(n+1) = tn;
    
    % Définir bn (terme source) à chaque temps tn:
    S = K^-1*fSourceM(X,Y,Z,tn, pui);
    SSS(:,:,:,n+1) = S;
    vecS = I*reshape(S,[N,1]);
    bn = b0 + dt*a^-1*vecS;
    
    % Résolution
    b = M*Tn + bn;
    Tnplus1 = Ainv*b;
    PPP(:, n+1) = Tnplus1;
    Tn = Tnplus1;
end
tres=toc(t3);
t_tot=toc(t4);

%% Sauver les données du workspace :

% filename = sprintf ( 'NU_CS_%d_0%d.mat', N, round(ee*10) );
% save(filename, 'N', 'pui', 'vect', 'PPP', 'SSS', 'X', 'Y', 'Z',...
%     'Hx', 'Hy', 'Hz', 'tmat', 'tres', 't_tot', 'Nx', 'Ny', 'Nz', 'x', 'y', 'z')

end



