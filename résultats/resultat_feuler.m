function [pui, PPP, Tmax, theat, twrite, tcool, vect, X, Y, Z, xcurie, ycurie, zcurie]= resultat_feuler(Nt, Nx, Ny, Nz, pui, ee)

% Arguments:
%   Nt : nbre d'itérations dans le temps
%   Nx, Ny, Nz : nombre de noeuds maillage spatial
%   pui : puissance moyenne du laser (W)
%   ee : facteur d'étirement 

% Retour:
%   pui : puissance moyenne du laser (W)
%   PPP : distribution de température à chaque temps t (NxNt)
%   Tmax : Température max à chaque temps t
%   theat, twrite, tcool : temps d'écritures caractéristiques (s)
%   vect : temps t
%   X, Y, Z : meshgrid des coordonnées
%   xcurie, ycurie, zcurie : coordonées max où T>Tcurie



%%
% Constantes &  Paramètres
rho = 8.9e3;    % densité du cobalt
cp = 420;       % capacité thermique du cobalt
K = 100;        % conductivité thermique du cobalt
E = 10;         % coefficient de transfert thermique air-cobalt (sans rotation du disque)
c = cp*rho;     % constante c=cp*rho
a = c/K;        % constante alpha = cp*rho/K

Ldiff = (K*3e-9/(cp*rho))^(1/2);    %Longueur de diffusion

Tc = 300;                   % Température Cobalt (300 K)
Tair = 300;                 % Température Air
Tcurie = 320 + 273.15;

%%% PARAMÈTRES %%%%%%%%%%
% Temps
t0 = 0e-9; % Temps initial
Lt = 20e-9;  dt = Lt/Nt;   %Pas de temps
% Longueur du domaine
Lx = 3*Ldiff/2; Ly = 3*Ldiff/2; Lz = 3*Ldiff;

N = Nx*Ny*Nz;
%% Définition du Maillage non uniforme


%i=1:x=0, i=Nx:x=Lx
ex = ee; %paramètre de scaling
syms kx; sx = symsum( (1+ex)^(kx-2), kx, 2, Nx);
h0x = double(Lx/sx);
Hx = [h0x, h0x]; x = [0,h0x]; %Hx:pas en x, x:coordonées en x
himoins1 = h0x;ximoins1=h0x;
for i=3:Nx
    hi =(1+ex)*himoins1; Hx(i) = hi;
    xi = ximoins1 + hi; x(i) = xi;
    himoins1=hi; ximoins1 = xi;
end

%j=1:y=0, j=Ny:y=Ly
ey = ee; %paramètre de scaling
syms ky; sy = symsum( (1+ey)^(ky-2), ky, 2, Ny); h0y = double(Ly/sy);
Hy = [h0y, h0y]; y = [0,h0y];
himoins1 = h0y;yimoins1=h0y;
for i=3:Ny
    hi =(1+ey)*himoins1; Hy(i) = hi;
    yi = yimoins1 + hi; y(i) = yi;
    himoins1=hi; yimoins1 = yi;
end

%k=1:z=0, k=Nz:z=Lz
ez = ee; %paramètre de scaling
syms kz; sz = symsum( (1+ez)^(kz-2), kz, 2, Nz); h0z= double(Lz/sz);
Hz = [h0z, h0z]; z = [0,h0z];
himoins1 = h0z; zimoins1=h0z;
for i=3:Nz
    hi =(1+ez)*himoins1; Hz(i) = hi; %H(i) = xi - xi-1
    zi = zimoins1 + hi; z(i) = zi;
    himoins1=hi; zimoins1 = zi;
end


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


%% Résolution dans le temps
Tn = T0;
% Ainv = inv(A);
SSS = zeros(Nx,Ny,Nz,Nt); %Terme source aux temps tn


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
    Tnplus1 = A\b;
    PPP(:, n+1) = Tnplus1;
    Tn = Tnplus1;
    
    %Tmax pour chaque tn
    Tmax(n+1) = max(Tnplus1, [], 'all');
    
    %Région supérieure à température de Curie
    T3D = reshape(Tnplus1, [Nx,Ny,Nz]); Tsup = T3D>Tcurie;
    if nnz(Tsup)    
        xcurie(n+1) = unique(max(X(Tsup), [], 'all'));
        ycurie(n+1) = unique(max(Y(Tsup), [], 'all')); 
        zcurie(n+1) = unique(max(Z(Tsup), [], 'all'));  
    else 
        xcurie(n+1) = 0; ycurie(n+1) =0; zcurie(n+1)=0;
    end
end


%% Calcul de diverses choses
Tint = 400;
[~, index_maxheat] = max(Tmax);

%temps d'écriture
theat = nnz( Tint<Tmax(1:index_maxheat)<Tcurie)*dt;
twrite = nnz( Tmax>=Tcurie)*dt; 
tcool = nnz( Tint<Tmax(index_maxheat:end)<Tcurie)*dt;


end



