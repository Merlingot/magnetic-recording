%% Constantes &  Param�tres
clear all
rho = 8.9e3;    % densit� du cobalt
cp = 420;       % capacit� thermique du cobalt
K = 100;        % conductivit� thermique du cobalt
E = 10;         % coefficient de transfert thermique air-cobalt (sans rotation du disque)
c = cp*rho;     % constante c=cp*rho
a = c/K;        % constante alpha = cp*rho/K

Ldiff = (K*3e-9/(cp*rho))^(1/2);    %Longueur de diffusion

Tc = 300;                   % Temp�rature Cobalt (300 K)
Tair = 300;                 % Temp�rature Air
TCurie = 320 + 273.15;      % Temp�rature de Curie du Cobalt (320 C)

%%% PARAM�TRES %%%%%%%%%%
% Temps
t0 = 0e-9; % Temps initial
Lt = 20e-9; Nt = 50;  dt = Lt/Nt;   %Pas de temps
% Longueur du domaine
Lx = 3*Ldiff; Ly = 3*Ldiff; Lz = 3*Ldiff;
% Puissance laser
pui = 1e5;
% Stretching 
ee=0.1;


tic
for nbre=5:5:30
    %% D�finition du nombre de points & Maillage
    Nx = nbre; Ny=nbre; Nz=nbre; %Nombre de points dans le maillage
    N = Nx*Ny*Nz;
    
    filename = sprintf ( 'NU_CS_%d_%d_%d.mat', nbre, round(ee) );
    
    %% D�finition du Maillage non uniforme
    tic
    
    %i=1:x=0, i=Nx:x=Lx
    ex = ee; %param�tre de scaling
    syms kx; sx = symsum( (1+ex)^(kx-2), kx, 2, Nx);
    h0x = double(Lx/sx);
    Hx = [h0x, h0x]; x = [0,h0x]; %Hx:pas en x, x:coordon�es en x
    himoins1 = h0x;ximoins1=h0x;
    for i=3:Nx
        hi =(1+ex)*himoins1; Hx(i) = hi;
        xi = ximoins1 + hi; x(i) = xi;
        himoins1=hi; ximoins1 = xi;
    end
    
    %j=1:y=0, j=Ny:y=Ly
    ey = ee; %param�tre de scaling
    syms ky; sy = symsum( (1+ey)^(ky-2), ky, 2, Ny); h0y = double(Ly/sy);
    Hy = [h0y, h0y]; y = [0,h0y];
    himoins1 = h0y;yimoins1=h0y;
    for i=3:Ny
        hi =(1+ey)*himoins1; Hy(i) = hi;
        yi = yimoins1 + hi; y(i) = yi;
        himoins1=hi; yimoins1 = yi;
    end
    
    %k=1:z=0, k=Nz:z=Lz
    ez = ee; %param�tre de scaling
    syms kz; sz = symsum( (1+ez)^(kz-2), kz, 2, Nz); h0z= double(Lz/sz);
    Hz = [h0z, h0z]; z = [0,h0z];
    himoins1 = h0z; zimoins1=h0z;
    for i=3:Nz
        hi =(1+ez)*himoins1; Hz(i) = hi; %H(i) = xi - xi-1
        zi = zimoins1 + hi; z(i) = zi;
        himoins1=hi; zimoins1 = zi;
    end
    tmail=toc;
    
    [X,Y,Z] = meshgrid(x,y,z);
    
    %% Initialisation des matrices et vecteurs
    
    %statiques:
    T0 = ones(N,1)*Tc;    %vecteur Nx1 : Distribution de temp�rature initiale
    M = sparse(N,N);      %sparse NxN identity matrix
    A = sparse(N,N);      %sparse NxN zeros matrix
    b0 = zeros(N,1);      %vecteur Nx1 : partie statique du vecteur bn
    I = speye(N,N);         %matrice qui multiplie le vecteur bn (�l�ments nuls pour x,y,z=Nx,Ny,Nz)
    %dynamiques:
    vect = zeros(0,Nt);   %vecteur des temps tn
    PPP = zeros(N, Nt);   %vecteur avec la solution (temp�rature) � tous les temps tn
    
    
    % D�finition de la matrice des coefficients
    tic
    for k = 1:Nz
        for j=1:Ny
            for i=1:Nx
                
                % remplir la ligne l :
                l = (k-1)*Nx*Ny + (j-1)*Nx + i;
                
                % 1. Conditions de Diriclet -------------
                if (i==Nx)
                    % noeud � la surface interne  x=Lx
                    A(l,l)=1; b0(l) = Tc; I(l,l)=0;
                    
                elseif (j==Ny)
                    % noeud � la surface interne  y=Ly
                    A(l,l)=1; b0(l) = Tc; I(l,l)=0;
                    
                elseif (k==Nz)
                    % noeud � la surface interne  z=Lz
                    A(l,l)=1; b0(l) = Tc; I(l,l)=0;
                % ----------------------------------------
                    
                % ----------------------------------------
                else
                    % Contribution du point i,j,k
                    
                    hx = 1/(Hx(i+1)*Hx(i)); hy = 1/(Hy(j+1)*Hy(j)); hz = 1/(Hz(k+1)*Hz(k));
                    A(l,l) = 1 + 2*dt*a^-1*( hx + hy + hz);
                    M(l,l) = 1;
                    % La contribution du terme source sera initialis�e dans la r�solution
                    % temporelle

                    % ----------------------------------------
                    if (i==1)  %r�flexion en x=0
                        
                        hx1 = Hx(i); A(l,l) = A(l,l) - dt*a^-1*hx1^-2*(4/3); 
                        % Contribution supl�mentaire du point i+1,j,k
                        c = index(i+1,j,k,Nx,Ny); A(l,c) = -1*dt*a^-1*hx1^-2*(1 -1/3);
                    else
                        % Contribution du point i-1,j,k
                        c = index(i-1,j,k,Nx,Ny); hx=2/(Hx(i)*(Hx(i+1)+Hx(i)));
                        A(l,c) = -1*dt*a^-1*hx;
                        % Contribution du point i+1,j,k
                        c = index(i+1,j,k,Nx,Ny); hx=2/(Hx(i+1)*(Hx(i+1)+Hx(i)));
                        A(l,c) = -1*dt*a^-1*hx;
                    end
                    
                    if (j==1) %r�flexion en y=0
                        
                        hy1 = Hy(j); A(l,l) = A(l,l) - dt*a^-1*hy1^-2*(4/3);
                        % Contribution supl�mentaire du point i,j+1,k
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
    tmat= toc;
    
    
    %% R�solution dans le temps
    bn = zeros([N,1]); Tn = T0;
    Ainv = inv(A);
    SSS = zeros(Nx,Ny,Nz,Nt); %Terme source aux temps tn
    
    tic
    for n =0:Nt
        % Temps
        tn = t0 + dt*n;
        vect(n+1) = tn; 
        
        % D�finir bn (terme source) � chaque temps tn:
        S = K^-1*fSourceM(X,Y,Z,tn, pui);
        SSS(:,:,:,n+1) = S;
        vecS = I*reshape(S,[N,1]);
        bn = b0 + dt*a^-1*vecS;
        
        % R�solution
        b = M*Tn + bn;
        Tnplus1 = Ainv*b;
        PPP(:, n+1) = Tnplus1;
        Tn = Tnplus1;
    end
    tres=toc;
    
    %% Sauver les donn�es du workspace :
    save(filename, 'nbre', 'pui', 'vect', 'PPP', 'SSS', 'X', 'Y', 'Z',...
    'Hx', 'Hy', 'Hz', 'tmail', 'tmat', 'tres', 't_tot', 'Nx', 'Ny', 'Nz', 'x', 'y', 'z')
end

t_tot=toc;
