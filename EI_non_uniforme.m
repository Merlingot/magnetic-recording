%% Constantes &  Paramètres
clear all
rho = 8.9e3;    % densité du cobalt
cp = 420;       % capacité thermique du cobalt
K = 100;        % conductivité thermique du cobalt
E = 10;         % coefficient de transfert thermique air-cobalt (sans rotation du disque)
c = cp*rho;     % constante c=cp*rho
a = c/K;        % constante alpha = cp*rho/K

Ldiff = (K*3e-9/(cp*rho))^(1/2);    %Longueur de diffusion

Tc = 300;                   % Température Cobalt (300 K)
Tair = 300;                 % Température Air
TCurie = 320 + 273.15;      % Température de Curie du Cobalt (320 C)

pui = 1e5; %Puissance laser

%%% paramètres :
Nx = 15; Lx = Ldiff;  
Ny = 10; Ly = Ldiff; 
Nz = 10; Lz = Ldiff; 
N = Nx*Ny*Nz;
Lt = 20e-9; Nt = 50;  dt = Lt/Nt;   %Pas de temps

%%% valeurs initiales:
T0 = ones(N,1)*Tc;  % Distribution de température initiale
t0 = 0e-9;          % Temps initial


%% Maillage non uniforme
tic
%i=1:x=0, i=Nx:x=Lx
ex = 0.2; %paramètre de scaling
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
ey = 0.2; %paramètre de scaling
syms ky; sy = symsum( (1+ey)^(ky-2), ky, 2, Ny); h0y = double(Ly/sy);
Hy = [h0y, h0y]; y = [0,h0y]; 
himoins1 = h0y;yimoins1=h0y;
for i=3:Ny
    hi =(1+ey)*himoins1; Hy(i) = hi;
    yi = yimoins1 + hi; y(i) = yi;
    himoins1=hi; yimoins1 = yi;
end

%k=1:z=0, k=Nz:z=Lz
ez = 0.2; %paramètre de scaling
syms kz; sz = symsum( (1+ez)^(kz-2), kz, 2, Nz); h0z= double(Lz/sz);
Hz = [h0z, h0z]; z = [0,h0z]; 
himoins1 = h0z; zimoins1=h0z;
for i=3:Nz
    hi =(1+ez)*himoins1; Hz(i) = hi; %H(i) = xi - xi-1
    zi = zimoins1 + hi; z(i) = zi;
    himoins1=hi; zimoins1 = zi;
end
t3=toc;

[X,Y,Z] = meshgrid(y,x,z); %Attention! il faut inverser x et y à cause de meshgrid

%% Initialisation des matrices M et A

M = sparse(N,N);      %sparse NxN identity matrix - Ne change pas dans le temps
A = sparse(N,N);      %sparse NxN zeros matrix - Ne change pas dans le temps
b0 = zeros(N,1);      %vecteur 1xN 
vect = zeros(0,Nt);   %vecteur des temps tn
PPP = zeros(N, Nt);   %vecteur avec les solutions à tous les temps tn
I = speye(N,N);         %matrice qui multiplie le vecteur b (éléments nuls pour x,y,z=Nx,Ny,Nz)

% Définition de la matrice des coefficients 
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
t1 = toc;


%% Résolution dans le temps
bn = zeros([N,1]); Tn = T0;
Ainv = inv(A);
SSS = zeros(Nx,Ny,Nz,Nt);

tic
for n =0:Nt-1
    
    tn = t0 + dt*n;
    vect(n+1) = tn;
    
    % Définir bn (terme source) à chaque temps tn:
    S = K^-1*fSourceM(X,Y,Z,tn, pui);
    SSS(:,:,:,n+1) = S; %pas optimal mais doit etre pareil au terme source
    vecS = I*reshape(S,[N,1]);
    bn = b0 + dt*a^-1*vecS;
    
    %Résolution
    b = M*Tn + bn; 
    Tnplus1 = Ainv*b;
    PPP(:, n+1) = Tnplus1; 
    Tn = Tnplus1;
end
t2=toc;

%% Visualisation
%%%POUR L'AFFICHAGE
an = [];
dim = [0.08 0.08 0.025 0.025];
q=int16(1);
Nx = int16(Nx); Ny=int16(Ny);Nz=int16(Nz);
Tmax = max(PPP,[], 'all');
Smax = max(SSS,[], 'all');
XX = X(1:Nx/q, 1:Ny/q, 1:Nz/q)*1e9; YY = Y(1:Nx/q, 1:Ny/q, 1:Nz/q)*1e9; ZZ = Z(1:Nx/q, 1:Ny/q, 1:Nz/q)*1e9; 
xslice = 0; yslice =0; zslice=0;
%%%

for  i = 1:Nt+1 
    Trr = reshape(PPP(:,i),size(X));
    Sr = SSS(1:Nx/q, 1:Ny/q, 1:Nz/q, i);
    Tr = Trr(1:Nx/q, 1:Ny/q, 1:Nz/q);
    subplot(1,3,1)
    slice(XX,YY,ZZ,Tr,xslice,yslice,zslice)
%     caxis([300 Tmax])
    colorbar
    colormap('hot')
    xlabel('x (nm)')
    ylabel('y (nm)')
    zlabel('z (nm)')
    title('Température (K)')
    ti = vect(i);
    str = ['t = ' num2str(ti) ' s'];
    delete(an); 
    an = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    subplot(1,3,2)
    Tbool = double(Tr>TCurie);
    slice(XX,YY,ZZ,Tbool,xslice,yslice,zslice)
    colorbar
    colormap('hot')
    xlabel('x (nm)')
    ylabel('y (nm)')
    zlabel('z (nm)')
    title('Régions de température supérieure à la température de Curie')
    caxis([0 1])
    
    subplot(1,3,3)
    slice(XX,YY,ZZ,Sr,[0],[0],[0])
%     caxis([0,Smax])
    colorbar
    colormap('hot')
    xlabel('x (nm)')
    ylabel('y (nm)')
    zlabel('z (nm)')
    title('Terme source')
    drawnow %%Pour voir frame par frame
end

fprintf('fin')
