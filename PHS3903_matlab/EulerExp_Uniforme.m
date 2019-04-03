clear all
%% Constantes
rho = 8.9*1.0000e3;     %densité du cobalt
cp = 420;               %capacité thermique du cobalt

k = 100;        %Conductivité thermique du cobalt
xi=10;          % coefficient de transfert thermique Air-Cobalt (sans rotation du disque)
c = cp*rho;     % constante c=cp*rho
a = c/k;        % constante alpha = cp*rho/k

Tc = 300;       % Température Cobalt (300 K)
Tair = 300;     % Température Air

%% Paramètres 
global Nx Ny
Nx = 10; Lx = 10e-9; hx = Lx/(Nx-1);
Ny = 10; Ly = 10e-9; hy = Ly/(Ny-1);
Nz = 10; Lz = 10e-9; hz = Lz/(Nz-1);
N = Nx*Ny*Nz;
Nt = 100;  dt = (a*hx^2); Lt = Nt*dt;  %Pas de temps défini à l'aide de la condition de stabilité
%%%

T0 = ones(Nx,Ny,Nz)*Tc;   %Distribution de température initiale
S0 = zeros(Nx,Ny,Nz);     %Terme source
Tnplus1 = ones(Nx,Ny,Nz);

%% Résolution dans le temps
Sn = S0; Tn = T0;

for n = 0:Nt
    % définir Sn (terme source) à chaque temps tn:
    for k = 1:Nz
        for j=1:Ny
            for i=1:Nx
                Sn(i,j,k) =  fSource([i*hx,j*hy,k*hz], dt*n + 7e-9);
            end
        end
    end
    
    % noeud qui est strictement à l'intérieur de la cellule de simulation
    for k = 2:Nz-1
        for j=2:Ny-1
            for i=2:Nx-1
                Tnplus1(i,j,k) = (dt/a)*( (Tn(i-1,j,k)+Tn(i+1,j,k))*hx^2 + ...
                    (Tn(i,j-1,k)+Tn(i,j+1,k))*hy^2 + (Tn(i,j,k-1)+Tn(i,j,k+1))*hz^2 ) ...
                    + Tn(i,j,k)*(1 - 2*dt*a^-1*(hx^-1 + hy^-1 + hz^-1) ) + Sn(i,j,k)*dt*a^-1*k^-1 ;  
            end
        end
    end
    
    
    for k=1:Nz
        for j=1:Ny
            for i=1:Nx
                
                %  Conditions de Neumann --------------
                if(i==1)% noeud sur le plan de symétrie x=0
                    Tnplus1(1,j,k) = (4*Tnplus1(2,j,k) - Tnplus1(3,j,k))/3;
                elseif(j==1) % noeud sur le plan de symétrie y=0
                    Tnplus1(i,1,k) = (4*Tnplus1(i,2,k) - Tnplus1(i,3,k))/3;
                    % Conditions de Diriclet --------------
                elseif (i==Nx) % noeud à la surface interne  x=Lx
                    Tnplus1(i,j,k) = Tc;
                elseif (j==Ny) % noeud à la surface interne  y=Ly
                    Tnplus1(i,j,k) = Tc;
                elseif (k==Nz) % noeud à la surface interne  z=Lz
                    Tnplus1(i,j,k) = Tc;
                    %Condition de Convection  --------------
                elseif (k==1) % noeud à la surface externe  z=0
                    Tnplus1(i,j,k) = (-2*xi*Tair*hz - 4*k*Tnplus1(i,j,2) + k*Tnplus1(i,j,3) )/(-2*xi*hz - 3*k);
                    
                end
            end
        end
    end
    
    % Tn+1n:
    Tn = Tnplus1;

end

%% Visualisation
Tr = reshape(Tnplus1,[Nx,Ny,Nz]);
x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
z = linspace(0,Lz,Nz);
[X,Y,Z] = meshgrid(x,y,z);
xslice = [];   
yslice = [];
zslice = [0, Lz/2, Lz];
slice(X,Y,Z,Tr,xslice,yslice,zslice)
colorbar
xlabel('x')
ylabel('y')
zlabel('z')
