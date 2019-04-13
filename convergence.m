% Longueur du domaine
clear all
rho = 8.9e3; cp = 420; K = 100; Ldiff = (K*3e-9/(cp*rho))^(1/2);
Lx = 3*Ldiff; Ly = 3*Ldiff; Lz = 3*Ldiff;
pui=1e5; ee=0.2; Lt = 20e-9; Ny = 20; Nz = 20;

%% Convergence spaciale :
Nt=50; dt = Lt/Nt;
% N = [10,11,12,13,14,15,16,17,18,19,20]; %Vecteurs avec le nombre de points en x
N = [10,11];
%Changement de dx:
for n=1:length(N)
    %2 fois plus de points
    Nx2=2*N(n); 
    [x2,y,z, Hx2, Hy, Hz] = pointsnu(Nx2,Ny,Nz,ee, 3*Ldiff);
    [sol2,t2, mem2] = feuler_conv(x2,y,z,Hx2,Hy,Hz,Nt, pui, ee);
    
    %normal
    x =  x2(1:2:end);
    Hx = [x(2)];
    for i=1:length(x)-1
        Hx(i+1) = x(i+1) - x(i);
    end
    [sol,t, mem] = feuler_conv(x,y,z,Hx,Hy,Hz,Nt, pui, ee);
    
    %enregistrement du temps et m�moire:
    tN(n) = t; %tN : vecteur avec les temps en fonction de Nx et Nx2
    memN(n) = mem;  %memN : vecteur avec la m�moire occup�e par la matrice des coefficients
    
    %calcul de l'erreur
    sol = reshape(sol,[Nx2/2,Ny,Nz,Nt]);
    sol2 = reshape(sol2,[Nx2,Ny,Nz,Nt]);
    sol2 = sol2(1:2:end, :, :, :);
    for l=1:Nt
        errl = sum( abs(sol(:,:,:,l) - sol2(:,:,:,l)), 'all');
        err(l)=errl;
    end
    errt(:,n) = err;
end
errmax = max(errt);

vect=(0:Nt-1)*dt;
filename='convg_spatiale.mat';
save(filename, 'N', 'errt', 'errmax','tN', 'memN', 'vect', 'Hx2','Hx')


%% Convergence temporelle
NN = [10,15]; 

Nx=20; [x,y,z,Hx,Hy,Hz] = pointsnu(Nx,Ny,Nz,ee, 3*Ldiff);

for n=1:length(NN)
    Nt = NN(n); dt = Lt/Nt; vect_dt(n)=dt;
    Nt2 = 2*Nt;
    
    [sol2,t2, mem2] = feuler_conv(x,y,z,Hx,Hy,Hz,Nt2, pui, ee); %2 fois plus de points
    [sol,t, mem] = feuler_conv(x,y,z,Hx,Hy,Hz,Nt, pui, ee);
    
    %enregistrement du temps et m�moire:
    tN(n) = t;      %tN : vecteur avec les temps en fonction de Nx
    memN(n) = mem;  %memN : vecteur avec la m�moire occup�e par la matrice des coefficients
    
    %calcul de l'erreur
    sol2 = sol2(:,1:2:end);
    for l=1:Nt
        errl = sum( abs(sol(:,l) - sol2(:,l)), 'all');
        err(l)=errl;
    end
%     errt(:,n) = err; %�a marche pas parce que le nombre de points change
%     (Nt)
    
end
errmax = max(errs); 

filename='convg_temp.mat';
save(filename, 'NN', 'errt', 'errmax', 'vect_dt', 'tN', 'memN')









