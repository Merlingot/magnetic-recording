% Longueur du domaine
clear all
rho = 8.9e3; cp = 420; K = 100; Ldiff = (K*3e-9/(cp*rho))^(1/2);
Lx = 3*Ldiff; Ly = 3*Ldiff; Lz = 3*Ldiff;
pui=1e5; ee=0.2; Lt = 20e-9; Ny = 20; Nz = 20;

%% Convergence spaciale :
Nt=50; dt = Lt/Nt;
N = [10,11,12,13,14,15,16,17,18,19,20]; %Vecteurs avec le nombre de points en x
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
    Hxn(n) = mean(Hx);
    [sol,t, mem] = feuler_conv(x,y,z,Hx,Hy,Hz,Nt, pui, ee);
    
    %enregistrement du temps et mémoire:
    tN(n) = t; %tN : vecteur avec les temps en fonction de Nx et Nx2
    memN(n) = mem;  %memN : vecteur avec la mémoire occupée par la matrice des coefficients
    
    %calcul de l'erreur
    sol = reshape(sol,[Nx2/2,Ny,Nz,Nt]);
    sol2 = reshape(sol2,[Nx2,Ny,Nz,Nt]);
    sol2 = sol2(1:2:end, :, :, :);
    for l=1:Nt
        errl = sum( abs(sol(:,:,:,l) - sol2(:,:,:,l)), 'all');
        err(l)=errl/sum( abs(sol(:,:,:,l)), 'all');
    end
     %err : vecteur avec les erreurs sur le maillage spatial pout chacun des points de temps
     %errt : vecteur err pour chacun des pas dx (Nx) essayé
    errt(:,n) = err;
end
errmoy = mean(errt); %vecteur avec les erreurs max pour chaque Nx essayé

vect=(0:Nt-1)*dt;
filename='convg_spatiale.mat';
save(filename, 'N', 'errt', 'errmoy','tN', 'memN', 'vect', 'Hx2','Hxn')


%% Convergence temporelle
NN = [10, 15,20,25,30,35,40,45,50]; 

Nx=20; [x,y,z,Hx,Hy,Hz] = pointsnu(Nx,Ny,Nz,ee, 3*Ldiff);

for n=1:length(NN)
    Nt = NN(n); dt = Lt/Nt; vect_dt(n)=dt;
    Nt2 = 2*Nt;
    
    [sol2,t2, mem2] = feuler_conv(x,y,z,Hx,Hy,Hz,Nt2, pui, ee); %2 fois plus de points
    [sol,t, mem] = feuler_conv(x,y,z,Hx,Hy,Hz,Nt, pui, ee);
    
    %enregistrement du temps et mémoire:
    tNN(n) = t;      %tN : vecteur avec les temps en fonction de Nt
    
    %memN : vecteur avec la mémoire occupée par la matrice des coefficients
    %(sert à rien car ne change par avec Nt)
    
    %calcul de l'erreur
    sol2 = sol2(:,1:2:end);
    errxyz = 0;div=0;
    for l=1:Nt
        errxyz = errxyz + abs( sol(:,l) - sol2(:,l));
        div = div + abs(sol(:,l));
    end
    errxyz = errxyz./div; %vecteur avec les erreurs pour chacun des points x,y,z
    errs(:,n) =  errxyz; 
end
errmoy = mean(errs); 

filename='convg_temp.mat';
save(filename, 'NN', 'errs', 'errmoy', 'vect_dt', 'tNN')










