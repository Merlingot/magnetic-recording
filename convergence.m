% Longueur du domaine

rho = 8.9e3; cp = 420; K = 100; Ldiff = (K*3e-9/(cp*rho))^(1/2);
Lx = 3*Ldiff; Ly = 3*Ldiff; Lz = 3*Ldiff;
pui=1e5; ee=0.2; Lt = 20e-9; Nt=50; dt = Lt/Nt;

%% Convergence spaciale :
Ny = 20; Nz = 20;
N = [10,11,12,13,14,15,16,17,18,19,20]; %Vecteurs avec le nombre de points en x
ts=[];
%Changement de dx:
for n=1:length(N)
    %2 fois plus de points
    Nx2=2*N(n); 
    [x2,y,z, Hx2, Hy, Hz] = pointsnu(Nx2,Ny,Nz,ee,Ldiff);
    [sol2,t2] = feuler2(x2,y,z,Hx2,Hy,Hz,Nt, pui, ee);
    %normal
    x =  x2(1:2:end); 
    Hx = [x(2)]; 
    for i=1:length(x)-1
        Hx(i+1) = x(i+1) - x(i);
    end
    [sol,t] = feuler2(x,y,z,Hx,Hy,Hz,Nt, pui, ee); ts(n) = t;
    sol = reshape(sol,[Nx2/2,Ny,Nz,Nt]);
    sol2 = reshape(sol2,[Nx2,Ny,Nz,Nt]);
    sol2 = sol2(1:2:end, :, :, :);
    for l=1:Nt
        errl = sum( abs(sol(:,:,:,l) - sol2(:,:,:,l)), 'all');
        err(l)=errl;
    end
    errs(:,n) = err;
end

%%
figure
hold on
vect=0:Nt-1;
for n=1:length(N)
    plot(vect*dt, errs(:,n))
end

figure
plot(N, max(errs), 'o')

%%
filename='convg_spatiale.mat';
save(filename, 'N', 'errs', 'vect','dt', 'Hx2','Hx', 'Hy', 'Hz', ...
     'Ny', 'Nz','pui','ee')
 
 %% Convergence temporelle 
 
 Nt = 
 
 
 
 
 
 
 
 
 
 
 
 
 
