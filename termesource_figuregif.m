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
Nx = 100; Lx = 5*Ldiff;  
Ny = 100; Ly = 5*Ldiff; 
Nz = 100; Lz = 5*Ldiff; 
N = Nx*Ny*Nz;
Lt = 20e-9; Nt = 50;  dt = Lt/Nt;   %Pas de temps

%%% valeurs initiales:
T0 = ones(N,1)*Tc;  % Distribution de température initiale
t0 = 0e-9;          % Temps initial


%% Maillage non uniforme
tic
%i=1:x=0, i=Nx:x=Lx
ex = 0.1; %paramètre de scaling
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
ey = 0.1; %paramètre de scaling
syms ky; sy = symsum( (1+ey)^(ky-2), ky, 2, Ny); h0y = double(Ly/sy);
Hy = [h0y, h0y]; y = [0,h0y]; 
himoins1 = h0y;yimoins1=h0y;
for i=3:Ny
    hi =(1+ey)*himoins1; Hy(i) = hi;
    yi = yimoins1 + hi; y(i) = yi;
    himoins1=hi; yimoins1 = yi;
end

%k=1:z=0, k=Nz:z=Lz
ez = 0.1; %paramètre de scaling
syms kz; sz = symsum( (1+ez)^(kz-2), kz, 2, Nz); h0z= double(Lz/sz);
Hz = [h0z, h0z]; z = [0,h0z]; 
himoins1 = h0z; zimoins1=h0z;
for i=3:Nz
    hi =(1+ez)*himoins1; Hz(i) = hi; %H(i) = xi - xi-1
    zi = zimoins1 + hi; z(i) = zi;
    himoins1=hi; zimoins1 = zi;
end
t3=toc;

[X,Y,Z] = meshgrid(y,x,z);

%%
for n =0:Nt-1
    
    tn = t0 + dt*n;
    vect(n+1) = tn;
    
    % Définir bn (terme source) à chaque temps tn:
    S = c^-1*fSourceM(X,Y,Z,tn, pui);
    SSS(:,:,:,n+1) = S; %pas optimal mais doit etre pareil au terme source
end

%%an = [];
%dim = [0.08 0.08 0.025 0.025];
Nx=30;Ny=30;Nz=30;
%Nx = int16(Nx); Ny=int16(Ny);Nz=int16(Nz);
%Tmax = max(PPP,[], 'all');
%Smax = max(SSS,[], 'all');
XX = X(1:Nx, 1:Ny, 1:Nz)*1e9; YY = Y(1:Nx, 1:Ny, 1:Nz)*1e9; ZZ = Z(1:Nx, 1:Ny, 1:Nz)*1e9; 
xslice = 0; yslice =0; zslice=0;
%%%
fig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'termesource.gif';
intensite=zeros(1,Nt);

for  i = 1:Nt 
    Sr = SSS(1:Nx, 1:Ny, 1:Nz, i);
    slice(XX,YY,ZZ,Sr,[0],[0],[ZZ(end)])
    caxis([0,1e14])
    colorbar
    colormap('parula')%parula
    xlabel('x [nm]')
    ylabel('y [nm]')
    zlabel('z [nm]')
    title('Terme source [W/m^3]')
    shading interp
    pause(0.01) %%Pour voir frame par frame
    
    frame = getframe(fig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
      
    intensite(i)=max(max(max(Sr)));
      
end
  %%
figure
Sr = SSS(1:Nx, 1:Ny, 1:Nz, 24);
    slice(XX,YY,ZZ,Sr,[0],[0],[ZZ(end)])
    caxis([0,1e14])
    colorbar
    colormap('parula')%parula
    xlabel('x [nm]')
    ylabel('y [nm]')
    zlabel('z [nm]')
    title('Terme source [W/m^3]')
    shading interp

fprintf('fin')

%%
plot(vect*1e9, intensite, 'linewidth',1.5)
xlabel('Temps[s]');
ylabel('Terme source [W/m^3]');
grid on