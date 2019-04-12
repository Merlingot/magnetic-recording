f = load('convg_spa_NU_30_0.1.mat');
PPP = f.PPP; SSS=f.SSS; 
X=f.X; Y=f.Y; Z=f.Z;
Nx=f.Nx; Ny=f.Ny; Nz=f.Nz;
vect= f.vect;

Nt = length(vect);

TCurie = 320 + 273.15;
%%
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

for  i = 1:Nt 
    Trr = reshape(PPP(:,i),[Nx,Ny,Nz]);
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