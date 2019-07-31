function gif(X,Y,Z, PPP, q, filename, Nx, Ny, Nz) 

XX = X(1:Nx/q, 1:Ny/q, 1:Nz/q)*1e9; 
YY = Y(1:Nx/q, 1:Ny/q, 1:Nz/q)*1e9; 
ZZ = Z(1:Nx/q, 1:Ny/q, 1:Nz/q)*1e9; 

P = reshape(PPP, [Nx,Ny,Nz,50]);

%%

fig = figure;
axis tight manual % this ensures that getframe() returns a consistent size

for  i = 1:Nt 
    Sr = P(1:Nx/q, 1:Ny/q, 1:Nz/q, i);
    slice(XX,YY,ZZ,Sr,[0],[0],[ZZ(end)])
    caxis([300,1000])
    colorbar
    colormap('parula')%parula
    xlabel('x (nm)')
    ylabel('y (nm)')
    zlabel('z (nm)')
    title('Température (K)')
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
           
end

