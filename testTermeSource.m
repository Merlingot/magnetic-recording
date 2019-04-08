
%%% test
L = 1e-9; N =15; dx = L/(N-1);
Lt = 20e-9; Nt = 100; dt = Lt/Nt;

x = linspace(-L,L,N);
z = linspace(0,L,N);
len = length(x);
[X,Y,Z] = meshgrid(x,x,z);

xslice = [];
yslice = [];
zslice = [0,L];

for i = 0:Nt
    
    title('matriciel');
    y = fSourceM(X,Y,Z,i*dt);
    slice(X,Y,Z,y,xslice,yslice,zslice)
    colorbar
    drawnow
    
end