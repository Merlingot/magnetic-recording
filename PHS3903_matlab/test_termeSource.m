
v = [1,1,1]*10e-9;


%%% test
L = 1e-9; N =50; dx = L/(N-1); 
Lt = 20e-9; Nt = 100; dt = Lt/Nt;

x = linspace(-L,L,N);
z = linspace(0,L,N);
len = length(x);
[X,Y,Z] = meshgrid(x,x,z);
s = zeros(N,N,N);
smaxt = [0:Nt];


for t = [0:Nt]
    for l = [1:len]
        for q = [1:len]
            for j = [1:len]
                xx = X(l,q,j);
                yy = Y(l,q,j);
                zz = Z(l,q,j);
                s(l,q,j) = fSource([xx,yy,zz], t*dt);
            end
        end
    end
    SSS(:,:,:,t+1) = s;
    smaxt(t+1) = max(s,[], 'all');
%     drawnow
%     plot(linspace(0,t*dt, t+1), smaxt(1:t+1));
end

%%

for i = 1:size(SSS,4)
xslice = [];   
yslice = [];
zslice = [0 L/2 L];
slice(X,Y,Z,SSS(:,:,:,i),xslice,yslice,zslice)
colorbar
caxis([0 1e10])
pause
end











