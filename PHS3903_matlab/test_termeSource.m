
v = [1,1,1]*10e-9;


%%% test
N = 10;
Nt = 1000;
x = [-N:N]*1e-8;
len = length(x);
[X,Y] = meshgrid(x,x);
s = zeros(N,N);
smaxt = [0:Nt];
for t = [0:Nt]
    for l = [1:len]
        for q = [1:len]
            xx = X(l,q);
            yy = Y(l,q);
            s(l,q) = fSource2D([xx,yy], t*1e-10 );
        end
    end
    smaxt(t+1) = max(s, [], 'all');
    drawnow
%     plot([0:t]*1e-10, smaxt(1:t+1) )
    contour(X,Y,s)
end
% %%
% plot([0:Nt]*1e-10,smaxt)









