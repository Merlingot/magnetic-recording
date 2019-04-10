clear all
rho = 8.9e3;    % densité du cobalt
cp = 420;       % capacité thermique du cobalt
K = 100;        % conductivité thermique du cobalt
l = (K*3e-9/(cp*rho))^(1/2);

N = 45;e = 0.15;

syms k  ; 
s = symsum( (1+e)^(k-2), k, 2, N);
h0  = double(2*l/s);

H = [0, h0]; X = [0,h0];
himoins1 = h0;ximoins1=h0;
for i=3:N
    hi =(1+e)*himoins1; H(i) = hi;
    xi = ximoins1 + hi; X(i) = xi;
    himoins1=hi; ximoins1 = xi;
end

%%
% plot(H, 'o')
plot(X, exp(-X/l), 'x')



