%%% Constantes


P = 1e3;        %Puissance du laser

ant = 10e-9;
k = 100;
tau0 = 1.8012e-9;
sigma = 1.0128e6;
ni = 4.7332;
Lambda = 780e-9;
w = 2*pi*3e8/Lambda;
t0 = 5*tau0;
eps0 = 8.85e-12;
rho = 8.9*1.0000e3;
cp = 420;

Einc = P*[1,0,0];
p = 4*pi*eps0*a^3*Einc;

xi=1;           % coefficient de transfert thermique
c = cp*rho;     % constante c=1/cp*rho
a = 1;          % constante a = k/cp*rho
Tc = 300;       % Température Cobalt (300 K)
Tair = 300;     % Température Air
%%%