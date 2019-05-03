function y = elecM(X,Y,Z, t, Pmoy)
% Calcul du champ électrique émis par l'antenne
% Arguments:
%   X,Y,Z (ndarray): meshgrid des coordonées
%   t (double): temps auquel évaluer le terme source (en seconde)
%   Pmoy (double) : puissance moyenne du laser (en W)
% Returns:
%   y (double): champ électrique 

    %Puissance du laser (crète) % Pcrete*dt=Pmoy*Tauxderépétition
    
    ant = 10e-9;            %Dimension de l'antenne
    k = 100;                %Conductivité thermique du cobalt
    tau0 = 1.8012e-9;       %tau0
    ni = 4.7332;            %Partie imaginaire de l'indice de réfraction du cobalt
    Lambda = 780e-9;        %Longueur d'onde
    w = 2*pi*3e8/Lambda;    %Fréquence
    t0 = 5*tau0;            %Délai tau0
    eps0 = 8.85e-12;        %Permitivité du vide
    mu0 = 4*pi*1e-7;        %Perméabilité du vide
    Z0 = sqrt(mu0/eps0);    %Impédance du vide
    dant = 1e-9;            %Distance entre l'antenne et le substrat (1nm). 
    d = 1e-6;               %Résolution du faiseau
    
    %Champ électrique incident
    Einc = sqrt(4*Z0*Pmoy/pi)/d;
    
    % dipole = p*[1,0,0];
    p = Einc*4*pi*eps0*ant^3;      %Norme du dipole
    px = 1*p; py = 0*p; pz = 0*p;       %Composantes du dipole
       
    ZZ = Z + dant;
    r = ( X.^2 + Y.^2 + ZZ.^2).^(1/2);
    uX = X./r; uY = Y./r; uZ = ZZ./r;
    Er = (4*pi*eps0*r.^3).*exp( -ni*k*r  );
    E = 0.5*exp( -(t-t0)^2/(2*tau0^2));
    
    D = uX.*px + uY.*py + uZ.*pz;
    pX = 3*uX.*D - px; pY = 3*uY.*D - py; pZ = 3*uZ.*D - pz;
    
    normP = ( pX.^2 + pY.^2 + pZ.^2).^(1/2);
    y = normP*E./Er;

end

