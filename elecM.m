function y = elecM(X,Y,Z, t)
% Champ �lectrique �mis par l'antenne

    

    Pui = 1e6;                %Puissance du laser (cr�te)
    % Pcrete*dt=Pmoy*Tauxder�p�tition
    
    ant = 10e-9;            %Dimension de l'antenne
    k = 100;                %Conductivit� thermique du cobalt
    tau0 = 1.8012e-9;       %tau0
    ni = 4.7332;            %Partie imaginaire de l'indice de r�fraction du cobalt
    Lambda = 780e-9;        %Longueur d'onde
    w = 2*pi*3e8/Lambda;    %Fr�quence
    t0 = 5*tau0;            %D�lai t0
    eps0 = 8.85e-12;        %Permitivit� du vide
    d = 1e-9;               %Distance entre l'antenne et le substrat.

    % dipole = p*[1,0,0];
    p = Pui*4*pi*eps0*ant^3;      %Norme du dipole
    px = 1*p; py = 0*p; pz = 0*p;       %Composantes du dipole
    
    
    ZZ = Z + d;
    r = ( X.^2 + Y.^2 + ZZ.^2).^(1/2);
    uX = X./r; uY = Y./r; uZ = ZZ./r;
    Er = (4*pi*eps0*r.^3).*exp( -ni*k*r  );
    E = 0.5*exp( -(t-t0)^2/(2*tau0^2));
    
    D = uX.*px + uY.*py + uZ.*pz;
    pX = 3*uX.*D - px; pY = 3*uY.*D - py; pZ = 3*uZ.*D - pz;
    
    normP = ( pX.^2 + pY.^2 + pZ.^2).^(1/2);
    y = normP*E./Er;
    
    
end

