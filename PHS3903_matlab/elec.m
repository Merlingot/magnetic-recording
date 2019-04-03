function y = elec(R, t)
% Champ �lectrique �mis par l'antenne

    P = 1e-3*1e2;           %Puissance du laser 
    
    ant = 10e-9;            %Dimension de l'antenne
    k = 100;                %Conductivit� thermique du cobalt
    tau0 = 1.8012e-9;       %tau0
    ni = 4.7332;            %Partie imaginaire de l'indice de r�fraction du cobalt
    Lambda = 780e-9;        %Longueur d'onde
    w = 2*pi*3e8/Lambda;    %Fr�quence
    t0 = 5*tau0;            %D�lai t0
    eps0 = 8.85e-12;        %Permitivit� du vide
    d= [0,0,1e-9];          %distance entre l'antenne et le substrat.

    Einc = P*[1,0,0];
    p = 4*pi*eps0*ant^3*Einc; %Dipole

    R = R + d;    
    r = norm(R);
    ru = R/norm(R);
    Er = ((3*ru*(ru*p')) - p)./(4*pi*eps0*r^3)*exp( -ni*k*r  );
    E = Er*0.5*exp( -(t-t0)^2/(2*tau0^2)  );
    y = norm(E);
end

