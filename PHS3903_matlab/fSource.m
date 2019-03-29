
function y = fSource(R, t)
    
    P = 1e-3;               %Puissance du laser 
    
    ant = 10e-9;            %Dimension de l'antenne
    k = 100;                %Conductivité thermique du cobalt
    tau0 = 1.8012e-9;       %tau0
    sigma = 1.0128e6;       %conductivité du cobalt
    ni = 4.7332;            %Partie imaginaire de l'indice de réfraction du cobalt
    Lambda = 780e-9;        %Longueur d'onde
    w = 2*pi*3e8/Lambda;    %Fréquence
    t0 = 5*tau0;            %Délai t0
    eps0 = 8.85e-12;        %Permitivité du vide

    Einc = P*[1,0,0];
    p = 4*pi*eps0*ant^3*Einc; %Dipole

    r = norm(R);
    ru = R/norm(R);
    Er = ((3*ru*(ru*p')) - p)./(4*pi*eps0*r^3)*exp( -ni*k*r  );
    E = Er*sin(w*(t-t0))*exp( -(t-t0)^2/(2*tau0^2)  );
    y = sigma*norm(E)^2;
end


















% function y = fE(R,t)
%     global w t0 tau0 
%     y = Er(R)*sin(w*(t-t0))*exp( -(t-t0)^2/(2*tau0^2)  );
% end 
% 
% function y = fEr(R)
%     global p eps0 ni k 
%     r = norm(R);
%     ru = R/norm(R);
%     y = ((3*ru*(ru*p')) - p)./(4*pi*eps0*r^3)*exp( -ni*k*r  );
% end
% 
% function y = fS(R,t)
%     global sigma 
%     y = sigma*norm(E(R,t))^2;
% end