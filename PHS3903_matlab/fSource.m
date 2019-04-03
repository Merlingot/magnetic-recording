
function y = fSource(R, t)
% Terme source :

    E = elec(R,t);
    sigma = 1.0128e6;       %conductivité du cobalt
    y = sigma*E^2;
    
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