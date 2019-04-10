
function y = fSourceM(X,Y,Z, t, pui)
% Terme source :

    E = elecM(X,Y,Z,t, pui);
    sigma = 1.0128e6;       %conductivit� du cobalt
    y = sigma*E.^2;
    
end
