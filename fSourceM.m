
function y = fSourceM(X,Y,Z, t, pui)
% Calcul du terme source
% Arguments:
%   X,Y,Z (ndarray): meshgrid des coordonées
%   t (double):      temps auquel évaluer le terme source (en seconde)
%   pui (double) :   puissance moyenne du laser (en W)
% Returns:
%   y (ndarray) : terme source

    E = elecM(X,Y,Z,t, pui);
    sigma = 1.0128e6;       %conductivité du cobalt
    y = sigma*E.^2;
    
end
