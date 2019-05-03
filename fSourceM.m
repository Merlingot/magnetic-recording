
function y = fSourceM(X,Y,Z, t, pui)
% Calcul du terme source
% Arguments:
%   X,Y,Z (ndarray): meshgrid des coordon�es
%   t (double):      temps auquel �valuer le terme source (en seconde)
%   pui (double) :   puissance moyenne du laser (en W)
% Returns:
%   y (ndarray) : terme source

    E = elecM(X,Y,Z,t, pui);
    sigma = 1.0128e6;       %conductivit� du cobalt
    y = sigma*E.^2;
    
end
