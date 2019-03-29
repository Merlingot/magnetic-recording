function [t, Q] = eulerexp(t0, N, dt, T0, f)

    %  Methode d'Euler pour les systemes d'equations 
    %  differentielles de la forme y' = f(t,y)

    %  Entree:
    %  1) f: function handle de la fonction f(t,y) - (string)
    %  2) t0: Temps de depart - (scalaire)
    %  3) T0: Distribution initale de Temp�rature - (Matrice 3D) 
    %  4) dt: Pas de temps - (scalaire)
    %  5) n:  Nombre de pas de temps -(scalaire)

    %  Retour:
    %  1) t est un vecteur contenant les temps ti.
    %  2) Q est un vecteur contenant la solution obtenue.
    %     La rang�e i de Q correspond a la solution de l'equation i (� un temps
    %     ti). 
    %       Q[i] = T(ti) 

    % Initialisation des vecteurs de solutions 
    % 1) t est un vecteur contenant les temps ti.
    t = [1:N].*dt;  
    t(1) = t0;

    % 2) Q est le vecteur contenant les solutions T 
    [dimx, dimy, dimz] = size(T0);
    Q = zeros(N, dimx, dimy, dimz);
    Q(1,:,:,:) = T0;                    %Solution � t0 (connue)

    Tn = T0;                            
    for n = [0:N]
           dT = feval( f, t(n), Tn );   %Calcul de la d�riv�e 
           Tnplus1 = dt*dT + Tn;        %Valeur de T_n+1
           Q(n+1, :, :, :) = Tnplus1;   %Ajout de la solution � Q
           Tn = Tnplus1;                %R�assignation pour prochaine boucle
    end
end



%% Notation

% - indices i,j,k = {0, 1, 2, ..., Nx} : discr�tisation spatialle
% - indice n = {0, 1, 2, ..., N} : discr�tisation temporelle 

% T : 
%   - Matrice 3D (nx, ny, nz)
%   - Distribution de temp�rature � un temps t
%   - T[i,j,k] = Valeur de temp�rature au noeud i,j,k du domaine

% T0 : 
%   - Matrice 3D (nx, ny, nz) 
%   - Distribution initiale de temp�rature ( � t=t0 )
% 
% Q :
%   - Vecteur 4D (N, nx, ny, nz)
%   - Vecteur contenant les T � chaque temps t 
%   - Q[n] = T := distribution de temp�rature au temps tn
% 
