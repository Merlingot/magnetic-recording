run('convergence.m');
ft = load('convg_temp.mat');
fx = load('convg_spatiale.mat');

%% curve fitting
errmoyt =  ft.errmoy; errmoyx =  fx.errmoy;
Hxn = fx.Hxmin; %pas d'espace le plus petit
fit1 = fit(ft.vect_dt',errmoyt','power1'); %erreur temps
fit2 = fit(Hxn'*1e9,errmoyx','power1'); %erreur espace
fit3 = fit(ft.NN', ft.tNN', 'power1');     %temps de calcul temps
fit4 = fit(fx.N', fx.tN', 'power1');       %temps de calcul espace-Nx
fit6 = fit(fx.N'*20*20, fx.tN', 'power1'); %temps de calcul espace-Ntot
% fit5 = fit(fx.N', fx.memN'/1024^3, 'power1'); %Mémoire
fit7 = fit(fx.N'*20*20, fx.tmN', 'power1'); %temps matrice des coeff espace-Ntot
%% figures
f1= figure; %erreur temps
hold on
loglog(ft.vect_dt, errmoyt,'o');
loglog(ft.vect_dt,fit1.a*ft.vect_dt.^fit1.b )
title('Erreur relative moyenne en fonction du pas de temps')
xlabel('Pas de temps (s)'); ylabel('Erreur relative moyenne')
legend('données', '$\sim O(h)$', 'Interpreter','latex') ;
saveas(f1,'erreur_temps.png', 'png')
hold off
%%
f2 = figure; %erreur spatiale
hold on
plot(Hxn*1e9, errmoyx, 'o');
plot(Hxn*1e9, fit2.a*(Hxn*1e9).^fit2.b)
% title({'Erreur relative moyenne en fonction du'; 'pas de discrétisation minimal dans la direction x'})
xlabel('Pas de discrétisation minimal en x (nm)'); ylabel('Erreur relative moyenne')
legend('données', '$\sim O(h^{0.5})$', 'Interpreter','latex')
saveas(f2,'erreur_spatiale.png', 'png')
hold off
%%
f3 = figure; %temps de calcul - temps
hold on
loglog(ft.NN, ft.tNN, 'o');
loglog(ft.NN, fit3.a*ft.NN.^fit3.b);
% title({'Temps de calcul en fonction du nombre';' de noeuds dans le maillage temporel'})
xlabel('Nombre de noeuds (N)'); ylabel('Temps de calcul (s)')
legend('données', '$\sim O(N)$', 'Interpreter','latex')
saveas(f3,'temps_t.png', 'png')
hold off
%%
f4 = figure; %temps de calcul - espace
hold on
loglog(fx.N, fx.tN, 'o');
loglog(fx.N, fit4.a*fx.N.^fit4.b);
% title({'Temps de calcul en fonction du nombre';" de noeuds dans la direction x "})
xlabel('Nombre de noeuds dans la direction x'); ylabel('Temps de calcul (s)')
legend('données', '$\sim O(N^2)$', 'Interpreter','latex')
saveas(f4,'temps_x.png', 'png')
hold off
%%
f6 = figure; %temps de calcul - espace avec nombre de points total
hold on
loglog(fx.N*20*20, fx.tN, 'o');
loglog(fx.N*20*20, fit6.a*(fx.N*20*20).^fit6.b);
% title({'Temps de calcul en fonction du nombre';' de noeuds dans le maillage spatial '})
xlabel('Nombre de noeuds dans le maillage spatial (N)'); ylabel('Temps de calcul (s)')
legend('données', '$\sim O(N^2)$', 'Interpreter','latex')
saveas(f6,'temps_xtot.png', 'png')
hold off

%%
f7 = figure; %temps coeff - espace avec nombre de points total
hold on
loglog(fx.N*20*20, fx.tmN, 'o');
loglog(fx.N*20*20, fit7.a*(fx.N*20*20).^fit7.b);
title({'Temps de calcul en fonction du nombre';' de noeuds dans le maillage spatial '})
xlabel('Nombre de noeuds dans le maillage spatial (N)'); ylabel('Temps de calcul (s)')
legend('données', '$\sim O(N^2)$', 'Interpreter','latex')
saveas(f7,'tempscoeff_xtot.png', 'png')
hold off
%% 
% f5 = figure; %mémoire
% hold on
% plot(fx.N, fx.memN/1024^3, 'o');
% plot(fx.N, fit5.a*fx.N.^fit5.b);
% title({'Mémoire utilisée en fonction du nombre'; 'de points dans une dimension (espace)'})
% xlabel('Nombre de points'); ylabel('Mémoire (Gb)')
% legend('données', '$O(N)$', 'Interpreter','latex')
% saveas(f5,'mem_x.png', 'png')
% hold off

