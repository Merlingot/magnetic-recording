ft = load('convg_temp.mat');
fx = load('convg_spatiale.mat');

%% curve fitting
errmoyt =  ft.errmoy; errmoyx =  fx.errmoy;
Hxn = fx.Hxn;
fit1 = fit(ft.vect_dt',errmoyt','power1');
fit2 = fit(fx.Hxn'*1e9,errmoyx','power1');
fit3 = fit(ft.NN', ft.tNN', 'power1');
fit4 = fit(fx.N', fx.tN', 'power1');
% fit5 = fit(fx.N', fx.memN'/1024^3, 'power1');

%% figures 
f1= figure; %erreur temps
hold on
plot(ft.vect_dt, errmoyt,'o');
plot(ft.vect_dt,fit1.a*ft.vect_dt.^fit1.b )
title('Erreur moyenne en fonction du pas de temps')
xlabel('Pas de temps'); ylabel('Erreur moyenne')
legend('données', '$O(h)$', 'Interpreter','latex') ;
saveas(f1,'erreur_temps.png', 'png')
hold off

f2 = figure; %erreur spatial
hold on
plot(Hxn*1e9, errmoyx, 'o');
plot(Hxn*1e9, fit2.a*(Hxn*1e9).^fit2.b)
title('Erreur moyenne en fonction du plus petit pas de discrétisation en x (nm)')
xlabel('Pas de discrétisation en x'); ylabel('Erreur moyenne')
legend('données', '$O(h^3)$', 'Interpreter','latex') 
saveas(f2,'erreur_spatiale.png', 'png')
hold off

f3 = figure; %temps de calcul - temps
hold on
plot(ft.NN, ft.tNN, 'o');
plot(ft.NN, fit3.a*ft.NN.^fit3.b);
title({'Temps de calcul en fonction du nombre';' de points pour le maillage temporel'})
xlabel('Nombre de points'); ylabel('Temps de calcul (s)')
legend('données', '$O(h^3)$', 'Interpreter','latex') 
saveas(f3,'temps_t.png', 'png')
hold off

f4 = figure; %temps de calcul - espace
hold on
plot(fx.N, fx.tN, 'o');
plot(fx.N, fit4.a*fx.N.^fit4.b);
title({'Temps de calcul en fonction du nombre';" de points dans la direction x "})
xlabel('Nombre de points dans la direction x'); ylabel('Temps de calcul (s)')
legend('données', '$O(N_x^3)$', 'Interpreter','latex') 
saveas(f4,'temps_x.png', 'png')
hold off

% f5 = figure;
% hold on
% plot(fx.N, fx.memN/1024^3, 'o');
% plot(fx.N, fit5.a*fx.N.^fit5.b);
% title({'Mémoire utilisée en fonction du nombre'; 'de points dans une dimension (espace)'})
% xlabel('Nombre de points'); ylabel('Mémoire (Gb)')
% legend('données', '$O(N)$', 'Interpreter','latex') 
% saveas(f5,'mem_x.png', 'png')
% hold off

