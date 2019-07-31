pmoy = linspace(1e-4,10e-4, 10);
N = 25; Nt=50; 
ee=0.2; Tcurie = 320 + 273.15;

%% Livrable 1 : courbes de puissances et Pmin
for n=(1:length(pmoy))
    [~,~, ~, ~,~ , ~, vect,~,~,~, xcurie, ycurie, zcurie ]= resultat_feuler(Nt, N, N, N, pmoy(n), ee);
    
    dimx(:,n) = xcurie; dimy(:,n) = ycurie; dimz(:,n) = zcurie;
    rx(n) = max(xcurie);  ry(n) = max(ycurie); rz(n) = max(zcurie);   
    
end

%%% Fit puissance et Pmin (R=1nm)
[P,S] = polyfit(pmoy',rx', 1);
pmin = 1e-9/P(1);


%%% Courbes rayon vs. puissance:

ColorSet=varycolor(3);
f1= figure; 
hold on
plot(pmoy*1e3, rx*1e9, '.')
plot(pmoy*1e3, ry*1e9, 'o')
plot(pmoy*1e3, rz*1e9, 'o')
xlabel('Puissance du laser (mW)'); ylabel('Dimension du bit (nm)')
legend('Largeur en x', 'Largeur en y', 'Profondeur en z') ;
% saveas(f1,'rayon_laser.png', 'png')
hold off

%%% Courbes rayon vs. temps en fonction de puissance:

ColorSet=varycolor(10);
f2= figure; 
hold on
for i = 1:length(pmoy)
 plot(vect*1e9, dimy(:,i)*1e9, 'Color',ColorSet(i,:))
end
xlabel('Temps (ns)'); ylabel('Rayon du bit (nm)')
leg = legend('0.1000',    '0.2000',    '0.3000',    '0.4000',    '0.5000',    '0.6000',    '0.7000',    '0.8000',  '0.9000', '1.0000') ;
title(leg,'Puissance (mW)')
% saveas(f2,'rayon_vect.png', 'png')
hold off

%% Livrable 2 : Densité de bits

r1 = 0.5*0.0254; r2=1.7*0.0254; r3=3.5*0.0254;
A = pi*(r2^2 - r1^2); ADD = pi*r3^2;

rbit = [1,5,10,15,20,25]*1e-9*2; 
Abit = pi*rbit.^2;
Nbit = A./Abit;
nbit = Nbit./ADD;

%% Livrable 3 : R=13 nm

p13 = 13e-9/P(1); %Puissance pour laquelle R=13

%%% Calcul des temps d'écriture
[pui, PPP, Tmax, theat, twrite, tcool, vect, X, Y, Z, xcurie, ycurie, zcurie ]= resultat_feuler(Nt, 30, 30, 30, p13, ee);

%%% Tracer le rayon du bit en fonction du temps pour p13
f3= figure; 
hold on
plot(vect*1e9, Tmax );
plot(vect*1e9, vect./vect*Tcurie, '--');
plot(vect*1e9, vect./vect*400, '--');
text(6.5, 375, '$\tau_{heat}$', 'Interpreter','latex');
text(10.5, 375, '$\tau_{cool}$', 'Interpreter','latex');
text(8.5, 375, '$\tau_{write}$', 'Interpreter','latex');
ha = annotation('doublearrow');
ha.Head1Width  = 3;ha.Head2Width  = 3;
ha.Head1Length = 3;ha.Head2Length = 3; 
ha.Parent = gca;           % associate the arrow the the current axes
ha.X = [6.5 7.2];          
ha.Y = [400 400]; 

ha2 = annotation('doublearrow');
ha2.Head1Width  = 3;ha2.Head2Width  = 3;
ha2.Head1Length = 3;ha2.Head2Length = 3; 
ha2.Parent = gca;           % associate the arrow the the current axes
ha2.X = [7.2 10.5];          
ha2.Y = [400 400]; 

ha3 = annotation('doublearrow');
ha3.Head1Width  = 3;ha3.Head2Width  = 3;
ha3.Head1Length = 3;ha3.Head2Length = 3; 
ha3.Parent = gca;           % associate the arrow the the current axes
ha3.X = [10.5 11.5];          
ha3.Y = [400 400]; 
xlabel('Temps (ns)'); ylabel('Température (K)')
legend('Température maximale du cobalt', 'Température de Curie', 'Température de base') ;
% saveas(f3,'Tmax.png', 'png')
hold off

%% Faire un gif:
% filename = 'resultat.gif'; q =1;
% gif(X,Y,Z, PPP, q, filename, N, N, N)
% filename = 'resultat_zoom.gif'; q =2;
% gif(X,Y,Z, PPP, q, filename, N, N, N)


