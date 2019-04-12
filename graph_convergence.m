

%%
F = [f5 f10 f15 f20 f25 f30];
N = [5^3, 10^3, 15^3, 20^3, 25^3, 30^3];
tmat=1:length(F); tres=tmat; tmail=tmat;

for i=1:length(F)
    tmat(i) = F(i).tmat;
    tres(i)= F(i).tres;
    tmail(i) = F(i).tmail;
end

t_tot = tmat+tres+tmail;
plot(N, t_tot);