f30=load('NU_CS_30_01.mat');
f25=load('NU_CS_25_01.mat');
f20=load('NU_CS_20_01.mat');
f15=load('NU_CS_15_01.mat');
f10=load('NU_CS_10_01.mat');
f5=load('NU_CS_5_01.mat');

%%
F = [f5 f10 f15 f20 f25 f30];
tmat=1:length(F); tres=tmat; tmail=tmat;

for i=1:length(F)
    tmat(i) = F(i).tmat;
    tres(i)= F(i).tres;
    tmail(i) = F(i).tmail;
end

%%