function [] = slider_plot(x,y,D,freq,range,type)

% [] = slider_plot(x,y,D,freq,range,type)
% x: position in x of the pixel
% y: position in y of the pixel
% D: 3D-Complex-valued FFT data (dim1:x, dim2:y, dim3:fft)
% freq: freq range spanning of the data(in the form -freqmax to +freqmax.
% range: how many points must be averaged?
 
% Plot different plots according to slider location.
S.fh = figure();
%               'units','pixels',...
%               'position',[300 300 300 300],...
%               'menubar','none',...
%               'name','slider_plot',...
%               'numbertitle','off',...
%               'resize','off');             
%S.ax = axes('unit','pix',...
%            'position',[20 80 260 210]);
       
S.sl = uicontrol('style','slide','units','pixels',...
                 'position',[525 30 40 360],...
                 'min',1,'max',length(D),'val',ceil(length(D)/2),...
                 'sliderstep',[1/length(D) .01],...
                 'callback',{@sl_call,S});  
 
             
u1=ceil(length(D)/2);
u2=u1+range;
for i = 1:size(D,1)
    for j=1:size(D,2)
        if strcmp(type,'abs')==1
            E(i,j) = mean(abs(D(i,j,u1:u2)),3);
        elseif strcmp(type,'real')==1
            E(i,j) = mean(real(D(i,j,u1:u2)),3);
        elseif strcmp(type,'imag')==1
            E(i,j) = mean(imag(D(i,j,u1:u2)),3);
        elseif strcmp(type,'angle')==1
            E(i,j) = mean(angle(D(i,j,u1:u2)),3);
        elseif strcmp(type,'logabs')==1
            E(i,j) = mean(log10(abs(D(i,j,u1:u2))),3);
        end  
    end
end
 
imagesc(x,y,(squeeze((E(:,:))))')
xlabel('Position (mm)')
ylabel('Position (mm)')
title(strcat(type,' image',num2str(freq(u1)),' to ','',strcat(num2str(freq(u2))),'THz (ind1= ',num2str(u1),' ind2= ',num2str(u2),')'))
colorbar
%caxis([0 150])
set(gca,'YDir','normal')
set(gcf,'color','w')
%axis equal
 
 
function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
cla
u1=round(get(h,'value'));
u2=u1+range;
for i = 1:size(D,1)
    for j=1:size(D,2)
        if strcmp(type,'abs')==1
            E(i,j) = mean(abs(D(i,j,u1:u2)),3);
        elseif strcmp(type,'real')==1
            E(i,j) = mean(real(D(i,j,u1:u2)),3);
        elseif strcmp(type,'imag')==1
            E(i,j) = mean(imag(D(i,j,u1:u2)),3);
        elseif strcmp(type,'angle')==1
            E(i,j) = mean(angle((D(i,j,u1:u2))));  %mean(atan(imag(D(i,j,u1:u2))./real(D(i,j,u1:u2))),3)
        elseif strcmp(type,'logabs')==1
            E(i,j) = mean(log10(abs(D(i,j,u1:u2))),3);
        end   
    end
end
 
[X,Y]=meshgrid(x,y);
p=pcolor(X,Y,(squeeze((E(:,:))))')
xlabel('Position (mm)')
ylabel('Position (mm)')
title(strcat(type,' image',num2str(freq(u1)),' to ','',strcat(num2str(freq(u2))),'THz (ind1= ',num2str(u1),' ind2= ',num2str(u2),')'))
colorbar
set(gca,'YDir','normal')
colormap jet
%caxis([-pi/2 pi/2])
set(gcf,'color','w')
p = set(p,'EdgeColor','none') 
 
end
end