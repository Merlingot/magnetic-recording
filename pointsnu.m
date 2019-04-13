
function [x, y, z, Hx, Hy, Hz] = pointsnu(Nx, Ny, Nz, ee, L)

Lx = L; Ly = L; Lz= L;

%i=1:x=0, i=Nx:x=Lx
    ex = ee; %paramètre de scaling
    syms kx; sx = symsum( (1+ex)^(kx-2), kx, 2, Nx);
    h0x = double(Lx/sx);
    Hx = [h0x, h0x]; x = [0,h0x]; %Hx:pas en x, x:coordonées en x
    himoins1 = h0x;ximoins1=h0x;
    for i=3:Nx
        hi =(1+ex)*himoins1; Hx(i) = hi;
        xi = ximoins1 + hi; x(i) = xi;
        himoins1=hi; ximoins1 = xi;
    end

    %j=1:y=0, j=Ny:y=Ly
    ey = ee; %paramètre de scaling
    syms ky; sy = symsum( (1+ey)^(ky-2), ky, 2, Ny); h0y = double(Ly/sy);
    Hy = [h0y, h0y]; y = [0,h0y];
    himoins1 = h0y;yimoins1=h0y;
    for i=3:Ny
        hi =(1+ey)*himoins1; Hy(i) = hi;
        yi = yimoins1 + hi; y(i) = yi;
        himoins1=hi; yimoins1 = yi;
    end

    %k=1:z=0, k=Nz:z=Lz
    ez = ee; %paramètre de scaling
    syms kz; sz = symsum( (1+ez)^(kz-2), kz, 2, Nz); h0z= double(Lz/sz);
    Hz = [h0z, h0z]; z = [0,h0z];
    himoins1 = h0z; zimoins1=h0z;
    for i=3:Nz
        hi =(1+ez)*himoins1; Hz(i) = hi; %H(i) = xi - xi-1
        zi = zimoins1 + hi; z(i) = zi;
        himoins1=hi; zimoins1 = zi;
    end
    
end
