
%make full grid of quantities
NSites = 80;
[X,Y] = meshgrid(0:NSites-1,0:NSites-1);
n_fn = @(x,y) 0.5 + 0.1*cos(2*pi/5*x);
n_r = n_fn(X,Y);

n_k = fftshift(fft2(ifftshift(n_r)));

[Kxs,Kys] = meshgrid(0:NSites-1,0:NSites-1);
Kxs = 2*pi/NSites*Kxs;
Kys = 2*pi/NSites*Kys;
Kxs = fftshift(Kxs);
Kys = fftshift(Kys);
%take care of potential numerical error problems
Kys(abs(Kys-pi)<1e-15) = pi;
Kxs(abs(Kxs-pi)<1e-15) = pi;
%finish shifting
Kys(Kys>=pi) = Kys(Kys>=pi)-2*pi;
Kxs(Kxs>=pi) = Kxs(Kxs>=pi)-2*pi;

%calculate quantities
Es = -2*(cos(Kxs)+cos(Kys));


Times = linspace(0,10,100);
n_k_evolved = [];
n_r_evolved = [];
for ii = 1:length(Times)
    n_k_evolved = cat(3,n_k_evolved,n_k.*exp(-1i*Es*Times(ii)));
    n_r_evolved = cat(3,n_r_evolved,fftshift(ifft2(ifftshift(n_k_evolved))));
end

