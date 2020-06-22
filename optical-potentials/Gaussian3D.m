function [V] = Gaussian3D(X,Y,Z)

Cx = 0;
Cy = 0;
Cz = 0;
Sx = 1;
Sy = 4;
Sz = 8;


V = exp(-(X-Cx).^2/(2*Sx^2)).*exp(-(Y-Cy).^2/(2*Sy^2)).*exp(-(Z-Cz).^2/(2*Sz^2));

end

