% Test expressions given in lattice.pdf for fourier components of potential
% by comparing the fourier expression with the direct expression.

E = 0.326;
chi = 1;
alpha = 0.6;
r = 0.6;
theta = 90 * pi/180;
k = 1;
depth = abs(chi) * E^2 * 4 * ( (1+r^2) * (cos(alpha)^2 + cos(theta) * sin(alpha)^2) + 2*r);

% Two options: take E = E_o*e^{i(kx-wt)} and
% either regarded actual field as Re{E}, or as E + E*. These
% conventions result in a factor of two difference in depth for given E_o.
% In the first case, <Re{E}^2>       = 0.5 * E E*
% In the second    , <(E+E*)(E+E*)>  = 2.0 * E E* 
% In any case, the `real' quantity is the intensity. So in some sense this
% is defining the real electric field.
E1 = @(x, y) E * [sin(alpha) * cos(theta/2), -sin(alpha) * sin(theta/2), cos(alpha)] * ...
     exp(1i * k * (x * sin(theta/2) + y * cos(theta/2)) );
E2 = @(x, y) E * -[sin(alpha) * cos(theta/2), sin(alpha) * sin(theta/2), cos(alpha)] * ...
     exp(1i * k * (x * sin(theta/2) - y * cos(theta/2)) );
E3 = @(x,y) E * -[sin(alpha) * cos(theta/2), sin(alpha) * sin(theta/2), cos(alpha)] * ...
     r * exp(-1i * k * (x * sin(theta/2) - y * cos(theta/2)) );
E4 = @(x, y) E * [sin(alpha) * cos(theta/2), -sin(alpha) * sin(theta/2), cos(alpha)] * ...
     r * exp(-1i * k * (x * sin(theta/2) + y * cos(theta/2)) );
E_tot = @(x,y) E1(x,y) + E2(x,y) + E3(x,y) + E4(x,y);
% This defines how I get the intensity. Saying that for a single beam the
% real electric field would be \sqrt{2} * E_o * cos(kx - omega*t), which is
% a little odd, but doesn't matter.
v_from_e = @(x,y) -chi * sum(E_tot(x,y) .* conj(E_tot(x,y)));

% fourier components
v00 = -chi * E^2 * 2 * (1+r^2);
v10 = chi * E^2 * 2 * r * (cos(alpha)^2 + cos(theta) * sin(alpha)^2);
v01 = chi * E^2 * (1 + r^2) * (cos(alpha)^2 + cos(theta) * sin(alpha)^2); 
v11 = chi * E^2 * (-r);
v_fourier = @(x,y) v00 + ...
                   2 * v10 * cos(2 * k * sin(theta/2) * x ) + ...
                   2 * v01 * cos(2 * k * cos(theta/2) * y ) + ...
                   2 * v11 * cos(2 * k * sin(theta/2) * x + 2 * k * cos(theta/2) * y) + ...
                   2 * v11 * cos(2 * k * sin(theta/2) * x - 2 * k * cos(theta/2) * y);

% plot results
pts = linspace(-5, 5, 200);
[xx, yy] = meshgrid(pts, pts);

v_final_test = zeros(size(xx));
v_efield = zeros(size(xx));
v_fourier_test = zeros(size(xx));
for ii = 1 :numel(xx)
%     v_final_test(ii) = v(xx(ii), yy(ii));
    v_fourier_test(ii) = v_fourier(xx(ii), yy(ii));
    v_efield(ii) = v_from_e(xx(ii), yy(ii));
end
depth_v_efield = max(v_efield(:)) - min(v_efield(:));
depth_v_fourier = max(v_fourier_test(:)) - min(v_fourier_test(:));

figure;
subplot(1, 2, 1)
imagesc(v_efield);
axis equal;
axis image;

subplot(1, 2, 2)
imagesc(v_fourier_test);
axis equal;
axis image;

fprintf('%0.3e is maximum different between direct calculation and fourier components\n',...
        max(abs(v_efield(:) - v_fourier_test(:))));
fprintf('calculated depth %0.3f\n', depth);
fprintf('V from fourier component depth %0.3f\n', depth_v_fourier);
fprintf('V from efield depth %0.3f\n', depth_v_efield);