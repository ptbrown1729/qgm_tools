%Constants
c = 299792458;
hbar = 1.0546e-34;
h = 2*pi*hbar;
mLi = 9.9883e-27;


lambda = 1064e-9;
klight = 2*pi/1064e-9;
P = 2;
Wo = 100e-6;
Int = 2*P/(pi*Wo^2);

omega0 = c*2*pi/671e-9;
omega = c*2*pi/1064e-9;
gamma = 2*pi*6e6;
alpha = ((3*pi*c^2)/(2*omega0^3))*(gamma/(omega0 - omega) + gamma/(omega0 + omega));
Vo = Int*alpha;


klatt = 2*pi/1064e-9*sqrt(2);
Er = hbar^2*klight^2/(2*mLi);
Pot = @(x) 2*Er*(1+cos(klatt*x));

[qvects,Es] = OneDBand(Pot,2*pi/klatt);

for ii = 1:3
    plot(Es(:,ii)/Er)
    hold on;
end