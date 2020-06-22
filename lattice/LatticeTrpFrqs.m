% Estimate lattice trapping frequencies and etc. from beam parameters.

%Fundamental constants
% Pi = 3.14159;
c = 3*10^8;
h = 6.63*10^(-34);
hbar = h/(2*pi);
kb = 1.38*10^(-23);
abohr = 5.2918e-11;

%Lithium
mLi = 9.9883e-27;
gamma = 2 * pi * 6e6;
lambda0 = 671e-9;
f0 = c / lambda0;
w0 = 2 * pi * f0; 

%lattice
lambda = 1064e-9; % * 10^(-9);
f = c / lambda;
w = 2 * pi * f;
k_laser = 2 * pi / lambda;
k_lattice = 2* pi / (lambda / sqrt(2));
Er_lattice = 0.25 * hbar^2 * k_lattice^2 / (2 * mLi); %Correct lattice recoil uses 
% k = pi/a, hence divide by 4 for a = 532*sqrt(2)nm compared to 1064nm.

%measured lattice vals
P = 0.045;
wx0 = 80e-6;
wy0 = wx0;%95e-6;
zx = pi * wx0^2 / lambda;
zy = pi * wy0^2 / lambda;

%potential depth
Int = 2 * P / (pi * wx0 * wy0);
Vdip = 1 * 3 * pi * c^2 / (2 * w0^3) * ( gamma/(w0 - w) + gamma/(w0 + w) ) * Int;
% Lattice depth factor of 16 deeper than single beam for vertical polarization
% factor 4 for horizontal polarization
Vlatt = 16 * abs(Vdip); 

%These expression are accurate for the 'Kai Sun' lattice, with vertical
%polarization.
omegax0_onsite = sqrt(Vlatt * k_laser^2 / mLi); 
omegay0_onsite = sqrt(Vlatt * k_laser^2 / mLi);
omegax0_confinement = sqrt(2 * Vlatt / (mLi*wx0^2));
omegay0_confinement = sqrt(2 * Vlatt / (mLi * wy0^2));
%estimate interaction using SHO ground state wavefunctions
g = 4 * pi * hbar^2 / mLi * 1000 * abohr;
OnSiteInteraction_FromSHOWaveFunction = g *sqrt(mLi * omegax0_onsite / pi / hbar / 2)...
                                          *sqrt(mLi * omegay0_onsite / pi / hbar / 2)...
                                          *sqrt(mLi * 2 * pi * 20e3 / pi / hbar / 2);

disp('Vertical Polarization/Kai Sun Lattice Params:')
fprintf('wx, on-site = (2pi)%0.1f KHz \n', omegax0_onsite / (2*pi) / 1e3)
fprintf('wy, on-site = (2pi)%0.1f KHz \n', omegay0_onsite / (2*pi) / 1e3)
fprintf('wx, confinement = (2pi)%0.1f Hz \n', omegax0_confinement / (2*pi))
fprintf('wy, confinement = (2pi)%0.1f Hz \n', omegay0_confinement / (2*pi))
fprintf('Trap Depth = %0.1f uk \n', Vlatt / kb * 1e6);
fprintf('Trap Depth = %0.1f KHz \n', Vlatt / h / 1e3);
fprintf('Trap Depth = %0.1f Er \n', Vlatt / Er_lattice);
fprintf('Er = %0.2f KHz \n', Er_lattice / h / 1e3);
fprintf('Estimated int from SHO wave function U = %0.1f KHz \n\n',...
         OnSiteInteraction_FromSHOWaveFunction / h / 1e3);

