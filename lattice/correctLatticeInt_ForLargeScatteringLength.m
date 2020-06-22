% When the scattering length a_s ~ a_SHO
% where a_SHO is the characteristic length scale of the lattice site SHO
% the interaction is strong enough to modify the shape of the eigenfunction
% (wannier state, in the case of the lattice).
% In this case the naive estimate U = g*int(psi^4)
% no longer gives the correct interaction.
% Isotropic harmonic oscillator: "Two Cold Atoms in a Harmonic Trap" Th. Busch,et al. (1997).
% anisotropic case: PHYSICAL REVIEW A 71, 050701sRd
% https://www.uni-ulm.de/fileadmin/website_uni_ulm/nawi.inst.220/publikationen/PRA71-050701-2005.pdf

hbar = 1.0546e-34;
h = hbar*2*pi;
abohr = 5.2918e-11;
c = 299792458;
mLi = 9.9883e-27;

%SHO Parameters
Omega = (2*pi)*92.6e3; %(2*pi)*48e3; %
aSHO = sqrt(hbar/(mLi*Omega));

%scattering parameters
a0 = 326.8554*abohr; % -4700*abohr;
g = 4*pi*hbar^2/mLi*a0;
EnergyEstimate = (g*sqrt(mLi*Omega/pi/hbar/2)^3)/(hbar*Omega); %dimensionsless
if EnergyEstimate<-1
    Estimate = -0.9;
elseif EnergyEstimate>1
    Estimate = 0.9;
else
    Estimate = EnergyEstimate;
end


%This function is modified slightly from the reference...
%They point out E ~ 3/2+2n*sqrt(2/pi)*(n+1/2)C(n)*a0
%so we actually want to use E -> E-3/2.
E_Eig_Fn = @(Edimensionless) sqrt(2)*gamma(-Edimensionless/2)./gamma(-Edimensionless/2-1/2)-aSHO/a0;
%choose lower bounds appropriate to second potential state (i.e. first
%level that is SHO level and not the bound state
ESolved = lsqnonlin(E_Eig_Fn,Estimate,-1,1);
fprintf('Naive U = %0.1f KHz \n',EnergyEstimate*Omega/(2*pi)/1e3);
fprintf('U = %0.1f KHz \n',ESolved*Omega/(2*pi)/1e3);
fprintf('U/Unaive = %0.3f \n \n',ESolved/EnergyEstimate);