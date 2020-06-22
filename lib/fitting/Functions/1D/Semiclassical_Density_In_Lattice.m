function [V,PNames,ArgStr] = Semiclassical_Density_In_Lattice(P,X)
PNames = {'wtrap','t','kT','mu','fudge'};
ArgStr = '1';
% P = [wtrap,t,kT,mu]
% returns radial profile of #/site with pixel (1pix=8.66e-7m) spacing. Profile is
% calculated using semiclassical approximation for non-interacting Fermi
% gas in lattice+Gaussian trap. Assumptions are LDA, tight binding
% (sinusoid bands), all population confined to ground band of lattice and
% ground and first excited state in axial confinement.
% [wtrap] = harmonic confinement at center of Gaussian trap (Hz)
% [t] = tunneling (Er)
% [kT] = temperature (Er)
% [mu] = chemical potential (Er)

if isempty(P)||isempty(X)
    V = 0;
    return
end

M = 6*1.66e-27; hbar = 1.055e-34; lambda = 1064e-9; 
a = lambda/(2^0.5); % Kai Sun lattice
Er = 0.5*(pi*hbar)^2/(M*a^2);
pixelSize = 8.66e-7;

wtrap = P(1)*2*pi; t = P(2); kT = P(3); mu = P(4); fudge = P(5);

trapradius = 90e-6; % 1/e^2 beam radius for trap
trapdepth = 0.25*M*wtrap^2*trapradius^2/Er;
trapParam = 2*pixelSize^2/(trapradius^2);
wa = 1.4; %axial trapping frequency

%hbar set to 1; all units in Er = 14.7kHz

qpointnumber = 60;
[qX,qY] = meshgrid(0.5:qpointnumber-0.5,0.5:qpointnumber-0.5);
qX = qX*1/qpointnumber; qY = qY*1/qpointnumber;
%qX's and qY's normalized to run between zero and one.

nsites = numel(X);
Density = zeros(nsites,1);
for ri = 0:nsites-1
        phaseSpDensity = exp( -(2*t*(2-cos(pi*qX)-cos(pi*qY))+trapdepth*(1-exp(-2*trapParam*ri*ri))-mu)/kT )./...
        (exp( -(2*t*(2-cos(pi*qX)-cos(pi*qY))+trapdepth*(1-exp(-2*trapParam*ri*ri))-mu)/kT ) +1)+...
        exp( -(2*t*(2-cos(pi*qX)-cos(pi*qY))+trapdepth*(1-exp(-2*trapParam*ri*ri))+wa-mu)/kT )./...
        (exp( -(2*t*(2-cos(pi*qX)-cos(pi*qY))+trapdepth*(1-exp(-2*trapParam*ri*ri))+wa-mu)/kT ) +1);
    Density(ri+1) = 4*sum(sum(phaseSpDensity));
end

DensityPerSite = Density*(pi/qpointnumber)^2/(2*pi)^2;
Ntot = 2*pi*sum(DensityPerSite.*[0.5:nsites-0.5]');
V = fudge*DensityPerSite;
end

