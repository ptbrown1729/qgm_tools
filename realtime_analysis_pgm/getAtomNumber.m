function [atom_num] = getAtomNumber(sum_od, lambda, pixel_area, detuning, gamma, intensity)
% getAtomNumber(sum_od, lambda, pixel_area, detuning, gamma, intensity)
%
% This function takes summed OD for an absorption image and detuning from
% resonance in Hz and returns the total atom number.
%
% Arguments:
% ------------------------
% sum_od: some of the optical depth over an image
%
% lambda: wavelength of the transition in meters.
%
% pixel_area: in meters.
%
% detuning: in radians/sec, i.e. should be specified as (2*pi)*frq
%
% gamma: in radians/sec
%
% intensity: in W/m^2

if ~exist('detuning', 'var') || isempty(detuning)
    detuning = 2 * pi * 0;
end

if ~exist('gamma', 'var') || isempty(gamma)
    gamma = 2 * pi * 6e6;
end

if ~exist('intensity', 'var') || isempty(intensity)
    intensity = 0;
end

% h = Constants.h;
% c = Constants.c;
% Isat = Constants.Isat_D2;
% gamma = Constants.gamma_D2;
% lambda = Constants.lambda_D2;
h = 6.62606957e-34; % J*s
c = 299792458; % m/s
Isat = pi * h * c / (3 * lambda^3) * gamma;

%define light scattering cross section
%resonance cross section, also = 3*lambda^2/(2*pi) = h * c * gamma / (lambda * 2 * Isat)
sigma_o = 3 * lambda^2 / (2*pi);
%light scattering scross section
sigma = sigma_o ./ (1 + (intensity / Isat) + 4 * detuning.^2 / gamma^2);

%Find the atom number. 
% number = pixel size * sum of od's / cross section
atom_num = pixel_area * sum_od ./ sigma';