function [intensity] = getIntensity(counts, image_time, wavelength, qe, counts_per_photon, pix_area)
%[intensity] = getIntensity(counts, image_time, wavelength, qe, counts_per_photon, pix_area)
% counts: number of camera electron counts
%
% image_time: exposure time in seconds
%
% wavelength: wavelength of light in meters
%
% qe: quantum efficiency of camera at given wavelength. This is the
% efficiency with which a photon is converted into counts
%
% counts_per_photon: the number of `counts' registered by the camera
% corresponding to a single photon
%
% pix_area: effective area of the pixel, i.e. accounting for magnification
% and etc.
%
%SI units, i.e. W/m^2.
%6Li D2 Isat = 25.4 W/m^2.
if ~exist('counters_per_photon', 'var') || isempty(counts_per_photon)
    counts_per_photon = 1;
end

hbar = 1.0546e-34;
c = 299792458;

w = 2 * pi * c / wavelength;

power = (counts * hbar * w) / qe / counts_per_photon;
intensity = power / pix_area / image_time;