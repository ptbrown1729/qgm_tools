function signal_convolved = convolve_gauss(signal, dw, sigma)
% Convolve a signal with a gaussian of given width
%
% signal: signal evaluated at equally spaced points
%
% dw: distance between signal points
%
% sigma: standard deviation of gaussian in the same units as dw.

% get normalized gaussian
gauss = @(w) 1 / sqrt( 2*pi * sigma^2 ) * exp(- w.^2 / (2*sigma^2) );
% get symmetric sample point (including zero)
w_sample = 0 : dw : 3 * sigma;
w_sample = horzcat(-flip(w_sample), w_sample(2:end));
gauss_sample = gauss(w_sample);
center_index = (length(gauss_sample) - 1 ) / 2 + 1;
% normalize so integral over this finite region = 1
% gauss_sample = gauss_sample / (trapz(gauss_sample) * dw);
gauss_sample = gauss_sample / (trapz(gauss_sample));

if ndims(signal) > 2 || size(signal, 2) > 1
    % if our signal has several dimensions, assume the first is the
    % frequency dimension. Then convolve along each axis independently.
    signal_convolved = zeros( size(signal) );
    for ii = 1 : numel( signal(1, :) )
        signal_convolved(:, ii) = convolve_gauss(signal(:, ii), dw, sigma);
    end
 
else
    % do convolution
    signal_convolved = conv(gauss_sample, signal);

    % clip so that signal_convolved is at the same points that signal was.
    signal_convolved = signal_convolved(center_index : end - center_index + 1);
end

end