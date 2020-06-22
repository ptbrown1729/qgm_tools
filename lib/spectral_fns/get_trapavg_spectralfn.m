function [spec_fn_trpavg, n_profile, areas] = get_trapavg_spectralfn(n_trap, areas_trap, spec_fn, n_samples, mode)
% densities_trap: real densities = nup + ndn
%
% areas_trap: area of the trap corresponding to each density
%
% spec_fn: an array representing a spectral function of arbitrary shape.
% Each dimension of this array corresponds to scanning some variable (e.g.
% temperature, interaction, k-vector, omega). The first dimension must correspond
% to density. Note that the spectral function array must be given at equal
% energies = omega + mu, as opposed to equal omegas.
%
% n_samples: Densities associated with one dimension of the
% spectral function
%
% mode: If "sum" then returns A_tr = \sum_n A(n) * Area(n). If "mean", then
% returns A_tr = \sum_n A(n) * Area(n) / \sum_n Area(n).
%
% returns:
%
% spec_fn_trap_avg:
%
% n_profile: 
%
% areas:

if ~exist('areas_trap', 'var') || isempty(areas_trap)
    areas_trap = ones( size(n_trap) );
end

if ~exist('mode', 'var') || isempty(mode)
    mode = 'sum';
end

% get density bins based on density samples provided.
if size(n_samples, 1) == 1
    n_end_pts = horzcat(0, n_samples, 2);
else
    n_end_pts = vertcat(0, n_samples, 2);
end

n_end_pts = 0.5 * (n_end_pts(2:end) + n_end_pts(1:end-1) );

% determine what area of our cloud falls between each density bin
[~, ~, area_fractions, ~, ~, npts, masks] = azAvg_General(areas_trap, [], n_trap, n_end_pts);
areas = area_fractions .* npts; 

% get the density profile of the the assumed cloud
n_samples_expanded = repmat(n_samples(:), [1, size(masks, 1), size(masks, 2)]);
n_samples_expanded = permute(n_samples_expanded, [2, 3, 1]);
n_profile = nansum( masks .* n_samples_expanded, 3);

% get the spectral function for the entire trap
spec_fn_size = size(spec_fn);
areas_expanded = repmat(areas, [1, spec_fn_size(2:end)]);
spec_fn_trpavg = zeros( spec_fn_size(2:end) );
spec_fn_trpavg(:) = nansum( areas_expanded .* spec_fn, 1);


if strcmp(mode, 'mean')
    spec_fn_trpavg = spec_fn_trpavg / nansum(areas);
elseif strcmp(mode, 'sum')
else
    error('mode was not one of the allowed values');
end

end