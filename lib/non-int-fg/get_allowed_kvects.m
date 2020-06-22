function [kxkx, kyky] = get_allowed_kvects(nsites_x, nsites_y, mode)
% Get the allowed k-vectors for a periodic lattice system with N sites.
% Return the k-vectors in the range (-pi, pi]

if ~exist('mode', 'var') || isempty(mode)
    mode = 'balanced';
end

allowed_modes = { 'balanced', 'positive'};
fn = @(c) strcmp(c, mode);
if ~any(cellfun(fn, allowed_modes))
    error('mode argument was not one of the allowed keywords.');
end

kxs = (2*pi) / nsites_x * (0 : nsites_x - 1 );
if strcmp(mode, 'balanced')
    kxs(kxs > pi) = kxs(kxs > pi) - 2*pi;
    kxs = sort(kxs);
end

kys = (2*pi) / nsites_y * (0 : nsites_y - 1 );
if strcmp(mode, 'balanced')
    kys(kys > pi) = kys(kys > pi) - 2*pi;
    kys = sort(kys);
end

[kyky, kxkx] = ndgrid(kys, kxs);

end