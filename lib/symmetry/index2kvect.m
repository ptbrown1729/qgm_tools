function [kxs, kys] = index2kvect(linear_indices, mode)
% Convert linear index giving position along GXMG or GYMG line to 2D
% k-vector index.
%
% linear_indices: an array of `linear indices' which are numbers in the
% range [0, 3], with 0 corresponding to the Gamma = (0,0) point, 1
% corresponding to the X = (pi, 0) point, 2 corresponding to M = (pi, pi),
% and 3 corresponding to Gamma again.
%
% mode: "gxmg" or "gymg"

if ~exist('mode', 'var')
    mode = 'gxmg';
end

allowed_modes = {'gxmg', 'gymg'};
fn = @(c) strcmp(c, mode);
if ~any(cellfun(fn, allowed_modes))
    error('mode argument was not one of the allowed keywords.');
end

kxs = zeros(size(linear_indices));
kys = zeros(size(linear_indices));

% recall that to map interval [a,b] -> [c,d], take
% f(x) = (x-a) * (d-c)/(b-a) + c

if strcmp(mode, 'gxmg')
    on_gx_cut = linear_indices >= 0 & linear_indices <= 1;
    kxs(on_gx_cut) = linear_indices(on_gx_cut) * pi;
    kys(on_gx_cut) = 0;
    
    on_xm_cut = linear_indices > 1 & linear_indices <=2;
    kxs(on_xm_cut) = pi;
    kys(on_xm_cut) = (linear_indices(on_xm_cut) - 1) * pi;
    
elseif strcmp(mode, 'gymg')
    on_gy_cut = linear_indices >= 0 & linear_indices <= 1;
    kxs(on_gy_cut) = 0;
    kys(on_gy_cut) = linear_indices(on_gy_cut) * pi;
    
    on_ym_cut = linear_indices > 1 & linear_indices <=2;
    kxs(on_ym_cut) = (linear_indices(on_ym_cut) - 1) * pi;
    kys(on_ym_cut) = pi;
else
    error('mode should be gxmg or gymg but was %s', mode);
end

% MG direction
on_mg_cut = linear_indices > 2 & linear_indices <= 3;
kxs(on_mg_cut) = (3 - linear_indices(on_mg_cut)) * pi;
kys(on_mg_cut) = kxs(on_mg_cut);


end