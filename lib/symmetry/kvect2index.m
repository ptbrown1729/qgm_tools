function [indices_gxmg, indices_gymg] = kvect2index(kxs, kys, tolerance, xmax, ymax)
% given arrays of kxs and kys, determine which points are on the
% Gamma-X-M-Gamma and Gamma-Y-M-Gamma lines within a certain tolerance.
% Return arrays of indices, where points on the Gamma-X line are assigned a
% value in the interval [0,1], points on the X-M line [1,2], and points on
% the M-Gamma line [2,3]
%
% kxs: array of kx vectors
%
% kys: array of ky vectors
%
% tolerance: Amount we allow the k-point to be off the line for us to count
% it. 
%
% xmax: maximum allowed value on the x-direction. Typically pi.
%
% ymax: maximum allowed value on the y-direction. Typically pi.

if ~exist('tolerance', 'var') || isempty(tolerance)
    tolerance = 0;
end

if ~exist('xmax', 'var') || isempty(xmax)
    xmax = pi;
end

if ~exist('ymax', 'var') || isempty(ymax)
    ymax = pi;
end


indices_gxmg = nan(size(kxs));
indices_gymg = nan(size(kxs));

on_gx_line = abs(kys - 0) <= tolerance & (kxs >= -tolerance) & (kxs <= xmax + tolerance);
on_xm_line = abs(kxs - pi) <= tolerance & (kys >= -tolerance) & (kys <= xmax + tolerance);

on_gy_line = abs(kxs - 0) <= tolerance & (kys >= -tolerance) & (kys <= ymax + tolerance);
on_ym_line = abs(kys - ymax) <= tolerance & (kxs >= -tolerance) & (kys <= ymax + tolerance);

on_mg_line = abs(kys - kxs) <= tolerance & (kys >= -tolerance) & (kys <= ymax + tolerance);

% do mg first so that 0 takes precedence over 3
indices_gxmg(on_mg_line) = 3 - kxs(on_mg_line) / xmax;
indices_gxmg(on_gx_line) = kxs(on_gx_line) / xmax; 
indices_gxmg(on_xm_line) = kys(on_xm_line) / ymax + 1;

indices_gymg(on_mg_line) = 3 - kxs(on_mg_line) / xmax;
indices_gymg(on_gy_line) = kys(on_gy_line) / ymax;
indices_gymg(on_ym_line) = kxs(on_ym_line) / xmax + 1;

end