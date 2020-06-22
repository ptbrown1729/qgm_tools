function overlap_area = get_pixel_overlap(xo, lo, x1, l1)
%   Compute the overlap between pixels of different sizes lo and l1
%   centered at xo and x1 respectively.
%
%   arguments:
%   -----------------
%   xo: position of the first pixel. This may be an array of arbitrary
%   size.
% 
%   lo: full width of the first pixel
%
%   x1: position of the second pixel. This may be n array of arbitrary
%   size.
%
%   l1: full width of the second pixel
%
%   return:
%   -----------------
%   overlap_area: area of the two pixels which overlaps.

if numel(xo) == 1 && numel(x1) > 1
    xo = xo * ones( size(x1) );
elseif numel(xo) > 1 && numel(x1) == 1
    x1 = x1 * ones( size(xo) );
end

if ~isequal( size(xo), size(x1))
    error();
end


half_width_o = 0.5 * lo;
half_width_1 = 0.5 * l1;

dx = abs(xo - x1);

% if dx > half_width_o + half_width_1
%     overlap = 0;
% else
%     max_size = 2 * min(half_width_o, half_width_1);
%     overlap = min(half_width_o + half_width_1 - dx,  max_size);
% end  
max_size = 2 * min(half_width_o, half_width_1);
overlap_area = min(half_width_o + half_width_1 - dx, max_size);
overlap_area( dx > half_width_o + half_width_1 ) = 0;

end