function [mx, my, weight] =  get_moment(img, order, cxs, cys)
% Compute arbitrary moments of a 2d distribution/image or a stack of such
% images. Assuming that the center of the pixel img(j, i) is located at the
% coordinate (x, y) = (i, j).
%
% [mx, my, weight] =  get_moment(img, order, cxs, cys)
%
% img: a single (2D) image, or a stack of images (ny x nx x nimgs matrix)
% this function computer the moments of these images separately.
%
% order: the order of the moments to compute. order = 1 is the center of
% mass, order = 2 is the variance, etc.
%
% cxs: a list of center coordinates along the x-axis to be used in each
% image in the stack. If only a single value is supplied, it will be
% applied to all images. If omitted, the center of masses will be used.
%
% cys: a list of center coordinates along the y-axis to be used in each
% image in the stack. If only a single value is supplied, it will be
% applied to all images. If omitted, the center of masses will be used.
%
% definition used to compute moments:
%
% mx = \frac{ sum_{ij} img(i, j) * ( j - cx )^(order) }{ weight }
%
% my = \frac{ sum_{ij} img(i, j) * ( i - cy )^(order) }{ weight }
%
% weight = \sum_{ij} img(i, j)

if ~exist('cxs', 'var') || isempty(cxs) || ~exist('cys', 'var') ||  isempty(cys)
    if order == 1
        cxs = 0;
        cys = 0;
    else    
        [cxs, cys, ~] = get_moment(img, 1, 0, 0);
    end
end

nimgs = size(img, 3);
nxs = size(img, 2);
nys = size(img, 1);
[xxs, yys, ~] = meshgrid(1:nxs, 1:nys, 1:nimgs);

% ensure cxs and cys are the correct size.
if nimgs > 1
    % if cxs is a list of numbers
    if length(cxs) > 1 && ndims(cxs) < 3
        cxs = repmat( cxs(:), [1, nys, nxs]);
        cxs = permute(cxs, [2, 3, 1]);
    % if cxs is a single number
    elseif length(cxs) == 1
        cxs = cxs * ones(nys, nxs, nimgs);
    else
        error('cxs wrong shape get_moment.m');
    end
    
    % if cxs is a list of numbers
    if length(cys) > 1 && ndims(cys) < 3
        cys = repmat( cys(:), [1, nys, nxs]);
        cys = permute(cys, [2, 3, 1]);
    % if cxs is a single number
    elseif length(cys) == 1
        cys = cys * ones(nys, nxs, nimgs);
    else
        error('cys wrong shape get_moment.m');
    end
end

xxs = xxs - cxs;
yys = yys - cys;

weight = squeeze( sum( sum(img, 2), 1) );
mx = squeeze( sum( sum(img .* (xxs .^ order), 2), 1) ) ./ weight;
my = squeeze( sum( sum(img .* (yys .^ order), 2), 1) ) ./ weight;

    
end