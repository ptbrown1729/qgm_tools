function [img_d4_avg, unc_d4, img_d4_std, img_d4_sem] = d4_average(img, img_sem)
% [img_d4_avg, unc_d4, img_d4_std, img_d4_sem] = d4_average(img, img_sem)
% d4_average Produces an average of a matrix assuming that matrix possess
% square symmetry. For example, the points at all four corners will be
% averaged together. We refer to this as 'd4 averaging' because the
% symmetry group of the square is the fourth dihedral group. It has 8
% symmetry operations: the identity, 90 degree rotation, 180 degree rotation,
% 270 degree rotation, reflection, reflection combined with 90 degree rotation,
% reflection combined with 180 degree rotation, and reflection combined with 270 degree
% rotation.
%
%   arguments:
%   -----------------
%   img: the 2D matrix to be average
%
%   img_sem: the standard errors of the mean associated with each point in
%   the image
%
%   outputs:
%   -----------------
%   img_d4_avg: the matrix average being D4 averaged. It is the same size
%   as the original matrix and therefore contains redundant information. An
%   eight of this matrix (i.e. half of a quadrant split along its diagonal)
%   contains the complete information.
%   
%   unc_d4: the weighted uncertainties associated with each pixel of
%   img_d4_avg. These are computed according to unc_I = 1 / sqrt( sum( 1 / sigma_i^2 ) )
%   where I refers to a pixel in the final image, and the i run over the
%   pixels which are averaged together
%
%   img_d4_std: The standard deviation of the four pixels which are being
%   averaged. This includes 'Bessels correction', i.e. we divide by 1/N-1 
%   instead of 1/N
%
%   img_d4_sem: The standard error of the mean of the eight pixels which are
%   bing averaged. This is equal to img_d4_std / sqrt(8).

if ~exist('img_sem', 'var') || isempty(img_sem)
    img_sem = zeros(size(img));
end

img_d4_avg = 0.125 * d4_sum(img);

% TODO: since some points are averaged with themselves, the uncertainty for
% those points should not be reduced.
% standard deviation and standard deviation of the mean
npix = 8;
img_d4_std = sqrt( d4_sum( (img_d4_avg - img).^2) / (npix - 1) );
img_d4_sem = img_d4_std / sqrt(npix);

% weighted average 
% ignore points with zero uncertainty for calculating uncertainty
unc_weight = 1 ./ img_sem .^2;
unc_weight( isinf(unc_weight) ) = 0;

unc_d4 =  1 ./ sqrt( d4_sum( unc_weight ) );

% handle case where none of the points in the ROI have finite uncertainty.
unc_d4( isinf(unc_d4 )) = 0;

end

function d4mat = d4_sum(mat)
    d4mat = mat + rot90(mat, 1) + rot90(mat, 2) + rot90(mat, 3) + ...
                transpose(mat) + rot90(transpose(mat), 1) + ...
                rot90(transpose(mat), 2) + rot90(transpose(mat), 3);
end



