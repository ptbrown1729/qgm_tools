function [qxs, qys, sfact, sfact_unc] = get_structure_fact(array, array_err, mode, use_fft, qxs, qys)
% compute structure factor for a given array, or a quarter of an array.
%
% array: If array has more than 2 dimensions, assume the first two dimensions are
% real space dimensions, and FFT these, leaving the others in place.
%
% array_err: uncertainty associated with array
% 
% mode: "full" or "quarter". If "quarter", then assumes you are supplying
% only the 'southeast' corner of the full matrix (i.e. x>=0 and y>=0
% portion). Further, assume that this quadrant includes the 'center'
% coordinates, x=0, y=0 (we are typically imaging array is a correlation
% matrix)

if ~exist('mode', 'var') || isempty(mode)
    mode = 'full';
end

if ~exist('use_fft', 'var') || isempty(use_fft)
    if ~exist('qxs', 'var') || ~exist('qys', 'var') || isempty(qxs) || isempty(qys)
        use_fft = 1;
    else
        use_fft = 0;
    end
end

if use_fft && (exist('qxs', 'var') || exist('qys', 'var'))
    error('use_fft=1, but qxs and qys were supplied. Either set use_fft=0, or do not supply qxs and qys.');
end

% get allowed qs
if strcmp(mode, 'full')
    nxqs = size(array, 2);
    nyqs = size(array, 1);

    cx = (nxqs + 1)/2;
    cy = (nyqs + 1)/2;
    [xx, yy] = meshgrid(1:nxqs, 1:nyqs);
    xx = xx - cx;
    yy = yy - cy;
elseif strcmp(mode, 'quarter')
    nxqs = 2 * size(array, 2) - 1;
    nyqs = 2 * size(array, 1) - 1;
    
    [xx, yy] = meshgrid(0:size(array, 2) - 1, 0:size(array,1)-1);
end

if ~exist('qxs', 'var') || ~exist('qys', 'var')
    qxs = 2*pi/nxqs * (0:nxqs-1);
    qys = 2*pi/nxqs * (0:nyqs-1);

    qxs(abs(qxs - pi) < 1e-15) = pi;
    qys(abs(qys - pi) < 1e-15 ) = pi;

    qxs(qxs >= pi) = qxs(qxs >= pi) - 2*pi;
    qys(qys >= pi) = qys(qys >= pi) - 2*pi;

    qxs = fftshift(qxs);
    qys = fftshift(qys);

    if strcmp(mode, 'quarter')
        qxs = qxs(qxs >= 0);
        qys = qys(qys >= 0);
    end
end

[qxqx, qyqy] = meshgrid(qxs, qys);

% get structure factor
sfact = zeros(size(qxqx));
if strcmp(mode, 'full')
    for ii = 1 : numel(array(1, 1, :))
        % ifftshift accounting for fact that we (x,y) = (0,0) to be the center
        % of the array.
        if use_fft
            sfact(:, :, ii) = fftshift(fft2(ifftshift(array(:, :, ii))));
        else
            for aa = 1 : size(qxqx, 1)
                for bb = 1 : size(qxqx, 2)
                    sfact(aa, bb, ii) = sum(sum(array(:, :, ii) .* exp(-1i * qxqx(aa, bb) * xx  - 1i * qyqy(aa, bb) * yy))); 
                end
            end
        end
    end
elseif strcmp(mode, 'quarter')
    xx_full = quad2full(xx, 'southeast', 1);
    xx_full(:, 1:size(xx_full,2)/2+1) = -xx_full(:, 1:size(xx_full,2)/2+1);
    yy_full = quad2full(yy, 'southeast', 1);
    yy_full(1:size(yy_full,1)/2+1, :) = -yy_full(1:size(yy_full,1)/2+1, :);
    
    for ii = 1 : numel(array(1, 1, :))
        [full_arr, ~, ~] = quad2full(array(:, :, ii), 'southeast', 1);
        if use_fft 
            sf = fftshift(fft2(ifftshift(full_arr)));
            [sfact(:, :, ii), ~, ~] = get_quadrant(sf, 'southeast', 'include_center');
        else
%             sf = zeros(size(qxqx));
            for aa = 1 : size(qxqx, 1)
                for bb = 1 : size(qxqx, 2)
                    sfact(aa, bb, ii) = sum(sum(full_arr(:, :, ii) .* exp(-1i * qxqx(aa, bb) * xx_full  - 1i * qyqy(aa, bb) * yy_full))); 
                end
            end
        end
      end
else 
    error('mode was not an allowed argument');
end

if exist('array_err', 'var') && ~isempty(array_err)
    % get uncertainty
    %Can write S(q) = sum_{ix,iy>=0) 2*cos(qx*ix + qy*iy) <n_i*n_j>_c
    %=> S(q)_Unc = sqrt(sum_{ix,iy>=0} 2*cos(qx*ix + qy*iy) * sigma_ij^2)
    % probably faster ways of doing this ... hopefully this is ok for now
    % [ixix, iyiy] = meshgrid( 1 : size(array, 2), 1: size(array, 1));
    sfact_unc = zeros(size(qxqx));
    
    for ii = 1: numel(array_err(1, 1, :))
        for aa = 1 : size(qxqx, 1)
            for bb = 1 : size(qxqx, 2)
                sfact_unc(aa, bb, ii) = sqrt( sum(sum(4*cos(qxqx(aa, bb) * xx + qyqy(aa, bb) * yy).^2 .* array_err(:, :, ii).^2 .* (yy >=0), 2), 1));
            end
        end
    end
else
    sfact_unc = zeros(size(sfact));
end
