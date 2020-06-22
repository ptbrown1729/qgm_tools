function color_map = get_color_map(npts, start_color, center_color, end_color)
% produce a colormap consisting of n-different colors, which interpolates
% from start_color, to center_color, and finally to end_color
%
% npts: number of colors to produce
%
% start_color: A 1 x 3 vector representing RGB value
%
% center_color: A 1 x 3 vector giving center color of the colormap
%
% end_color:
%
% Returns:
%
% color_map: A npts x 3 array where each row represents a color

if ~exist('npts', 'var') || isempty(npts)
    npts = 200;
end

if ~exist('start_color', 'var') || isempty(start_color)
    start_color = [0, 0, 1]; % blue
end

if ~exist('center_color', 'var') || isempty(center_color)
    center_color = [1, 1, 1]; % white
end

if ~exist('end_color', 'var') || isempty(end_color)
    end_color = [1, 0, 0]; % red
end

fn_StartToCenter = @(t) (1-t)*start_color + t*center_color;
fn_CenterToEnd = @(t) (1-t)*center_color + t*end_color;
fn = @(t) fn_StartToCenter(2*t)*(t<0.5) + fn_CenterToEnd(2*(t-0.5))*(t>0.5);
ts = linspace(0, 1, npts);

color_map = zeros(npts, 3);
for ii = 1:npts
    color_map(ii, :) = fn(ts(ii));
end

end