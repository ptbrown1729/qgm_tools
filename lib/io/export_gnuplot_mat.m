function export_gnuplot_mat(fname, mat, x, y)
% Export a matrix to a text file in an appropriate format to be plotted
% using gnuplot. The syntax would be "plot fname matrix nonuniform with image"
% In this format the first column gives the y-coordinates of the image, the
% top row gives the x-coordinates of the image, and the remaining entries
% define the image.
%
% fname: name of the file
%
% mat: a 2D matrix
%
% x: x-coordinates, same size as mat
%
% y: y-coordinates, same size as mat

% TODO: shouldn't I use a meshgrid in this case?
if (~exist('x', 'var') || isempty(x)) && (~exist('y', 'var') || isempty(y))
    dlmwrite(fname, mat, '\t');
else
    
    N = length(x) + 1;
    data = horzcat(N, x(:)');
    data = vertcat(data, horzcat(y(:), mat));
    dlmwrite(fname, data, 'delimiter', '\t');
end

end