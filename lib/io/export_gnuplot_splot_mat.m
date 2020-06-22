function export_gnuplot_splot_mat(fname, mat, x, y, z)
% Export a matrix to a text file in an appropriate format to be plotted
% using gnuplot. 
% 
% This file is expected to be lines of x, y, z, f(x,y,z) form with a
% newline between each row. This format is appropriate for data which
% should be displayed like a matrix, but which may not consist of
% uniformaly spaced points.
%
% The syntax would be e.g. "splot fname using 1:2:3:4 with lines".
%
% fname: name of the file
%
% mat: a 2D matrix
%
% x: x-coordinates, same size as mat
%
% y: y-coordinates, same size as mat
%
% z: z-coordinates, same size as mat

% TODO: is this the correct way to handle? Shouldn't I set x,y to a
% meshgrid in this case?
if (~exist('x', 'var') || isempty(x)) && (~exist('y', 'var') || isempty(y))
    dlmwrite(fname, mat, '\t');
else
    
if ~exist('z', 'var') || isempty(z) 
    z = zeros(size(mat));
end
    
% loop over rows.
for ii = 1 : size(x, 2)
    data = horzcat(x(:, ii), y(:, ii), z(:, ii), mat(:, ii));
    dlmwrite(fname, data, 'delimiter', ' ', '-append');

% newline between rows, but not after last row
    if ii < size(x, 2)
        fid = fopen(fname, 'a');
        fprintf(fid, '\n');
        fclose(fid);
    end
end
    
end

end