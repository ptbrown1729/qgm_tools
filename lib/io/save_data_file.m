function save_data_file(array, titles_cell, delimiter, header_comment,...
                        folder, fname, append_file)
%   save_data_file(array, names_cell, delimiter, header_comment, folder, fname, append_file)
% writes a data file with an arbitrary header string, column headings, and
% columns of data delimited with delimiter. The column headings are described
% by the elements of titles_cell, which is a cell array of strings. If
% titles_cell contains only empty strings, the titles will not be written,
% and no extra line will be added to the file.
% The file is saved in folder
% with name fname. This file can be read using dlmread(fname, delimiter, n, 0)
% where n is the number of rows to be skipped at the beginning of the file.
% This will be the number of lines of the header_comment plus one line for
% titles.

if ischar(titles_cell)
    %catch simple case where only want one field and enter a
    %string.
    titles_cell = {titles_cell};
end

if size(array,2) ~= length(titles_cell)
    error('Number of columns in array was not equal to length of NamesCell');
end

if ~exist('delimiter', 'var') || isempty(delimiter)
    delimiter = '\t';
end

if ~exist('header_comment', 'var')
    header_comment = '';
end

if ~exist('folder','var') || isempty(folder)
    folder = pwd;
end

if ~exist('fname','var') || isempty(fname)
    root = 'test';
    fname = fullfile([root,'.txt']);
end

if ~exist('append_file','var') || isempty(append_file)
    append_file = 0;
end

fpath = fullfile(folder, fname);

if append_file
    FHandle = fopen(fpath, 'a');
else
    FHandle = fopen(fpath, 'w');
end

%write header comment
if ~isempty(header_comment)
    fprintf(FHandle, '%s', header_comment);
    if ~strcmp(header_comment(end-1:end), '\n')
        fprintf(FHandle, '\n');
    end
end

%write column headings
if ~all(cellfun('isempty', titles_cell))
    for ii = 1:length(titles_cell)-1
        string_lit = sprintf('%s%s', '%s', delimiter);
        fprintf(FHandle, string_lit, titles_cell{ii});
        %                     fprintf(FHandle,'%s\t', titles_cell{ii});
    end
    fprintf(FHandle,'%s\n', titles_cell{end});
end
fclose(FHandle);

%write data
dlmwrite(fpath, array, 'delimiter', delimiter, 'precision', 10, '-append');
end