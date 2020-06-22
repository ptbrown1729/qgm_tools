function [images, layers, stamp, vals, keys] = readimg(file_path)
% Read .aia or .fits file and return images and metadata.
%
% Images: A stack of images stored in the file
%
% Layers: The number of images. This is redundant because it can be
% obtained from size(Images, 3), but leaving it in for backwards
% compatibility
%
% Stamp: file name.
%
% Vals: Metadata is assumed to be a set of value, key pairs. This is an
% array of the vallues
%
% Keys: Cell array of keys associated with values.

[~, name, ext] = fileparts(file_path);
% if strcmp(FileName(end - 2:end), 'aia')
if strcmp(ext, '.aia')
    
    %Open file and parse filename. Exception handling to close file in case
    %of error.
    try
        
        fid = fopen(file_path);
        stamp = name;
        
        if length(stamp) ~= 19
            % TODO: this is an awful way to try and do this check. Use a
            % regular expression, or better yet, don't bother with this
            % check...
            disp('DateTime Stamp is not of the format yyyy-mm-dd;hh;mm;ss.')
        end
        
        %Parse .aia file
        %First three are 'AIA'
        %add check to verify that...
		fread(fid,1);
        fread(fid,1);
        fread(fid,1);
		fread(fid,1);
        fread(fid,1); %unused .aia lines
        r = fread(fid, 1, 'uint16'); %rows
        c = fread(fid, 1, 'uint16'); %columns
        layers = fread(fid, 1, 'uint16'); %layers
        %3D array holding images.
        images = zeros(r, c, layers);

        for ii = 1:layers
            for jj = 1:r
                images(jj, :, ii) = fread(fid, c, 'uint16');
            end
        end
 
        % read annotation at the end of file. This is in format ASCII 
        % character terminated by null character, then 64bit unsigned 
        % integer (double float).
        keys= cell(0, 0);
        vals = [];
        while ~feof(fid)
            curr = '';
            phrase = '';
            while ~strcmpi(char(0), curr) && ~feof(fid) %read until a null character occurs (char(0))
                curr = fread(fid, 1, '*char');
                phrase = strcat(phrase, curr);
            end
            if isempty(phrase)
                % to fix problem where after reading last number but before 
                % eof enters second while loop and puts empty strings on 
                % the end of everything.
            else
                cellphrase = cell(1 ,1); %create cell array then put phrase in it. Maybe there's a better way of dealing with these.
                cellphrase{1} = phrase;
                keys=cat(2, keys, cellphrase);
                newval=fread(fid, 1, 'double', 0); %then read off the number
                vals=cat(2, vals, newval);
            end
        end
        
        %close file
        fclose(fid);
        
        if length(keys) ~= length(vals)
            error('Keywords and vals of different lengths.')
        end
        
    catch err
        fclose(fid);
        disp(err.message)
    end
    
elseif strcmp(ext, '.fits')
    try
        data = fitsread(file_path);
        %Andor software number y-axis of arrays from the bottom. Matlab starts
        %from the top.
        layers = size(data, 3);
        images = flipdim(data, 1);
        stamp = file_path;
%%%get metadata...this part works, but not sure that I want to print all of this info all of the time
%so want to fix that first.
%         info = fitsinfo(FileName);
%         Keys = info.PrimaryData.Keywords(:,1);
%         Vals = info.PrimaryData.Keywords(:,2);
%         IsNumeric = cellfun(@(x) isnumeric(x),Vals);
%         Keys = Keys(IsNumeric);
%         Vals = cell2mat(Vals(IsNumeric));
    catch err
        disp(err.message)
        layers = 0;
        images = NaN;
        stamp = '';
    end
    vals=[];
    keys={};
end
