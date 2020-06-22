function [ FitsPath ] = aia2Fits(AIAPath)
%AIA2FITS takes an AIA file and produces a Fits file.
%[FitsPath] = aia2Fits(AIAPath)
% TODO: add support for metadata

[Path, Name, Ext] = fileparts(AIAPath);

if strcmp(Ext, '.aia')
    
    [Images, ~, ~, Vals, Keys] = readimg(AIAPath);
    Images = flip(Images, 1); %so orientation when read by matlab agrees.
    FitsPath = fullfile(Path, cat(2, Name, '.fits'));
    fitswrite(Images, FitsPath);
    
    %should write Vals, Keys to fits metadata if possible.
    
else
    error('File did not have .aia extension');
    %throw exception
end

end

