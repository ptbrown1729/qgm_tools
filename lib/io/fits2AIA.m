function [AIAPath] = fits2AIA(FitsPath)
%FITS2AIA takes a fits file and produces an AIA file.
%[AIAPath] = fits2AIA(FitsPath)
% TODO: add support for metadata

[Path, Name, Ext] = fileparts(FitsPath);

if strcmp(Ext, '.fits')

[Images, ~, ~, ~, ~] = readimg(FitsPath);
AIAPath = fullfile(Path,cat(2, Name, '.aia'));
saveAIA(Images, AIAPath, [], {});

else
    error('File did not have .fits extension');
   %throw exception 
end



end

