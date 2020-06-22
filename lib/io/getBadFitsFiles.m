function [BadFNames] = getBadFitsFiles(Path)
%function [FNames] = rmBadFitsFiles(FolderPath)
%Edit this function to check for correct exception...
Files = dir(fullfile(Path,'*.fits'));

BadFNames = {};
for ii =1:length(Files)
    FName = Files(ii).name;
    try
        fitsread(fullfile(Path,FName));
        %readimg(fullfile(Path),FName);
    catch MException
        %if strcmp(MException.identifer,'MATLAB:imagesci:fits:libraryError')&&strcmp(MException.message,'CFITSIO library error (108): error reading from FITS file')
            BadFNames = cat(2,BadFNames,FName);
        %end
    end
end


end

