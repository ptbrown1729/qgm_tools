function [YSize,XSize,Layers] = getImageSize(FilePath)
%[YSize,XSize,Layers] = getImageSize(FilePath)
% TODO: can I get rid of thsi function?
    [Path,FName,Ext] = fileparts(FilePath);
    if exist(FilePath,'file')
    
    if strcmp(Ext,'.aia')
        try
            fid = fopen(FilePath);
            fread(fid,1);fread(fid,1);fread(fid,1);fread(fid,1);fread(fid,1); %unused .aia lines
            YSize=fread(fid,1,'uint16'); %rows
            XSize=fread(fid,1,'uint16'); %columns
            Layers=fread(fid,1,'uint16'); 
        catch
            if ~exist('YSize','var') || ~exist('XSize','var') || ~exist('Layers','var')
                YSize = 0; XSize = 0; Layers = 0;
            end     
        end
        fclose(fid);
        
    elseif strcmp(Ext,'.fits')
        try
            Info = fitsinfo(FilePath);
            Sizes = Info.PrimaryData.Size;
            YSize = Sizes(1); XSize = Sizes(2); Layers = Sizes(3);
        catch
            YSize = 0; XSize = 0; Layers = 0;
            disp('Problem with fitsinfo...')
        end
        
    else
        YSize = 0; XSize = 0; Layers = 0;
        disp('File extension not supported');
    end
    else
        disp('FilePath did not exist...')
        YSize = 0; XSize = 0; Layers = 0;
    end
end