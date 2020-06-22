function [NewPath] = timeStampFitsFile(FitsPath)
%[NewPath] = timeStampFitsFile(FitsPath)
%rename fits file with File modified date from FITS file metadata

try
    info = fitsinfo(FitsPath);
    TimeStamp = getAIAStyleTimeStamp(info.FileModDate);
    
    [Pathstr,~,Ext] = fileparts(FitsPath);
    NewPath = fullfile(Pathstr,cat(2,TimeStamp,Ext));
    
    if strcmp(FitsPath,NewPath)
    else
        movefile(FitsPath,NewPath)
    end
catch Exception
    MsgString = getReport(Exception);
    disp(MsgString);
    
end
end

