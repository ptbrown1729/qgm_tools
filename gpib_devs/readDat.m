function [Data] = readDat(GPIB_Dev, CommandStr, IsNumber)
    if ~exist('IsNumber','var')
        IsNumber = 0;
    end
    
    %write command to device
    fprintf(GPIB_Dev, CommandStr); 
    fprintf('Sent %s \n', CommandStr);
    %read response
    readasync(GPIB_Dev) 
    stat = 'read';

    while ~strcmp(stat, 'idle')
        stat = GPIB_Dev.TransferStatus;
%         fprintf('still reading \n')
        pause(1);
    end

    Data = fscanf(GPIB_Dev);
    fprintf('Received %d bytes \n', length(Data));
    if IsNumber
        Data = str2double(strsplit(Data, ','));
        Data = Data(~isnan(Data));
    end
	
    clrdevice(GPIB_Dev);
end