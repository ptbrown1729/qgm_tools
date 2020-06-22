function writeLog(LogPath, Vals, Keys, Stamp)
%writeLog(LogPath,Vals,Keys,Stamp)
%Writes log file. One line per call. Format is
%'timestamp,keyword,value,keyword,value,keyword \n'
%%%Arguments%%%
%Vals is a list of numers
%Keys is a cell array of strings
%Stamp is a string


if length(Vals)~=length(Keys)
    error('Length of vals and keywords does not match.');
end

try
    fid=fopen(LogPath,'a');
    fprintf(fid, '%s,%s,', char(Stamp), 'stamp');
    if ~isempty(Vals)
        for ii=(1:length(Keys) - 1)
            fprintf(fid, '%e,%s,', Vals(ii), char(Keys(ii)));
        end
        fprintf(fid, '%e,%s\n', Vals(end), char(Keys(end)));
        fclose(fid);
    else
        fprintf(fid, '\n');
        disp('No vals, keywords data pairs')
    end
    
catch Exception
    MsgString = getReport(Exception);
    disp('Writing to log file failed in writeLog.m');
    disp(MsgString);
    
    try
        %fclose can still fail, so wrap it in try...catch
        fclose(fid);
    catch Exception
        MsgString = getReport(Exception);
        disp('Closing log file failed in writeLog.m');
        disp(MsgString);
    end
end
end

