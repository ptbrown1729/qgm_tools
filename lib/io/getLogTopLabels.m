function getLogTopLabels(fpath)
%getLogTopLabels(fpath)
%This is not the preferred way to save log files, since it cannot deal with
%input information of different lengths. Useful for creating txt files that
%can be easily imported to spreadsheet programs.

[pathstr, ~, ~] = fileparts(fpath);
spath=fullfile(pathstr, 'logTopLabels.txt');

[valsArray, keywords] = parseLog(fpath);
[rows, cols] = size(valsArray);

try
    %write Labels on top
    fid=fopen(spath, 'w');
    for jj=1:cols - 1
        fprintf(fid, '%s,', char(keywords(jj)));
    end
    fprintf(fid, '%s\n', char(keywords(end)));
    
    %write values
    for ii=1:rows
        for jj=1:cols - 1
            fprintf(fid, '%e,', valsArray(ii, jj));
        end
        fprintf(fid, '%e\n', valsArray(ii, end));
    end
   
    fclose(fid);
catch Exception
    MsgString = getReport(Exception);
    disp(MsgString);
    disp('Writing to log file failed in getLogTopLabels.m.');
    
    try
        fclose(fid);
    catch Exception
        MsgString = getReport(Exception);
        disp(MsgString);
        disp('Closing log file failed in getLogTopLabels.m.');
    end
    
end

% [NumDataSets,NumDataPnts] = size(vals);
%
% NumShots = 0;
% for i = 2:NumDataPnts
%     if strcmp(keywords(i),'Cx')
%         NumShots = NumShots+1;
%     end
% end
%
% KeyWordsCondensed = repmat({''},1,9);
%
% for i = 1:9
%     KeyWordsCondensed(i) = keywords(i+1+(i-1)*NumShots);
% end

end