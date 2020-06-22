function [h] = watchReconstruction(Paths,Index,Lims,StartPt_StepSizePairs)
%[h] = watchLog(Paths,Index,Lims,StartPt_StepSizePairs)
%Paths is a cell array of different paths to watch.
%Monitors a log file, which is a collection of value, keyword pairs.
%Continuously check for the pairs listed in Index variable.
%
%StartPt_StepSizePairs is a nx2 array.
close all;

FolderNum = [2016,9,7,27];
dataset = {FolderNum(1),FolderNum(2),FolderNum(3),FolderNum(4),1,2};

NImgs = length(Index);
NRows = ceil(sqrt(NImgs));
NCols = ceil(NImgs/NRows); 

if ischar(Paths)
    Paths = {Paths};
end

if exist('StartPt_StepSizePairs','var')
    if size(StartPt_StepSizePairs,1) == 1
        StartPt_StepSizePairs = repmat(StartPt_StepSizePairs,length(Index));
    end
end


Continue = 1;
while(Continue)
    
    for ii = 1:length(Paths)
        if iscell(Paths)
            Path = char(Paths{ii});
        elseif ismatrix(Paths)
            [filenames,~] = getFilenames(getDataSetPath(datenum(date),Paths(ii)));
%             Path = fullfile(getDataSetPath(datenum(date),Paths(ii)),'log.txt');
        else
            Path = '';
        end
        FigureName = strcat('Watching Log File ',Path);
        fhandle = findobj('type','figure','name',FigureName);
        if isempty(fhandle)
            figure('name',FigureName);
        elseif length(fhandle)>1
            close(fhandle)
            figure('name',FigureName);
        else
            figure(fhandle);
        end
        
        %stop button to exit without ctrl+c
        %uicontrol('parent',fhandle,'Style','Pushbutton','String','Stop','units','normalized','position',[0,0.5,0.4,0.2],'Visible','on','callback',@setContinueFalse);
        %c = uicontrol('Style','pushbutton','String','Stop','Position',[20,20,60,20],'Callback','Continue = 0');
        
        if exist(Path,'file') == 2
            [Vals,Keys] = parseLog(Path);
            
            for jj = 1:length(Index)
                
                WatchedVar = Vals(:,Index(jj));
                WatchedVarName = Keys(Index(jj));
                
                if exist('StartPt_StepSizePairs','var')
                    StartPt = StartPt_StepSizePairs(jj,1);
                    StepSize = StartPt_StepSizePairs(jj,2);
                    XPts = StartPt:StepSize:(StartPt+StepSize*(length(WatchedVar)-1));
                else
                    XPts = 1:length(WatchedVar);
                end
                
                Avg = mean(WatchedVar);
                Std = std(WatchedVar);
                subplot(NCols,NRows,jj);
                plot(XPts,WatchedVar,'r');
                hold on;
                plot(XPts,WatchedVar,'rx');
                plot([XPts(1),XPts(end)],[Avg,Avg],'b');
                plot([XPts(1),XPts(end)],[Avg+Std,Avg+Std],'g-')
                plot([XPts(1),XPts(end)],[Avg-Std,Avg-Std],'g-')
                grid on;
                hold off;
                ylabel(WatchedVarName);
                title(sprintf('Mean: %.2E, SD: %.2E, NPts: %u, SD/Mean: %.2f%%', Avg,Std,length(XPts),100*Std/Avg));
                
                
                if exist('Lims','var')
                    ylim(Lims(jj,:))
                end
                
                if length(XPts)>1
                    xlim([XPts(1),XPts(end)]);
                end
                
                dcm_obj = datacursormode(gcf);
                set(dcm_obj,'UpdateFcn',@dataTipUpdateFn);
                
            end
        else
        end
    end
    pause(2);
end

end
