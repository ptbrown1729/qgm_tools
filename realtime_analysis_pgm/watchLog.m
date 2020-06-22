function [h] = watchLog(Paths, Index, Lims, StartPt_StepSizePairs)
%[h] = watchLog(Paths,Index,Lims,StartPt_StepSizePairs)
%Monitors a log file, which is a collection of value, keyword pairs.
%Continuously check for the pairs listed in Index variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paths is either a cell array of different paths to watch or a number
% specifying the data set under todays date to analyze.
%
% Index is an array of indices to watch in the log file.
%
% Lims is an nx2 array of plot limits for each log file index.
%
% StartPt_StepSizePairs is a nx2 array describing the
%

% check arguments
if ischar(Paths)
    Paths = {Paths};
end

if exist('StartPt_StepSizePairs', 'var')
    if size(StartPt_StepSizePairs, 1) == 1
        StartPt_StepSizePairs = repmat(StartPt_StepSizePairs, length(Index));
    end
end

% useful information
NImgs = length(Index);
NRows = ceil(sqrt(NImgs));
NCols = ceil(NImgs/NRows);

% structure to store figure handles
fighandles = cell(1, length(Paths));
figaxes = cell(1, length(Paths));
vals_data = cell(1, length(Paths));
keys_data = cell(1, length(Paths));

% create figures
for ii = 1:length(Paths)
    fighandles{ii} = figure;
end

Continue = 1;
while(Continue)
    
    for ii = 1:length(Paths)
        
        % get paths
        if iscell(Paths)
            Path = char(Paths{ii});
        elseif ismatrix(Paths)
            Path = fullfile(getDataSetPath(today, Paths(ii)), 'log.txt');
        else
            Path = '';
        end
        
        [dir_name, ~, ~] = fileparts(Path);
        
%         FigureName = strcat('Watching Log File ',Path);
%         fhandle = findobj('type','figure','name',FigureName);
%         if isempty(fhandle)
%             figure('name',FigureName);
%         elseif length(fhandle)>1
%             close(fhandle)
%             figure('name',FigureName);
%         else
%             figure(fhandle);
%         end
        
        %stop button to exit without ctrl+c
        %uicontrol('parent', fhandle, 'Style', 'Pushbutton', 'String',...
        %'Stop', 'units', 'normalized', 'position', [0,0.5,0.4,0.2],...
        %'Visible', 'on', 'callback', @setContinueFalse);
        %c = uicontrol('Style', 'pushbutton', 'String', 'Stop',...
        %'Position', [20,20,60,20], 'Callback', 'Continue = 0');
        
        if exist(Path, 'file') == 2
            [Vals, Keys] = parseLog(Path);
            
            % if data has not changed, do not update our figure
            if ~isequal(Vals, vals_data{ii})
                vals_data{ii} = Vals;
                keys_data{ii} = Keys;

                % make our desired figure active
                if ~isvalid(fighandles{ii})
                    fighandles{ii} = figure;
                else
                    set(0, 'currentfigure', fighandles{ii});
                end
                
                for jj = 1:length(Index)

                    WatchedVar = Vals(:, Index(jj));
                    WatchedVarName = Keys(Index(jj));

                    if exist('StartPt_StepSizePairs', 'var')
                        StartPt = StartPt_StepSizePairs(jj, 1);
                        StepSize = StartPt_StepSizePairs(jj, 2);
                        XPts = StartPt:StepSize:(StartPt + StepSize * (length(WatchedVar) - 1));
                    else
                        XPts = 1:length(WatchedVar);
                    end

                    % plot each variable
                    subplot(NCols ,NRows, jj);
                    plot(XPts, WatchedVar, 'b');
                    hold on;
                    plot(XPts, WatchedVar, 'bo-');
                    % average and standard evaition
                    avg = mean(WatchedVar);
                    stdev = std(WatchedVar);
                    plot([XPts(1), XPts(end)], [avg,avg], 'r');
                    plot([XPts(1), XPts(end)], [avg + stdev, avg + stdev], 'r--')
                    plot([XPts(1), XPts(end)], [avg - stdev, avg - stdev], 'r--')
                    % labels
                    grid on;
                    hold off;
                    ylabel(WatchedVarName);
                    title(sprintf('%s\n Mean: %.2E, SD: %.2E, NPts: %u, SD/Mean: %.2f%%',...
                        dir_name, avg, stdev, length(XPts), 100*stdev/avg));

                    if exist('Lims','var')
                        ylim(Lims(jj, :))
                    end

                    if length(XPts) > 1
                        xlim([XPts(1), XPts(end)]);
                    end

                    dcm_obj = datacursormode(fighandles{ii});
                    set(dcm_obj, 'UpdateFcn', @dataTipUpdateFn);

                end
            end
        else
            pause(2);
        end
    end
    pause(2);
end

end

%when have time, set this up so this is the ext update function for the
%data tip....it's the same as the normal one except with extra decimal places showing...
function output_txt = dataTipUpdateFn(obj, event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
output_txt = {['X: ',num2str(pos(1),8)],...
    ['Y: ',num2str(pos(2),8)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),8)];
end
end

 function setContinueFalse(source, callbackdata)
    
 end