function [avg_vals, stddev_vals, nrepeats, keys] = getLogStats(LogPath, NCompare, NSkip)
%function [AvgVals,StdDevVals,NRepeats,Keys] = getStats(LogPath,NCompare)
%Does statistics on log file. Identifies duplicate entries by comparing
%final N fields.
%Even does 'average' of given values. So from AvgVals can extract all
%parameters of .aia file.
%
% %%ARGUMENTS%%
%
% LogPath is the path of the file to be processed
%
% NCompare is the number of fields to use when comparing different entries
% in the log file. For example, if NCompare = 2, only the last two vals of
% each entry will be considered. If these two vals are equal for different
% log entries, these entries will be averaged.
%
% NSkip is the number of lines to skip at the top of the log-file.
%
% %%OUTPUTS%%
%
% avg_vals
%
% stddev_vals
%
% nrepeats
%
% keys

[Vals, keys, ~] = parseLog(LogPath);

%modify so default bhavior is to average everything.
if ~exist('NCompare', 'var') || isempty(NCompare)
    NCompare = size(Vals, 2);
elseif NCompare == 0
    NCompare = size(Vals, 2);
end

if ~exist('NSkip', 'var') || isempty(NSkip)
    NSkip = 0;
elseif NSkip >= size(Vals, 1)
    NSkip = 0;
    warning('NSkip larger than number of files. Set to zero instead.');
end

CmpVals = Vals(1 + NSkip:end, end - NCompare + 1:end);
%ic contains an index for each row. If two rows have the same index, they
%are identical.
[~, ~, ic] = unique(CmpVals, 'rows');
%find number of unique indices.
UnqIndices = unique(ic);
NUnique = length(UnqIndices);
avg_vals = zeros(NUnique, size(Vals, 2));
stddev_vals = zeros(NUnique, size(Vals, 2));
nrepeats = zeros(NUnique, 1);

% now get statistics
for jj = 1:NUnique
    SingleKindVals = Vals(ic == UnqIndices(jj), :);
    avg_vals(jj, :) = mean(SingleKindVals, 1);
    stddev_vals(jj, :) = std(SingleKindVals, 0, 1);
    nrepeats(jj) = size(SingleKindVals, 1);
end

end

