function addAnalysis2Path(RootAnalysisPath)
%addAnalysis2Path(RootAnalysisPath)
%Add folders that ProgramClass() relies on to the matlab path.
SubPaths = {'','azimuthalavg','display','experiment','fitting',...
    'fitting/functions','fitting/functions/1D','fitting/functions/2D',...
    'Helper Functions','Helper Functions/Timestamps','Helper Functions/coordconversion','imageproc',...
    'physdata','statistics'};

%Folders on path.
PathCellArrayOfStrings = regexp(path,pathsep,'split');
%PathCellArrayOfStrings = vertcat(PathCellArrayOfCells{:});

%Check if directory is already on path, because adding to path takes a long
%time.
for ii = 1:length(SubPaths)
    CurrentPath = fullfile(RootAnalysisPath,SubPaths{ii});
    if ~any(strcmpi(CurrentPath,PathCellArrayOfStrings));
        addpath(CurrentPath);
    end
end

end

