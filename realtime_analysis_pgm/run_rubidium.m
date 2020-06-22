function p = run_rubidium(FolderNames, run_subfolders)
    % convenience function for rubidium experiment
    
    if ~exist('FolderNames', 'var')
        FolderNames = '';
    end
    
    if ~exist('run_subfolders', 'var')
        run_subfolders = 1;
    end
    
    root_dir = fullfile('//', '128.112.86.75', 'lithium',...
                        'MoleculeExperiment', 'IMAGING DATA (MOLECULES');
    constants = ConstantsClass('rubidium-87');
    p = ProgramClass(FolderNames, run_subfolders, constants, root_dir);

end
