% Collect quest DQMC output file data in matlab format.
% Ultimately, produce a set of interpolating functions which
% can be used for temperature/interaction fitting.
%
% 1. Read all quest output files and aggregate into an array
%    This is handled by the function "aggregate_dqmc_output_files"
%    Individual quest output files are read by "convert_qmc_file"
% 2. Reshape the quest data into a regular grid. This is handled by
%    the function "convert_dqmc_to_grid"
% 3. Convert from QMC natural variables (mu, T, U) to variables (n, T U)
%    and then create interpolating functions in those variables. This is
%    handled by "create_dqmc_interpolating_fn"

saving = 1;

file_pattern = '*.out';
dir_paths = {'Z:\QuestQMC\dqmc_data\2019_05_12\out',...
             'Z:\QuestQMC\dqmc_data\2019_05_13\out',...
             'Z:\QuestQMC\dqmc_data\2019_05_14\out'};

% dir_paths = {'Z:\QuestQMC\dqmc_data\VaryU\n=8,T=0.2to0.8,U=-6.5to-5.5,mu=-2to0.5,npass=50000\out'};
% dir_paths = {'Z:\QuestQMC\dqmc_data\VaryFilling\n=8,U=8,T=0.5to15_mu=-120to0,npass=50000\out'};
% dir_paths = {'Z:\QuestQMC\dqmc_data\VaryFilling\n=8,U=8,T=0.25to0.45,mu=-4.5to0,npass=50000\out'};
% dir_paths = {'Z:\QuestQMC\dqmc_data\VaryImbalance\n=8,T=0.5,U=-6to-4,mu=-8to8\Outputs',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\n=8,T=0.5,U=-8to-6,mu=-8to8\Outputs',...};
% dir_paths = {'Z:\QuestQMC\dqmc_data\VaryImbalance\new_n=8,T=0.7,U=-8to-4,mu=-8to8\output'};

% attractive hubbard model
% dir_paths = {'Z:\QuestQMC\dqmc_data\VaryImbalance\new_n=8,T=0.3,U=-8to-4,mu=-8to8\output',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\new_n=8,T=0.7,U=-8to-4,mu=-8to8\output',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\new_n=8,T=0.3to0.7,U=-2to-3.8,mu=-8to8\out',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\n=8,T=1.0to5.0,U=-2to-8,mu=-8to8\out',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\new_n=8,T=1.0to5.0,U=-2to-8,mu=-8to8\out',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\n=8,T=0.3to8.0,U=-2to-8,mu=-12to-8\out',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\n=8,T=10,U=-2to-8,mu=-12to8\out',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\2018_09_06\out',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\2018_09_07\out',...
%     'Z:\QuestQMC\dqmc_data\VaryImbalance\2018_11_19\out'};

% dir_paths = {...
% 'Z:\QuestQMC\dqmc_data\VaryFilling\n=8,U=8,T=0.25to0.45,mu=-4.5to0,npass=50000\out',...
% 'Z:\QuestQMC\dqmc_data\VaryFilling\n=8,U=8,T=0.5to15_mu=-120to0,npass=50000\out'};

% dir_paths = {...
% 'Z:\QuestQMC\dqmc_data\repulsive_hubbard_dqmc_fitting_results\2018_08_21\out'};

% first extract variables from all dqmc files into lists
% dqmc_results = aggregate_dqmc_output_files(dir_paths, file_pattern, [], [], [-10,10], [], [0,0]);
dqmc_results = aggregate_dqmc_output_files(dir_paths, file_pattern, [], [], [], [], []);
% reorganize dqmc fields into U, T, mu grids
dqmc_results_grid = convert_dqmc_to_grid(dqmc_results);

if saving
    name_root_dqmc_grid = sprintf('%s_gridded', dqmc_results_grid.identifier);
    while exist(fullfile([name_root_dqmc_grid, '.mat']), 'file')
        name_root_dqmc_grid = sprintf('%s_new', name_root_dqmc_grid);
    end
    fname_output_dqmc = fullfile([name_root_dqmc_grid, '.mat']);
    save(fname_output_dqmc, '-struct', 'dqmc_results_grid');
    fprintf('Saved gridded dqmc results to %s\n', fname_output_dqmc);
end

% create dqmc interpolating functions from dqmc parameter grids
interp_fns = create_dqmc_interpolating_fn(dqmc_results_grid);
if saving
    name_root_dqmc_interp = sprintf('%s_interp_fns', dqmc_results_grid.identifier);
    while exist(fullfile([name_root_dqmc_interp, '.mat']), 'file')
        name_root_dqmc_interp = sprintf('%s_new', name_root_dqmc_interp);
    end
    fname_output_interp = fullfile([name_root_dqmc_interp, '.mat']);
    save(fname_output_interp, '-struct', 'interp_fns');
    fprintf('Saved dqmc interpolating functions to %s\n', fname_output_interp);
end