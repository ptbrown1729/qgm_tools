classdef AttHubbClass < handle
    
    properties
        identifier
        date
        folders
        hold_time
        modulation_time
        using_flattening
        
        ones
        threes
        singles
        doubles
        
        dens
        dens_unc
        
        mag
        mag_unc
        
        pol
        pol_unc
        
        % 2d fit to total density
        gaussian_fit_params
        gaussian_fit_std_errs
        
        % azimuthal average of total density
        dens_azavg
        dens_azavg_unc
        dens_azavg_dist
        
        % sums of azimuthal averages of partial densities.
        dens13_azavg
        dens13_azavg_unc
        denssd_azavg
        denssd_azavg_unc
        
        dqmc_data_path
    end
    
    methods
        function obj = AttHubbClass(date_cell, folders, bin_endpts, exclude_pics, az_avg_mode, az_avg_grid, dqmc_data_path)
            
            if ~exist('az_avg_mode', 'var')
                az_avg_mode = '';
            end
            
            if ~exist('az_avg_grid', 'var')
                az_avg_grid = [];
            end
            
            if ~exist('dqmc_data_path', 'var') || isempty(dqmc_data_path)
                dqmc_data_path = fullfile('\\128.112.86.75\lithium\Publications\arpes_and_att_hubb_spin_suc\dqmc',...
                'dqmc_nsites=8x8_T=0.3-10.0_U=-14.5--2.0_muavg=-10.0-8.0_mudelta=0.0-0.0_npass=0-20000_interp_fns.mat');
            %                     'dqmc_nsites=8x8_T=0.3-10.0_U=-10.0--2.0_muavg=-10.0-8.0_mudelta=0.0-0.0_npass=0-20000_interp_fns.mat');
%                 'dqmc_nsites=8x8_T=0.3-10.0_U=-8.0--2.0_muavg=-12.0-8.0_mudelta=0.0-0.0_npass=0-20000_interp_fns.mat');
%                   'dqmc_nsites=8x8_T=0.3-8.0_U=-8.0--2.0_muavg=-8.0-8.0_mudelta=0.0-0.0_npass=5000-20000_interp_fns.mat');
              %                     'dqmc_nsites=8x8_T=0.3-5.0_U=-8.0--2.0_muavg=-8.0-8.0_mudelta=0.0-0.0_npass=5000-20000_interp_fns.mat');
            end
            
            obj.dqmc_data_path = dqmc_data_path;
            
            components = DataFolder( expandDataset({date_cell{:}, folders, 1, 1}),...
                                    exclude_pics, bin_endpts, az_avg_mode, az_avg_grid);
            ones = components(1);
            threes = components(2);
            singles = components(3);
            doubles = components(4);
            
            % magnetization
            diff = threes.Occs_ImgAvg - ones.Occs_ImgAvg;
            dens = threes.Occs_ImgAvg + ones.Occs_ImgAvg;
            unc = sqrt(threes.Occs_ImgAvgUnc .^2 + ones.Occs_ImgAvgUnc .^2);

            % azimuthal average of density
            init_params_gauss2d = [size(dens, 2)/2 , size(dens, 1)/2, 20, 20, 1, 0, 0];
            [fit_params_gauss2d, ~, ~, std_errs_gauss2d] = fit2D([], [],...
            dens, [], {'gaussian2D'}, init_params_gauss2d, zeros(size(init_params_gauss2d)));
            cx = fit_params_gauss2d(1);
            cy = fit_params_gauss2d(2);
            sx = fit_params_gauss2d(3);
            sy = fit_params_gauss2d(4);
            theta = fit_params_gauss2d(6);
            
             if sx > sy
                % ensure that the smaller direction is always the x
                % direction
                sy = fit_params_gauss2d(3);
                sy_err = std_errs_gauss2d(3);
                sx = fit_params_gauss2d(4);
                sx_err = std_errs_gauss2d(4);
                theta = theta + pi/2;
                
                fit_params_gauss2d(3) = sx;
                fit_params_gauss2d(4) = sy;
                fit_params_gauss2d(6) = theta;
                std_errs_gauss2d(3) = sx_err;
                std_errs_gauss2d(4) = sy_err;   
            end
            
            
            obj.gaussian_fit_params = fit_params_gauss2d;
            obj.gaussian_fit_std_errs = std_errs_gauss2d;
            
            az_avg_grid = ellipticalGrid(size(dens), [cx, cy, sx, sy, theta]);
            [dist, dist_sdm, dens_az_avg, dens_sdm, ~, dens_npts_az_avg, mask_stack] = ...
                azAvg_General(dens, [], az_avg_grid, ones.BinEdges);
            obj.dens_azavg_dist = dist;
            obj.dens_azavg = dens_az_avg;
            obj.dens_azavg_unc = dens_sdm;
            
            obj.folders = folders;
            obj.date = date_cell;
            
            obj.ones = ones;
            obj.threes = threes;
            obj.singles = singles;
            obj.doubles = doubles;
            
            obj.dens13_azavg = ones.Occs_AzAvg + threes.Occs_AzAvg;
            obj.dens13_azavg_unc = sqrt( ones.Occs_AzAvgUnc .^2 + threes.Occs_AzAvgUnc .^2);
            obj.denssd_azavg = singles.Occs_AzAvg + 2 * doubles.Occs_AzAvg;
            obj.denssd_azavg_unc = sqrt(singles.Occs_AzAvgUnc .^2 + 4 * doubles.Occs_AzAvgUnc.^2);
            
            obj.dens = dens;
            obj.dens_unc = unc;
            obj.mag = diff;
            obj.mag_unc = unc;
            obj.pol = diff ./ dens;
            obj.pol_unc = obj.pol * sqrt( (unc / diff) .^2 + (unc / dens) .^2);
            
            obj.identifier = sprintf('%04d-%02d-%02d-folders=%03d-%03d',...
                                     date_cell{1}, date_cell{2}, date_cell{3},...
                                     min(folders), max(folders));
            
        end
        
        function [fitp, std_err, mean_sqr_diff, fig_handle_mag, dens_mean, dens_std] = ...
                get_magnetization(obj, crop_radius, angle_grad, bin_width, bin_2d_size)
            if ~exist('crop_radius', 'var')
                crop_radius = 10;
            end
            
            if ~exist('angle_grad', 'var')
                angle_grad = -60 * pi / 180;
            end
            
            if ~exist('bin_width', 'var')
                bin_width = 2;
            end
            
            if ~exist('bin_2d_size', 'var')
                bin_2d_size = 1;
            end
            
            mag = binImg(obj.mag, bin_2d_size, bin_2d_size, 'Normal');
            mag_unc = sqrt(binImg(obj.mag_unc .^2 , bin_2d_size, bin_2d_size, 'Normal')) / bin_2d_size;
            dens = binImg(obj.dens, bin_2d_size, bin_2d_size, 'Normal');
            dens_unc = binImg(obj.dens_unc, bin_2d_size, bin_2d_size, 'Normal');
            pol = mag ./ dens;
            
                    % %crop off everything outside of a certain radius
%             cx =  0.5 * (obj.threes.Cx_AzAvg + obj.ones.Cx_AzAvg);
%             cy = 0.5 * (obj.threes.Cy_AzAvg + obj.ones.Cy_AzAvg);
%             aspect_rat = abs(obj.ones.GaussFitParams(3) / obj.ones.GaussFitParams(4));
%             angle_gauss = obj.ones.GaussFitParams(6);
            cx = obj.gaussian_fit_params(1);
            cy = obj.gaussian_fit_params(2);
            aspect_rat = abs(obj.gaussian_fit_params(3) / obj.gaussian_fit_params(4));
            angle_gauss = obj.gaussian_fit_params(6);
           
            
            
            [xx, yy] = meshgrid(1:size(obj.mag, 2), 1:size(obj.mag, 1));
            xx_cloud = (xx - cx) * cos(angle_gauss) - (yy - cy) * sin(angle_gauss);
            yy_cloud = (yy - cy) * cos(angle_gauss) + (xx - cx) * sin(angle_gauss);
            
            xx_cloud = binImg(xx_cloud, bin_2d_size, bin_2d_size, 'Normal');
            yy_cloud = binImg(yy_cloud, bin_2d_size, bin_2d_size, 'Normal');
            
            xx_bin = binImg(xx, bin_2d_size, bin_2d_size, 'Normal');
            yy_bin = binImg(yy, bin_2d_size, bin_2d_size, 'Normal'); 
            [~, cx_bin] = min(abs(xx_bin(1, :) - cx)); 
            [~, cy_bin] = min(abs(yy_bin(:, 1) - cy));
            
            radii = sqrt((xx_cloud).^2 / aspect_rat + (yy_cloud).^2 * aspect_rat);

            crop_mask = double(radii < crop_radius);
            nan_mask = crop_mask;
            nan_mask(nan_mask==0)=nan;
            mag_cropped = mag .* nan_mask;
            mag_unc_cropped = mag_unc .* nan_mask;
            dens_cropped = dens .* nan_mask;
            dens_unc_cropped = dens_unc .* nan_mask;

            dens_mean = nanmean(dens_cropped(:));
            dens_std = std(dens_cropped(~isnan(dens_cropped)));
            
            % sum along direction orthogonal to gradient
            xx_rot = (xx - cx) * cos(angle_grad) - (yy - cy) * sin(angle_grad);
            yy_rot = (yy - cy) * cos(angle_grad) + (xx - cx) * sin(angle_grad);
            
            xx_rot = binImg(xx_rot, bin_2d_size, bin_2d_size, 'Normal');
            yy_rot = binImg(yy_rot, bin_2d_size, bin_2d_size, 'Normal');

            % bins for averaging perpendicular to direction of applied
            % gradient
            bin_end_pts_x = [floor(min(xx_rot(:))):bin_width:ceil(max(xx_rot(:)))];
%             bin_centers = 0.5*(bin_end_pts_x(2:end) + bin_end_pts_x(1:end-1));
            bin_end_pts_y = [floor(min(yy_rot(:))):bin_width:ceil(max(yy_rot(:)))];
            % magnetization averages
            [x_pos, ~, mag_binned, mag_binned_unc, ~ ,npts_mag_binned, mask_stack] = ...
                azAvg_General(mag_cropped, [], xx_rot, bin_end_pts_x);
            [y_pos, ~, mag_binned_y, mag_binned_y_unc, ~, ~, mask_stack_y] = ...
                azAvg_General(mag_cropped, [], yy_rot, bin_end_pts_y);
            % average densities
            [~, ~, dens_binned, dens_binned_unc, ~, npts_dens_binned, ~] = ...
                azAvg_General(dens_cropped, [], xx_rot, bin_end_pts_x);
            [~, ~, dens_binned_y, dens_binned_y_unc, ~, ~, ~] = ...
                azAvg_General(dens_cropped, [], yy_rot, bin_end_pts_y);
            %Remove Nans
            npts_mag_binned = npts_mag_binned(~isnan(mag_binned));
            x_pos = x_pos(~isnan(mag_binned));
            mag_binned_unc = mag_binned_unc(~isnan(mag_binned));
            mag_binned_unc(mag_binned_unc == 0) = 1;
            mag_binned = mag_binned(~isnan(mag_binned));
            dens_binned = dens_binned(~isnan(dens_binned));
            dens_binned_unc = dens_binned_unc(~isnan(dens_binned_unc));
            
            y_pos = y_pos(~isnan(dens_binned_y));
            mag_binned_y_unc = mag_binned_y_unc(~isnan(mag_binned_y));
            mag_binned_y_unc(mag_binned_y_unc == 0) = 1;
            mag_binned_y = mag_binned_y(~isnan(mag_binned_y));
            dens_binned_y_unc = dens_binned_y_unc(~isnan(dens_binned_y));
            dens_binned_y = dens_binned_y(~isnan(dens_binned_y));

            use_bootstrap = 1;
            % linear fit
            initp = [0, -0.3];
                fixedp = [0, 0];
            if ~use_bootstrap
                [fitp, ~, ffh, std_err, chi_sqr] = fit1D(x_pos, mag_binned, 1 ./ mag_binned_unc.^2, {'line1D'},...
                    initp, fixedp);
            else
                % bootstrap
                [fitp, std_err, fitp_distributions] = bootstrap1d(x_pos, mag_binned,...
                                                        1 ./ mag_binned_unc.^2, {'line1D'},...
                                                        initp, fixedp, [], [],...
                                                        length(initp), 'lsqnonlin', 300, 'fit');
                ffh = @(x) line1D(fitp, x);
            end
           
            mean_sqr_diff = mean((mag_binned - ffh(x_pos)).^2 ./ mag_binned_unc .^2);
            x_interp = linspace(min(x_pos), max(x_pos), 300);
            line_fit_str = sprintf('slope = %0.2e +/- %0.2e', fitp(1), std_err(1));
            
            
            % figure plotting results
            nrows = 4;
            ncols = 6;
            
            fig_handle_mag = figure;
            subplot(nrows, ncols, [1, 2])
            % colormap(df.ColorMap);
            imagesc(mag, [-0.15, 0.15])
            hold on;
            scatter(cx_bin, cy_bin, 'gx');
            axis equal;
            axis image;
            colorbar;
            ax = gca;
            ax.YTickLabel = '';
            ax.XTickLabel = '';
            title('n3 - n1');
            %
            
            subplot(nrows, ncols, [7,8]);
            imagesc(dens, [0, 1]);
            hold on;
            scatter(cx_bin, cy_bin, 'gx');
            axis equal;
            axis image;
            ax = gca;
            ax.YTickLabel = '';
            ax.XTickLabel = '';
            colorbar;
            title('n3 + n1');

            subplot(nrows, ncols, [3,4])
            pol(isnan(pol)) = 0;
            imagesc(pol, [-1, 1]);
            hold on;
            scatter(cx_bin, cy_bin, 'gx');
            axis equal;
            axis image;
            colorbar;
            ax = gca;
            ax.YTickLabel = '';
            ax.XTickLabel = '';
            title('(n3 - n1)/(n3+n1)');

            subplot(nrows, ncols, [9,10])
            imagesc(abs(mag) ./ mag_unc, [0, 5]);
            hold on;
            scatter(cx_bin, cy_bin, 'gx');
            axis equal;
            axis image;
            ax = gca;
            ax.YTickLabel = '';
            ax.XTickLabel = '';
            colorbar;
            title('|n3-n1| / unc');



            subplot(nrows, ncols, [13,14,19,20]);
            errorbar(x_pos, mag_binned, mag_binned_unc, 'o');
            hold on;
            plot(x_interp, ffh(x_interp), 'r');
            max_slope_fitp = fitp;
            max_slope_fitp(1) = max_slope_fitp(1) + std_err(1);
            min_slope_fitp = fitp;
            min_slope_fitp(1) = min_slope_fitp(1) - std_err(1);
            plot(x_interp, line1D(max_slope_fitp, x_interp), 'r--');
            plot(x_interp, line1D(min_slope_fitp, x_interp), 'r--');            
            xlabel('position (sites)');
            ylabel('n3 - n1, perp to grad');
            ylim([-0.1, 0.1]);
            grid on;
            title(line_fit_str);
            
            subplot(nrows, ncols, [15,16,21,22]);
            errorbar(y_pos, mag_binned_y, mag_binned_y_unc, 'o');
            grid on;
            xlabel('position (sites)');
            ylabel('n3 - n1');
            ylim([-0.1, 0.1]);
            title('n3 - n1, along grad');

            subplot(nrows, ncols, [17,18,23,24]);
            errorbar(x_pos, dens_binned, dens_binned_unc, 'o');
            hold on;
            errorbar(y_pos, dens_binned_y, dens_binned_y_unc, 'o');
            grid on;
            ylim([0,1]);
            xlabel('position (site)');
            ylabel('n3 + n1');
            title('density');
            legend({'perp to grad', 'along grad'});
            
            subplot(nrows, ncols, [5,6]);
            mag_cropped(isnan(mag_cropped)) = 0;
            imagesc(mag_cropped, [-0.1, 0.1]);
            hold on;
            scatter(cx_bin, cy_bin, 'gx');
            axis equal;
            axis image;
            colorbar;
            ax = gca;
            ax.YTickLabel = '';
            ax.XTickLabel = '';
            title('cropped region');

            % display bins
            display_bins_x = zeros(size(mask_stack, 1), size(mask_stack,2));

            nmasks = size(mask_stack,3);
            for ii = 1 : nmasks
                curr_mask = mask_stack(:,:,ii);
                curr_mask(curr_mask==1) = ii;
                display_bins_x = display_bins_x + curr_mask;
            end
            
            subplot(nrows, ncols, 11);
            imagesc(display_bins_x .* nan_mask);
%             colorbar;
            ax = gca;
            ax.YTickLabel = '';
            ax.XTickLabel = '';
            title('bins, perp to grad');
            axis equal;
            axis image;

            % display bins y
            display_bins_y = zeros(size(mask_stack_y, 1), size(mask_stack_y,2));

            nmasks = size(mask_stack_y,3);
            for ii = 1 : nmasks
                curr_mask = mask_stack_y(:,:,ii);
                curr_mask(curr_mask==1) = ii;
                display_bins_y = display_bins_y + curr_mask;
            end
            
            subplot(nrows, ncols, 12);
            imagesc(display_bins_y .* nan_mask);
%             colorbar;
            ax = gca;
            ax.YTickLabel = '';
            ax.XTickLabel = '';
            title('bins, along grad');
            axis equal;
            axis image;

            n3 = mean(obj.threes.AtomNumbers);
            n1 = mean(obj.ones.AtomNumbers);
            p_global = (n3 - n1) / (n3 + n1);
            suptitle(sprintf('%s\n angle = %0.0f deg, mean dens = %0.2f with std = %0.2f, global polarization = %0.2f',...
                obj.identifier, angle_grad * 180/pi, dens_mean, dens_std, p_global)); 

            obj.ones.getColorMap;
            colormap(obj.ones.ColorMap);

        end
        
        function fig_handle = get_mag_angle(obj, crop_radius, angles, bin_width, bin_2d_size)
            fitp = zeros(length(angles), 2);
            std_err = zeros(length(angles), 2);
            mean_sqr_diff = zeros(length(angles), 1);
            
            for ii = 1:length(angles)
                [fitp(ii, :), std_err(ii, :), mean_sqr_diff(ii), ~] = obj.get_magnetization(crop_radius, angles(ii), bin_width, bin_2d_size);
            end
            
            fig_handle = figure;
            subplot(1, 3, 1);
            errorbar(angles * 180/pi, fitp(:, 1), std_err(:, 1));
            xlabel('angle (deg)');
            ylabel('slope (mag/site)');
            grid on;
            
            subplot(1, 3, 2);
            errorbar(angles * 180/pi, fitp(:, 2), std_err(:, 2));
            xlabel('angle (deg)');
            ylabel('offset');
            grid on;
            
            subplot(1, 3, 3);
            plot(angles * 180/pi, mean_sqr_diff);
            xlabel('angle (deg)');
            ylabel('mean square diff');
            grid on;
        end
        
        function [fit_params_dqmc, std_errs_dqmc, redchi_sqr, fig_handle_dqmc] = ...
            fit_dqmc(obj, init_params, fixed_params, use_to_fit, density_limits)

            if ~exist('init_params', 'var') || isempty(init_params)
                init_params = [-4.5, 0.5, 0, 0.97, 0.9];
            end
            
            if ~exist('fixed_params', 'var') || isempty(fixed_params)
                fixed_params = [0, 0, 1, 1, 1];
            end

            if ~exist('use_to_fit', 'var') || isempty(use_to_fit)
                use_to_fit = [];
            end
            
            if ~exist('density_limits', 'var') || isempty(density_limits)
                density_limits = [0.1, 2];
            end
                      
            dqmc = load(obj.dqmc_data_path);
            [fit_params_dqmc, std_errs_dqmc, redchi_sqr, mu_struct, fit_struct] = ...
                fit_dqmc_attractive(obj.ones, obj.threes, obj.singles, obj.doubles,...
                                    use_to_fit, density_limits, init_params, fixed_params, dqmc);
			fig_handle_dqmc = display_dqmc_fit(fit_struct, mu_struct);

            % get suptitle and add our identifier string to it
            tags = cell(length(fig_handle_dqmc.Children), 1);
            for ii = 1:length(fig_handle_dqmc.Children)
                tags{ii} = fig_handle_dqmc.Children(ii).Tag;
            end
            fn = @(c) strcmp(c, 'suptitle');
            index = find(cellfun(fn, tags));
            fig_handle_dqmc.Children(index).Children.String = ...
                vertcat(obj.identifier, fig_handle_dqmc.Children(index).Children.String);
            %sprintf('%s\n%s',...
             %   obj.identifier, char(fig_handle_dqmc.Children(index).Children.String));
        end
        
        function [fig_handle, group_identifier] = showSets(obj, instance_stack)
            if ~exist('instance_stack', 'var')
                instance_stack = [];
            end
            
            instance_stack = vertcat(obj, instance_stack);
            
            fig_handle = figure;
            leg_singles = {};
            leg_doubles = {};
            leg_doubles_corr = {};
            leg_dens = {};
            leg_corr = {};
            leg_component_dens = {};
            
            nrows = 2;
            ncols = 3;
            for ii = 1:length(instance_stack)
                inst = instance_stack(ii);
                
                interp_ns = linspace(0, max(inst.dens13_azavg));
                
                % singles versus density
                subplot(nrows, ncols, 1);
                errorbar(inst.dens13_azavg, inst.singles.Occs_AzAvg,...
                         inst.singles.Occs_AzAvgUnc, inst.singles.Occs_AzAvgUnc,...
                         inst.dens13_azavg_unc, inst.dens13_azavg_unc);
                hold on;
                plot(interp_ns, interp_ns - 2*(interp_ns/2).^2, 'g');
                hold on;
                xlabel('<n>');
                ylabel('<n^s>');
                grid on;
                title('singles vs density');
                leg_singles = horzcat(leg_singles, sprintf('%03d', inst.singles.Dataset{4}));

                % doubles vs density
                subplot(nrows, ncols, 4);
                
                errorbar(inst.dens13_azavg, inst.doubles.Occs_AzAvg,...
                         inst.doubles.Occs_AzAvgUnc, inst.doubles.Occs_AzAvgUnc,...
                         inst.dens13_azavg_unc, inst.dens13_azavg_unc);
               hold on;   
               
               % determine number of doubles expected based on density of
               % 1s and 3s
                expected_doubles = (inst.dens13_azavg - inst.singles.Occs_AzAvg) / 2;
                expected_doubles_unc = 0.5 * sqrt(inst.dens13_azavg_unc .^2 + inst.singles.Occs_AzAvgUnc .^2);

               % infer doublon imaging efficiency
                doublon_efficiency = inst.doubles.Occs_AzAvg ./ expected_doubles;
                doublon_efficiency_unc = doublon_efficiency .* ...
                    sqrt( (inst.doubles.Occs_AzAvgUnc ./ inst.doubles.Occs_AzAvg).^2 + ...
                          (expected_doubles_unc ./ expected_doubles).^2 );
                
                
                errorbar(inst.dens13_azavg, expected_doubles, expected_doubles_unc);
                 
                plot(interp_ns, (interp_ns/2).^2, 'g');
                errorbar(inst.dens13_azavg, doublon_efficiency, doublon_efficiency_unc);
                
                xlabel('<n>');
                ylabel('<d>');
                grid on;
                title('doubles vs density');
                leg_doubles = horzcat(leg_doubles, sprintf('%03d', inst.doubles.Dataset{4}));   
                leg_doubles = horzcat(leg_doubles, 'expected doubles'); 
                
                ylim([-0.1, 1]);
                
                subplot(nrows, ncols, 2);
                corr1 = squeeze(inst.ones.Density_Corr_AzAvg(inst.ones.CenterIndex_CorrMatrix, inst.ones.CenterIndex_CorrMatrix + 1, :));
                corr1_unc = squeeze(inst.ones.Density_Corr_AzAvgUnc(inst.ones.CenterIndex_CorrMatrix, inst.ones.CenterIndex_CorrMatrix + 1, :));
               
                corr3 = squeeze(inst.threes.Density_Corr_AzAvg(inst.threes.CenterIndex_CorrMatrix, inst.threes.CenterIndex_CorrMatrix + 1, :));
                corr3_unc = squeeze(inst.threes.Density_Corr_AzAvgUnc(inst.threes.CenterIndex_CorrMatrix, inst.threes.CenterIndex_CorrMatrix + 1, :));
               
                errorbar(inst.dens13_azavg, corr1, corr1_unc, corr1_unc,...
                    inst.dens13_azavg_unc, inst.dens13_azavg_unc);
                hold on;
                errorbar(inst.dens13_azavg, corr3, corr3_unc, corr3_unc,...
                    inst.dens13_azavg_unc, inst.dens13_azavg_unc);
               
                xlabel('<n>');
                ylabel('<nn>');
                grid on;
                title('density corr vs density');
                leg_corr = horzcat(leg_corr, sprintf('<11> %03d', inst.ones.Dataset{4}), sprintf('<33> %03d', inst.threes.Dataset{4}));

                subplot(nrows, ncols, 5);
                dd = squeeze(inst.doubles.Density_Corr_AzAvg(inst.doubles.CenterIndex_CorrMatrix, inst.doubles.CenterIndex_CorrMatrix + 1, :));
                dd_unc = squeeze(inst.doubles.Density_Corr_AzAvgUnc(inst.doubles.CenterIndex_CorrMatrix, inst.doubles.CenterIndex_CorrMatrix + 1, :));
                
                errorbar(inst.dens13_azavg, dd, dd_unc, dd_unc,...
                    inst.dens13_azavg_unc, inst.dens13_azavg_unc);
                hold on;
                xlabel('<n>');
                ylabel('<dd>');
                grid on;
                title('doublon corr vs density');
                leg_doubles_corr = horzcat(leg_doubles_corr, sprintf('%03d', inst.doubles.Dataset{4}));

                subplot(nrows, ncols, 3);
                errorbar(inst.ones.RadialPos, inst.dens13_azavg,...
                    inst.dens13_azavg_unc, inst.dens13_azavg_unc,...
                    inst.ones.RadialPosUnc, inst.ones.RadialPosUnc);
                hold on;
                errorbar(inst.ones.RadialPos, inst.denssd_azavg,...
                    inst.denssd_azavg_unc, inst.denssd_azavg_unc,...
                    inst.ones.RadialPosUnc, inst.ones.RadialPosUnc);
                errorbar(inst.dens_azavg_dist, inst.dens_azavg, inst.dens_azavg_unc);

                xlabel('position (sites)');
                ylabel('<n>');
                grid on;
                title('total dens vs position');
                leg_dens = horzcat(leg_dens, sprintf('%03d, %03d, 13', inst.ones.Dataset{4}, inst.threes.Dataset{4}));
                leg_dens = horzcat(leg_dens, sprintf('%03d, %03d, sd', inst.singles.Dataset{4}, inst.doubles.Dataset{4}));
                leg_dens = horzcat(leg_dens, sprintf('%03d, %03d, 13 then az avg', inst.ones.Dataset{4}, inst.threes.Dataset{4}));
                
                subplot(nrows, ncols, 6);
                errorbar(inst.ones.RadialPos, inst.ones.Occs_AzAvg,...
                    inst.ones.Occs_AzAvgUnc, inst.ones.Occs_AzAvgUnc,...
                    inst.ones.RadialPosUnc, inst.ones.RadialPosUnc);
                hold on;
                errorbar(inst.threes.RadialPos, inst.threes.Occs_AzAvg,...
                    inst.threes.Occs_AzAvgUnc, inst.threes.Occs_AzAvgUnc,...
                    inst.threes.RadialPosUnc, inst.threes.RadialPosUnc);
                xlabel('position (sites)');
                ylabel('<n_\sigma>');
                grid on;
                title('component dens vs position');
                leg_component_dens = horzcat(leg_component_dens,...
                                     sprintf('%03d 1s', inst.ones.Dataset{4}));
                leg_component_dens = horzcat(leg_component_dens,...
                                     sprintf('%03d 3s', inst.threes.Dataset{4}));
            
            end
            subplot(nrows, ncols, 1);
            legend(leg_singles);
            subplot(nrows, ncols, 4);
            legend(leg_doubles);
            subplot(nrows, ncols, 2);
            legend(leg_corr);
            subplot(nrows, ncols, 5);
            legend(leg_doubles_corr);
            subplot(nrows, ncols, 3);
            legend(leg_dens);
            subplot(nrows, ncols, 6);
            legend(leg_component_dens);
           
            
            group_identifier = sprintf('%04d-%02d-%02d_folders=%03d-%03d',...
                instance_stack(1).date{1}, instance_stack(1).date{2}, instance_stack(1).date{3},...
                min([instance_stack.folders]), max([instance_stack.folders]));
            suptitle(group_identifier); 
        end
    
        function [fig_handle] = showFlatness(obj)
            fig_handle = figure;
            nrows = 2;
            ncols = 2;
            
            subplot(nrows, ncols, 1);
            imagesc(obj.dens, [0, 1]);
            axis equal;
            axis image;
            title('Density');
            
            subplot(nrows, ncols, 3);
            errorbar(obj.dens_azavg_dist, obj.dens_azavg, obj.dens_azavg_unc);
            grid on;
            xlabel('Position');
            ylabel('<n>');
            
            subplot(nrows, ncols, 2);
            imagesc(obj.ones.Occs_ImgAvg, [0, 0.5]);
            axis equal;
            axis image;
            title('Ones');
            
            subplot(nrows, ncols, 4);
            imagesc(obj.threes.Occs_ImgAvg, [0, 0.5]);
            axis equal;
            axis image;
            title('Threes');
            
            suptitle(sprintf('%s', obj.identifier));
            
        end
        
    end
end