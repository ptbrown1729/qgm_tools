function [fit_params] = min_img_xform(fn, init_pt)
    % simple gradient descent solver to try and understand what is going
    % wrong in fitting image transformation.

    iters = 1000;
    
    convergence_tol = 1e-8;
    
    % parameters to store values
    fn_vals = zeros(1, iters);
    pts = zeros(iters, length(init_pt));
    directions = zeros(iters, length(init_pt));
    step_sizes = zeros(1, iters);
    jacobians = zeros(iters, length(init_pt));
    jacobians_norm = zeros(1, iters);
    angles = zeros(1, iters);
    
    % initial step size and start point
    step_sizes(1) = 0.1;
    pts(1, :) = init_pt;
    
    for ii = 1 : iters
    
        % calculate jacobian
        [jacobians(ii, :), fn_vals(ii)] = get_jacobian(fn, pts(ii, :),...
                                                        step_sizes(ii) * ones(1, length(init_pt)));
        % store information about jacobian
        jacobians_norm(ii) = sqrt(sum(jacobians(ii, :).^2));
        directions(ii, :) = jacobians(ii, :) / jacobians_norm(ii);
        
        % gradient descent
        pts(ii + 1, :) = pts(ii, :) - step_sizes(ii) * directions(ii, :);
        
        % calculate angle between last two directions
        if ii > 1
            angles(ii) = sum( jacobians(ii - 1, :) .* jacobians(ii, :) )...
                        / jacobians_norm(ii) / jacobians_norm(ii-1);
        end
        
        % increment step size
        step_sizes(ii + 1) = step_sizes(ii);
        % if oscillating, decrease step size
        if ii > 2
            if all(angles(ii-2:ii) < -0.2)
                step_sizes(ii + 1) = step_sizes(ii) / 2;
            elseif ii > 5 && all(angles(ii-5:ii) > 0.8)
                step_sizes(ii + 1) = step_sizes(ii) * 2;
            end
        end
        
        % end calculation if too many are the same
        if ii > 10 && all( abs(fn_vals(ii - 10 : ii) - fn_vals(ii)) < convergence_tol)
            disp('ending early because converged');
            break;
        end
        
    end
    
    last_iteration = ii;
    
%     [~, I] = min(fn_vals(1:last_iteration));
%     fit_params = pts(I, :);
    fit_params = pts(last_iteration, :);
end