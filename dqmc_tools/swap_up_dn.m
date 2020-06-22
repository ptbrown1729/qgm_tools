function [dqmc_results_stack_out]  = swap_up_dn(dqmc_results_stack_in)
% Given a dqmc structure produced by convert_qmc_file(), if the structure
% has spin imbalance (mu_up != mu_dn), then flip produce a new structure
% with up and down spins switched.
%
% dqmc_results_stack_in is a vertical array of dqmc structures.
%
% dqmc_results_stack_out is dqmc_results_stack_in with these additional
% structures appended.

dqmc_results_stack_out = dqmc_results_stack_in;

for ii = 1:length(dqmc_results_stack_in)
    
    
    if (dqmc_results_stack_in(ii).mu_dn ~= dqmc_results_stack_in(ii).mu_up) &&...
       (dqmc_results_stack_in(ii).t_dn  == dqmc_results_stack_in(ii).t_up) 
   
        fprintf('flipping spin of %d/%d\n', ii, length(dqmc_results_stack_in));
   
        dqmc_result_flipped_spins = dqmc_results_stack_in(ii);
   
        dqmc_result_flipped_spins.mu_up = dqmc_results_stack_in(ii).mu_dn;
        dqmc_result_flipped_spins.mu_dn = dqmc_results_stack_in(ii).mu_up;
        
        dqmc_results_flipped_spins.up_equal_time_greens_function = ...
                dqmc_results_stack_in(ii).down_equal_time_greens_function;
        dqmc_results_flipped_spins.down_equal_time_greens_function = ...
            dqmc_results_stack_in(ii).up_equal_time_greens_function;
        % not correctly handling density_density_correlation_fn_up_up
        
        % equal time
        dqmc_result_flipped_spins.equal_time_measurements.up_spin_occupancy = ...
            dqmc_results_stack_in(ii).equal_time_measurements.down_spin_occupancy;
        dqmc_result_flipped_spins.equal_time_measurements.down_spin_occupancy = ...
            dqmc_results_stack_in(ii).equal_time_measurements.up_spin_occupancy;
        % signs of equal time
        dqmc_result_flipped_spins.sign_of_equal_time_measurements.avg_up_sign = ...
            dqmc_results_stack_in(ii).sign_of_equal_time_measurements.avg_dn_sign;
        dqmc_result_flipped_spins.sign_of_equal_time_measurements.avg_dn_sign = ...
            dqmc_results_stack_in(ii).sign_of_equal_time_measurements.avg_up_sign;
            
        dqmc_results_stack_out = vertcat(dqmc_results_stack_out, dqmc_result_flipped_spins);
    else
    end
end

end